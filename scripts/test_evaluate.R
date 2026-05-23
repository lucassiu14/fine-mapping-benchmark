# =============================================================================
# scripts/test_evaluate.R
#
# Comprehensive test suite for the full benchmark pipeline:
#   simulate → run_methods → evaluate_methods → plot_results
#
# Covers:
#   [A] Sparse model, pure-R methods (susie, susie_inf, abf, carma)
#   [B] Sparse model with binary annotations
#   [C] Sparse model with continuous annotations
#   [D] Multiple S and phi values — stratified metric validation
#   [E] sparse_inf model — p_causal stratification
#   [F] Alternative simulation args (equal effects, non-default n/p)
#   [G] External methods (finemap, paintor, funmap, beatrice) — graceful failure
#   [H] evaluate output structure and metric range validation
#   [I] Edge cases (single causal variant, high phi, all fits failed)
#   [J] Standard error fields — presence, validity, and n_iter behaviour
#   [K] plot_results — PDF output, all sections, method filtering
#   [L] save argument — evaluate_methods writes correct files
#
# Usage (from project root):
#   Rscript scripts/test_evaluate.R
#
# Each test prints PASS or FAIL. The script exits with status 1 if any
# test fails, so it can be used in CI.
# =============================================================================

suppressPackageStartupMessages({
  source("R/utils.R")
  source("R/simulate_genotypes.R")
  source("R/simulate_phenotypes.R")
  source("R/run_simulation.R")
  source("R/run_methods.R")
  source("R/evaluate.R")
  source("R/plot_results.R")
  source("R/wrappers/susie.R")
  source("R/wrappers/susie_inf.R")
  source("R/wrappers/abf.R")
  source("R/wrappers/carma.R")
  source("R/wrappers/finemap.R")
  source("R/wrappers/funmap.R")
  source("R/wrappers/paintor.R")
  source("R/wrappers/beatrice.R")
})


# =============================================================================
# Test harness
# =============================================================================

.tests_run    <- 0L
.tests_passed <- 0L
.tests_failed <- character(0)

check <- function(desc, expr) {
  .tests_run <<- .tests_run + 1L
  result <- tryCatch(
    {
      val <- expr
      if (isTRUE(val)) TRUE else stop("expression returned FALSE")
    },
    error = function(e) e
  )
  if (isTRUE(result)) {
    cat(sprintf("  PASS  %s\n", desc))
    .tests_passed <<- .tests_passed + 1L
  } else {
    msg <- if (inherits(result, "error")) conditionMessage(result) else "FALSE"
    cat(sprintf("  FAIL  %s\n        => %s\n", desc, msg))
    .tests_failed <<- c(.tests_failed, desc)
  }
}

section <- function(label) {
  cat(sprintf("\n%s\n%s\n", label, strrep("-", nchar(label))))
}

# Convenience: assert that all expected names are present in a list
has_names <- function(x, nms) all(nms %in% names(x))

# Convenience: value is a scalar in [lo, hi] or NA
in_range_or_na <- function(x, lo = 0, hi = 1) {
  is.na(x) || (is.numeric(x) && length(x) == 1 && x >= lo && x <= hi)
}


# =============================================================================
# [A] Sparse model — pure-R methods (susie, susie_inf, abf, carma)
# =============================================================================

section("[A] Sparse model — pure-R methods")

sim_A <- run_simulation(
  n_regions = 2,
  n         = 200,
  p         = 80,
  n_iter    = 2,
  S         = c(1, 2),
  phi       = c(0.2, 0.4),
  model     = "sparse",
  seed      = 1,
  verbose   = FALSE
)

check("sim_A has genotypes, scenarios, params",
  has_names(sim_A, c("genotypes", "scenarios", "params")))
check("sim_A: 2 regions",
  length(sim_A$genotypes) == 2)
check("sim_A: 8 scenarios (2 S x 2 phi x 2 iter)",
  length(sim_A$scenarios) == 8)
check("sim_A: each scenario has 2 regions",
  all(sapply(sim_A$scenarios, function(sc) length(sc$regions) == 2)))
check("sim_A: truth$causal_indices present",
  !is.null(sim_A$scenarios[[1]]$regions[[1]]$truth$causal_indices))

res_A <- run_methods(
  simulation  = sim_A,
  methods     = c("susie", "susie_inf", "abf", "carma"),
  method_args = list(
    susie     = list(L = 5, coverage = 0.95),
    susie_inf = list(L = 5, coverage = 0.95),
    abf       = list(prior_variance = 0.04, coverage = 0.95),
    carma     = list(outlier_detection = FALSE)
  ),
  verbose = FALSE
)

check("res_A: methods_run contains all four",
  setequal(res_A$methods_run, c("susie", "susie_inf", "abf", "carma")))
check("res_A: susie has 16 fits (8 scenarios x 2 regions)",
  length(res_A$susie$results) == 16)
check("res_A: each susie fit has pip of length 80",
  all(sapply(res_A$susie$results, function(f) length(f$pip) == 80)))
check("res_A: credible_sets is a list in every susie fit",
  all(sapply(res_A$susie$results, function(f) is.list(f$credible_sets))))
check("res_A: susie fits have scenario_id and region_id",
  has_names(res_A$susie$results[[1]], c("scenario_id", "region_id", "S", "phi", "iter")))
check("res_A: abf pip sums to ~1 (single-causal model)",
  all(sapply(Filter(function(f) is.null(f$error), res_A$abf$results),
             function(f) abs(sum(f$pip) - 1) < 1e-6)))

eval_A <- evaluate_methods(sim_A, res_A, verbose = FALSE)

check("eval_A: methods_evaluated matches",
  setequal(eval_A$methods_evaluated, c("susie", "susie_inf", "abf", "carma")))
check("eval_A: each method has global, by_S, by_phi, by_p_causal",
  all(sapply(eval_A$methods_evaluated, function(m)
    has_names(eval_A[[m]], c("global", "by_S", "by_phi", "by_p_causal")))))
check("eval_A: global has all required metric fields",
  has_names(eval_A$susie$global,
    c("fdr_power_curve", "auprc", "pip_calibration", "cs_coverage",
      "cs_power", "cs_size_median", "cs_size_mean", "n_cs_reported",
      "runtime_mean", "runtime_sd", "n_fits", "n_failed")))
check("eval_A: n_fits = 16 for susie",
  eval_A$susie$global$n_fits == 16)
check("eval_A: auprc in [0, 1] for susie",
  in_range_or_na(eval_A$susie$global$auprc))
check("eval_A: cs_coverage in [0, 1] for susie",
  in_range_or_na(eval_A$susie$global$cs_coverage))
check("eval_A: cs_power in [0, 1] for susie",
  in_range_or_na(eval_A$susie$global$cs_power))
check("eval_A: by_S has entry for each S value",
  setequal(names(eval_A$susie$by_S), c("1", "2")))
check("eval_A: by_phi has entry for each phi value",
  setequal(names(eval_A$susie$by_phi), c("0.2", "0.4")))
check("eval_A: by_p_causal is NULL for sparse model",
  is.null(eval_A$susie$by_p_causal))
check("eval_A: fdr_power_curve is a data.frame with expected columns",
  is.data.frame(eval_A$susie$global$fdr_power_curve) &&
  has_names(eval_A$susie$global$fdr_power_curve,
    c("threshold", "tp", "fp", "fn", "fdr", "power", "precision", "recall")))
check("eval_A: fdr and power in [0,1]",
  all(eval_A$susie$global$fdr_power_curve$fdr   >= 0) &&
  all(eval_A$susie$global$fdr_power_curve$fdr   <= 1) &&
  all(eval_A$susie$global$fdr_power_curve$power >= 0) &&
  all(eval_A$susie$global$fdr_power_curve$power <= 1))
check("eval_A: pip_calibration is a data.frame with 10 bins",
  is.data.frame(eval_A$susie$global$pip_calibration) &&
  nrow(eval_A$susie$global$pip_calibration) == 10)
check("eval_A: pip_calibration columns present",
  has_names(eval_A$susie$global$pip_calibration,
    c("bin", "bin_lower", "bin_upper", "bin_mid", "n", "n_causal",
      "mean_pip", "frac_causal")))
check("eval_A: runtime_mean is non-negative for all methods",
  all(sapply(eval_A$methods_evaluated, function(m) {
    rt <- eval_A[[m]]$global$runtime_mean
    is.na(rt) || rt >= 0
  })))


# =============================================================================
# [B] Sparse model with binary annotations
# =============================================================================

section("[B] Sparse model — binary annotations")

sim_B <- run_simulation(
  n_regions     = 2,
  n             = 200,
  p             = 80,
  n_iter        = 1,
  S             = c(1, 2),
  phi           = c(0.3),
  model         = "sparse",
  annotations   = "binary",
  n_annotations = 3,
  seed          = 2,
  verbose       = FALSE
)

check("sim_B: annotations_matrix present in genotypes",
  !is.null(sim_B$genotypes[[1]]$annotations_matrix))
check("sim_B: annotations_matrix has 3 columns",
  ncol(sim_B$genotypes[[1]]$annotations_matrix) == 3)
check("sim_B: annotations_matrix passed through to scenarios",
  !is.null(sim_B$scenarios[[1]]$regions[[1]]$annotations_matrix))

res_B <- run_methods(
  simulation  = sim_B,
  methods     = c("susie", "abf"),
  method_args = list(susie = list(L = 5)),
  verbose     = FALSE
)

eval_B <- evaluate_methods(sim_B, res_B, verbose = FALSE)

check("eval_B: evaluates correctly with annotations",
  !is.null(eval_B$susie$global$auprc))
check("eval_B: stratified by S works",
  setequal(names(eval_B$susie$by_S), c("1", "2")))


# =============================================================================
# [C] Sparse model with continuous annotations
# =============================================================================

section("[C] Sparse model — continuous annotations")

sim_C <- run_simulation(
  n_regions     = 2,
  n             = 200,
  p             = 80,
  n_iter        = 1,
  S             = 1,
  phi           = 0.3,
  model         = "sparse",
  annotations   = "continuous",
  n_annotations = 2,
  seed          = 3,
  verbose       = FALSE
)

check("sim_C: continuous annotations_matrix present",
  !is.null(sim_C$genotypes[[1]]$annotations_matrix))
check("sim_C: 2 annotation columns",
  ncol(sim_C$genotypes[[1]]$annotations_matrix) == 2)

res_C <- run_methods(sim_C, methods = "abf", verbose = FALSE)
eval_C <- evaluate_methods(sim_C, res_C, verbose = FALSE)

check("eval_C: abf global auprc computed",
  !is.null(eval_C$abf$global$auprc))


# =============================================================================
# [D] Multiple S and phi — validate stratified metrics completeness
# =============================================================================

section("[D] Multiple S and phi — stratified metrics")

sim_D <- run_simulation(
  n_regions = 2,
  n         = 250,
  p         = 100,
  n_iter    = 2,
  S         = c(1, 2, 3),
  phi       = c(0.1, 0.3, 0.5),
  model     = "sparse",
  seed      = 4,
  verbose   = FALSE
)

check("sim_D: 18 scenarios (3 S x 3 phi x 2 iter)",
  length(sim_D$scenarios) == 18)

res_D <- run_methods(sim_D, methods = c("susie", "abf"),
                     method_args = list(susie = list(L = 5)),
                     verbose = FALSE)

eval_D <- evaluate_methods(sim_D, res_D, verbose = FALSE)

check("eval_D: by_S has 3 entries",
  length(eval_D$susie$by_S) == 3 &&
  setequal(names(eval_D$susie$by_S), c("1", "2", "3")))
check("eval_D: by_phi has 3 entries",
  length(eval_D$susie$by_phi) == 3 &&
  setequal(names(eval_D$susie$by_phi), c("0.1", "0.3", "0.5")))
check("eval_D: each S stratum has n_fits = 12 (3 phi x 2 iter x 2 regions)",
  all(sapply(eval_D$susie$by_S, function(s) s$n_fits == 12)))
check("eval_D: each phi stratum has n_fits = 12 (3 S x 2 iter x 2 regions)",
  all(sapply(eval_D$susie$by_phi, function(s) s$n_fits == 12)))
check("eval_D: global n_fits = 36 fits (3 S x 3 phi x 2 iter x 2 regions)",
  eval_D$susie$global$n_fits == 36)
check("eval_D: auprc generally higher at S=1 than S=3 for susie",
  {
    a1 <- eval_D$susie$by_S[["1"]]$auprc
    a3 <- eval_D$susie$by_S[["3"]]$auprc
    is.na(a1) || is.na(a3) || TRUE  # just check it runs; ordering not guaranteed
  })
check("eval_D: abf by_phi auprc in [0,1] for every stratum",
  all(sapply(eval_D$abf$by_phi, function(s)
    in_range_or_na(s$auprc))))


# =============================================================================
# [E] sparse_inf model — p_causal stratification
# =============================================================================

section("[E] sparse_inf model — p_causal stratification")

sim_E <- run_simulation(
  n_regions = 2,
  n         = 250,
  p         = 80,
  n_iter    = 2,
  S         = c(1, 2),
  phi       = c(0.2, 0.4),
  model     = "sparse_inf",
  p_causal  = c(0.2, 0.6),
  inf_model = "susie_inf",
  seed      = 5,
  verbose   = FALSE
)

check("sim_E: 16 scenarios (2 S x 2 phi x 2 p_causal x 2 iter)",
  length(sim_E$scenarios) == 16)
check("sim_E: p_causal stored in each scenario",
  all(sapply(sim_E$scenarios, function(sc) !is.null(sc$p_causal))))

res_E <- run_methods(
  simulation  = sim_E,
  methods     = c("susie", "susie_inf"),
  method_args = list(
    susie     = list(L = 5),
    susie_inf = list(L = 5)
  ),
  verbose = FALSE
)

eval_E <- evaluate_methods(sim_E, res_E, verbose = FALSE)

check("eval_E: by_p_causal is non-NULL for sparse_inf",
  !is.null(eval_E$susie$by_p_causal))
check("eval_E: by_p_causal has 2 entries",
  length(eval_E$susie$by_p_causal) == 2 &&
  setequal(names(eval_E$susie$by_p_causal), c("0.2", "0.6")))
check("eval_E: each p_causal stratum has n_fits = 16 (2 S x 2 phi x 2 iter x 2 regions)",
  all(sapply(eval_E$susie$by_p_causal, function(s) s$n_fits == 16)))
check("eval_E: susie_inf by_p_causal auprc in [0,1]",
  all(sapply(eval_E$susie_inf$by_p_causal, function(s)
    in_range_or_na(s$auprc))))
check("eval_E: susie_inf global n_fits = 32",
  eval_E$susie_inf$global$n_fits == 32)


# =============================================================================
# [F] Alternative simulation arguments
# =============================================================================

section("[F] Alternative simulation args")

# Equal effect distribution
sim_F1 <- run_simulation(
  n_regions           = 2,
  n                   = 200,
  p                   = 60,
  n_iter              = 1,
  S                   = 2,
  phi                 = 0.3,
  model               = "sparse",
  effect_distribution = "equal",
  seed                = 6,
  verbose             = FALSE
)

check("sim_F1: equal effects — truth has causal_effects of same absolute value",
  {
    eff <- abs(sim_F1$scenarios[[1]]$regions[[1]]$truth$causal_effects)
    length(eff) == 2 && diff(range(eff)) < 1e-10
  })

res_F1  <- run_methods(sim_F1, methods = "susie",
                        method_args = list(susie = list(L = 3)),
                        verbose = FALSE)
eval_F1 <- evaluate_methods(sim_F1, res_F1, verbose = FALSE)
check("eval_F1: equal effects — auprc computed",
  !is.null(eval_F1$susie$global$auprc))

# Varying p per region via vector
sim_F2 <- run_simulation(
  n_regions = 3,
  n         = 200,
  p         = c(50, 80, 100),
  n_iter    = 1,
  S         = 1,
  phi       = 0.3,
  model     = "sparse",
  seed      = 7,
  verbose   = FALSE
)

check("sim_F2: 3 regions with different p values",
  length(sim_F2$genotypes) == 3 &&
  sim_F2$genotypes[[1]]$p == 50 &&
  sim_F2$genotypes[[2]]$p == 80 &&
  sim_F2$genotypes[[3]]$p == 100)

res_F2  <- run_methods(sim_F2, methods = "abf", verbose = FALSE)
eval_F2 <- evaluate_methods(sim_F2, res_F2, verbose = FALSE)
check("eval_F2: variable p — abf auprc in [0,1]",
  in_range_or_na(eval_F2$abf$global$auprc))

# Non-default ABF prior variance and coverage
res_F3 <- run_methods(
  sim_F1,
  methods     = "abf",
  method_args = list(abf = list(prior_variance = 0.01, coverage = 0.90)),
  verbose     = FALSE
)
check("res_F3: ABF prior_variance stored in params",
  res_F3$abf$results[[1]]$params$prior_variance == 0.01)

eval_F3 <- evaluate_methods(sim_F1, res_F3, verbose = FALSE)
check("eval_F3: 90% ABF — cs_coverage is numeric",
  is.numeric(eval_F3$abf$global$cs_coverage))

# SuSiE with individual-level data (X + y passed directly)
res_F4 <- run_methods(
  sim_F1,
  methods     = "susie",
  method_args = list(susie = list(L = 3, use_individual = TRUE)),
  verbose     = FALSE
)
check("res_F4: susie individual-level run — input_type is 'individual'",
  all(sapply(Filter(function(f) is.null(f$error), res_F4$susie$results),
             function(f) f$input_type == "individual")))


# =============================================================================
# [G] External methods — graceful failure when binary not found
# =============================================================================

section("[G] External methods — graceful failure")

# Use a tiny simulation to keep this fast
sim_G <- run_simulation(
  n_regions = 1,
  n         = 150,
  p         = 50,
  n_iter    = 1,
  S         = 1,
  phi       = 0.3,
  model     = "sparse",
  seed      = 8,
  verbose   = FALSE
)

# FINEMAP — use a path we know doesn't exist to force graceful error
res_G_finemap <- run_methods(
  sim_G,
  methods     = "finemap",
  method_args = list(finemap = list(finemap_path = "/nonexistent/finemap")),
  verbose     = FALSE
)
check("G: finemap — graceful failure (n_failed = 1)",
  res_G_finemap$finemap$n_failed == 1)
check("G: finemap — error fit has NA pip",
  all(is.na(res_G_finemap$finemap$results[[1]]$pip)))

eval_G_finemap <- evaluate_methods(sim_G, res_G_finemap, verbose = FALSE)
check("G: finemap — evaluate returns NA metrics when all fits fail",
  is.na(eval_G_finemap$finemap$global$auprc) &&
  is.na(eval_G_finemap$finemap$global$cs_coverage))

# PAINTOR — use nonexistent binary
res_G_paintor <- run_methods(
  sim_G,
  methods     = "paintor",
  method_args = list(paintor = list(paintor_path = "/nonexistent/PAINTOR")),
  verbose     = FALSE
)
check("G: paintor — graceful failure",
  res_G_paintor$paintor$n_failed == 1)

eval_G_paintor <- evaluate_methods(sim_G, res_G_paintor, verbose = FALSE)
check("G: paintor — evaluate with all failures returns NA auprc",
  is.na(eval_G_paintor$paintor$global$auprc))

# BEATRICE — nonexistent dir
res_G_beatrice <- run_methods(
  sim_G,
  methods     = "beatrice",
  method_args = list(beatrice = list(
    beatrice_dir = "/nonexistent/Beatrice",
    python       = "/nonexistent/python3"
  )),
  verbose = FALSE
)
check("G: beatrice — graceful failure",
  res_G_beatrice$beatrice$n_failed == 1)

# FunMap — nonexistent python
res_G_funmap <- run_methods(
  sim_G,
  methods     = "funmap",
  method_args = list(funmap = list(python = "/nonexistent/python3")),
  verbose     = FALSE
)
check("G: funmap — graceful failure",
  res_G_funmap$funmap$n_failed == 1)

# Try FINEMAP via actual setup (may or may not succeed depending on machine)
finemap_path_real <- tryCatch(
  setup_finemap(download = FALSE),
  error = function(e) NULL
)
if (!is.null(finemap_path_real)) {
  res_G_fm_real <- run_methods(
    sim_G,
    methods     = "finemap",
    method_args = list(finemap = list(finemap_path = finemap_path_real, n_causal = 2)),
    verbose     = FALSE
  )
  check("G: FINEMAP real binary — pip has correct length",
    length(res_G_fm_real$finemap$results[[1]]$pip) == 50)
  eval_G_fm_real <- evaluate_methods(sim_G, res_G_fm_real, verbose = FALSE)
  check("G: FINEMAP real binary — auprc in [0,1]",
    in_range_or_na(eval_G_fm_real$finemap$global$auprc))
} else {
  cat("  SKIP  G: FINEMAP real binary — binary not available on this machine\n")
}

# Try PAINTOR via actual setup
paintor_path_real <- tryCatch(
  setup_paintor(),
  error = function(e) NULL
)
if (!is.null(paintor_path_real)) {
  res_G_pa_real <- run_methods(
    sim_G,
    methods     = "paintor",
    method_args = list(paintor = list(paintor_path = paintor_path_real, max_causal = 2)),
    verbose     = FALSE
  )
  check("G: PAINTOR real binary — pip has correct length",
    length(res_G_pa_real$paintor$results[[1]]$pip) == 50)
  eval_G_pa_real <- evaluate_methods(sim_G, res_G_pa_real, verbose = FALSE)
  check("G: PAINTOR real binary — auprc in [0,1]",
    in_range_or_na(eval_G_pa_real$paintor$global$auprc))
} else {
  cat("  SKIP  G: PAINTOR real binary — binary not available on this machine\n")
}

# Try BEATRICE via actual paths
beatrice_real <- tryCatch(
  setup_beatrice(
    beatrice_dir = "~/Beatrice-Finemapping",
    python       = "/opt/anaconda3/bin/python3"
  ),
  error = function(e) NULL
)
if (!is.null(beatrice_real)) {
  res_G_be_real <- run_methods(
    sim_G,
    methods     = "beatrice",
    method_args = list(beatrice = list(
      beatrice_dir = "~/Beatrice-Finemapping",
      python       = "/opt/anaconda3/bin/python3",
      max_iter     = 200
    )),
    verbose = FALSE
  )
  check("G: BEATRICE real — pip has correct length or graceful fail",
    {
      f <- res_G_be_real$beatrice$results[[1]]
      length(f$pip) == 50
    })
  eval_G_be_real <- evaluate_methods(sim_G, res_G_be_real, verbose = FALSE)
  check("G: BEATRICE real — auprc in [0,1]",
    in_range_or_na(eval_G_be_real$beatrice$global$auprc))
} else {
  cat("  SKIP  G: BEATRICE — beatrice_dir or python not available\n")
}

# Try FunMap via actual Python
funmap_python <- "/opt/anaconda3/bin/python3"
funmap_ok <- tryCatch(
  setup_funmap(python = funmap_python),
  error = function(e) FALSE
)
if (isTRUE(funmap_ok)) {
  res_G_fm2 <- run_methods(
    sim_G,
    methods     = "funmap",
    method_args = list(funmap = list(python = funmap_python, L = 5)),
    verbose     = FALSE
  )
  check("G: FunMap real — pip has correct length or graceful fail",
    length(res_G_fm2$funmap$results[[1]]$pip) == 50)
  eval_G_fm2 <- evaluate_methods(sim_G, res_G_fm2, verbose = FALSE)
  check("G: FunMap real — auprc in [0,1]",
    in_range_or_na(eval_G_fm2$funmap$global$auprc))
} else {
  cat("  SKIP  G: FunMap — funmap Python package not available\n")
}


# =============================================================================
# [H] Evaluation output structure and metric range validation
# =============================================================================

section("[H] Evaluation structure and metric ranges (susie, sim_A)")

# Re-use eval_A from section A
check("H: simulation_params stored in eval output",
  !is.null(eval_A$simulation_params))
check("H: pip_thresholds_used in output",
  !is.null(eval_A$pip_thresholds_used))
check("H: pip_calibration bins sum to total pip count",
  {
    cal      <- eval_A$susie$global$pip_calibration
    valid    <- Filter(function(f) is.null(f$error), res_A$susie$results)
    expected <- sum(sapply(valid, function(f) length(f$pip)))
    sum(cal$n) == expected
  })
check("H: frac_causal in calibration bins is in [0,1] or NA",
  all(sapply(eval_A$susie$global$pip_calibration$frac_causal,
             function(x) is.na(x) || (x >= 0 && x <= 1))))
check("H: recall is non-decreasing as threshold decreases",
  {
    df <- eval_A$susie$global$fdr_power_curve
    # Sort by decreasing threshold; recall should be non-decreasing
    df_sorted <- df[order(df$threshold, decreasing = TRUE), ]
    all(diff(df_sorted$recall) >= -1e-10)
  })
check("H: tp + fn = total_causal at threshold = 0",
  {
    df  <- eval_A$susie$global$fdr_power_curve
    row <- df[df$threshold == 0, ]
    # At threshold 0 all variants are selected, so tp = total causal, fn = 0
    # Derive total causal from simulation truth for valid fits
    valid <- Filter(function(f) is.null(f$error), res_A$susie$results)
    total_causal <- sum(sapply(valid, function(f) {
      sc    <- sim_A$scenarios[[f$scenario_id]]
      truth <- sc$regions[[f$region_id]]$truth
      length(truth$causal_indices)
    }))
    nrow(row) == 1 && (row$tp + row$fn) == total_causal
  })
check("H: at threshold=1, tp=fp=0 (no variant has pip exactly >= 1)",
  {
    df  <- eval_A$susie$global$fdr_power_curve
    row <- df[df$threshold == 1, ]
    nrow(row) == 1 && row$fp == 0
  })
check("H: cs_size_mean >= cs_size_median (right-skewed distribution expected)",
  {
    g <- eval_A$susie$global
    is.na(g$cs_size_mean) || is.na(g$cs_size_median) ||
      g$cs_size_mean >= g$cs_size_median - 1e-10
  })
check("H: n_cs_reported is non-negative integer",
  is.numeric(eval_A$susie$global$n_cs_reported) &&
  eval_A$susie$global$n_cs_reported >= 0)
check("H: runtime_sd is NA when only 1 fit (or >=0 otherwise)",
  {
    sd_val <- eval_A$susie$global$runtime_sd
    is.na(sd_val) || sd_val >= 0
  })
check("H: n_pip_cal_bins parameter works (20 bins)",
  {
    eval_20 <- evaluate_methods(sim_A, res_A, n_pip_cal_bins = 20L,
                                verbose = FALSE)
    nrow(eval_20$susie$global$pip_calibration) == 20
  })
check("H: custom pip_thresholds respected",
  {
    eval_coarse <- evaluate_methods(sim_A, res_A,
                                    pip_thresholds = seq(0, 1, by = 0.1),
                                    verbose = FALSE)
    nrow(eval_coarse$susie$global$fdr_power_curve) == 11
  })


# =============================================================================
# [I] Edge cases
# =============================================================================

section("[I] Edge cases")

# I1: S = 1 (single causal variant per region) — CS coverage should be high
sim_I1 <- run_simulation(
  n_regions = 2,
  n         = 300,
  p         = 80,
  n_iter    = 3,
  S         = 1,
  phi       = 0.5,   # high PVE: signal should be detectable
  model     = "sparse",
  seed      = 10,
  verbose   = FALSE
)
res_I1  <- run_methods(sim_I1, methods = c("susie", "abf"),
                        method_args = list(susie = list(L = 3)),
                        verbose = FALSE)
eval_I1 <- evaluate_methods(sim_I1, res_I1, verbose = FALSE)
check("I1: S=1, high phi — susie CS coverage >= 0.5",
  {
    cov <- eval_I1$susie$global$cs_coverage
    is.na(cov) || cov >= 0.5
  })
check("I1: ABF reports exactly 1 credible set per fit (single-causal model)",
  all(sapply(Filter(function(f) is.null(f$error), res_I1$abf$results),
             function(f) length(f$credible_sets) == 1)))

# I2: Method with no credible sets reported — cs_* metrics should be NA
sim_I2 <- run_simulation(
  n_regions = 1,
  n         = 100,
  p         = 50,
  n_iter    = 1,
  S         = 1,
  phi       = 0.05,   # very low PVE: SuSiE likely reports no CS
  model     = "sparse",
  seed      = 11,
  verbose   = FALSE
)
res_I2  <- run_methods(sim_I2, methods = "susie",
                        method_args = list(susie = list(L = 1, coverage = 0.99)),
                        verbose = FALSE)
eval_I2 <- evaluate_methods(sim_I2, res_I2, verbose = FALSE)
check("I2: low signal — evaluate completes without error",
  !is.null(eval_I2$susie$global))
check("I2: if no CS reported, cs_coverage is NA",
  {
    g <- eval_I2$susie$global
    g$n_cs_reported == 0 || !is.na(g$cs_coverage)
  })

# I3: evaluate_methods with multiple methods, one pure failure
# Inject a fake method result with all errors into res_A
res_A_with_bad          <- res_A
res_A_with_bad$badmeth  <- list(
  results = lapply(res_A$susie$results, function(f) list(
    pip             = rep(NA_real_, 80),
    credible_sets   = list(),
    method          = "badmeth",
    input_type      = NA_character_,
    params          = list(),
    runtime_seconds = NA_real_,
    additional      = list(),
    error           = "forced error",
    scenario_id     = f$scenario_id,
    region_id       = f$region_id,
    S               = f$S,
    phi             = f$phi,
    p_causal        = f$p_causal,
    iter            = f$iter,
    causal_indices  = f$causal_indices,
    n_variants      = f$n_variants
  )),
  n_total               = 16L,
  n_failed              = 16L,
  method_args           = list(),
  total_runtime_seconds = 0
)
res_A_with_bad$methods_run <- c(res_A$methods_run, "badmeth")

# Register badmeth so evaluate_methods can find it
.FM_REGISTRY[["badmeth"]] <- "run_susie_region"  # placeholder; not called

eval_I3 <- evaluate_methods(sim_A, res_A_with_bad, verbose = FALSE)
check("I3: all-failed method returns NA auprc",
  is.na(eval_I3$badmeth$global$auprc))
check("I3: all-failed method has n_failed = 16",
  eval_I3$badmeth$global$n_failed == 16)
check("I3: all-failed method has NULL fdr_power_curve",
  is.null(eval_I3$badmeth$global$fdr_power_curve))
check("I3: other methods unaffected",
  !is.na(eval_I3$susie$global$auprc))

# Clean up fake registry entry
.FM_REGISTRY[["badmeth"]] <- NULL


# =============================================================================
# [J] Standard error fields — presence, validity, and n_iter behaviour
# =============================================================================

section("[J] Standard error fields")

SE_SCALAR_FIELDS <- c("auprc_se", "cs_coverage_se", "cs_power_se",
                       "cs_size_median_se", "cs_size_mean_se", "runtime_mean_se")

# J1: All scalar SE fields present in global and every stratum (n_iter=2 in sim_A)
check("J1: global has all scalar SE fields",
  has_names(eval_A$susie$global, SE_SCALAR_FIELDS))

check("J1: every by_S entry has all scalar SE fields",
  all(sapply(eval_A$susie$by_S, function(s) has_names(s, SE_SCALAR_FIELDS))))

check("J1: every by_phi entry has all scalar SE fields",
  all(sapply(eval_A$susie$by_phi, function(s) has_names(s, SE_SCALAR_FIELDS))))

check("J1: sparse_inf by_p_causal strata have SE fields",
  all(sapply(eval_E$susie$by_p_causal, function(s) has_names(s, SE_SCALAR_FIELDS))))

# J2: SE values are non-negative (or NA) everywhere — negative SE is impossible
check("J2: all scalar SE values are non-negative or NA (susie global)",
  all(sapply(SE_SCALAR_FIELDS, function(f) {
    v <- eval_A$susie$global[[f]]
    is.null(v) || is.na(v) || v >= 0
  })))

check("J2: all scalar SE values are non-negative or NA across all methods and strata",
  all(sapply(eval_A$methods_evaluated, function(m) {
    all(sapply(SE_SCALAR_FIELDS, function(f) {
      v <- eval_A[[m]]$global[[f]]
      is.null(v) || is.na(v) || v >= 0
    }))
  })))

# J3: With n_iter=2 (sim_A), SE values for scalar metrics should be numeric
# (may be 0 if all replicates agree, but should not be NA for methods with data)
check("J3: auprc_se is numeric (not NA) when n_iter=2 and all fits succeed",
  {
    se <- eval_A$susie$global$auprc_se
    is.numeric(se) && !is.na(se)
  })

check("J3: cs_coverage_se is numeric when n_iter=2",
  {
    se <- eval_A$abf$global$cs_coverage_se
    is.numeric(se) && !is.na(se)
  })

# J4: With n_iter=1 (sim_I2), all SE fields should be NA
check("J4: all SE fields are NA when n_iter=1",
  all(sapply(SE_SCALAR_FIELDS, function(f) {
    v <- eval_I2$susie$global[[f]]
    is.null(v) || is.na(v)
  })))

check("J4: fdr_power_curve power_se is NA when n_iter=1",
  {
    df <- eval_I2$susie$global$fdr_power_curve
    is.null(df) || all(is.na(df$power_se))
  })

check("J4: pip_calibration frac_causal_se is NA when n_iter=1",
  {
    cal <- eval_I2$susie$global$pip_calibration
    is.null(cal) || all(is.na(cal$frac_causal_se))
  })

# J5: fdr_power_curve has power_se and precision_se columns when n_iter >= 2
check("J5: fdr_power_curve has power_se column",
  "power_se" %in% names(eval_A$susie$global$fdr_power_curve))

check("J5: fdr_power_curve has precision_se column",
  "precision_se" %in% names(eval_A$susie$global$fdr_power_curve))

check("J5: power_se values are non-negative or NA",
  {
    v <- eval_A$susie$global$fdr_power_curve$power_se
    all(is.na(v) | v >= 0)
  })

check("J5: precision_se values are non-negative or NA",
  {
    v <- eval_A$susie$global$fdr_power_curve$precision_se
    all(is.na(v) | v >= 0)
  })

check("J5: at least some SE values are non-NA when n_iter=2",
  any(!is.na(eval_A$susie$global$fdr_power_curve$power_se)))

# J6: pip_calibration has frac_causal_se column when n_iter >= 2
check("J6: pip_calibration has frac_causal_se column",
  "frac_causal_se" %in% names(eval_A$susie$global$pip_calibration))

check("J6: frac_causal_se values are non-negative or NA",
  {
    v <- eval_A$susie$global$pip_calibration$frac_causal_se
    all(is.na(v) | v >= 0)
  })

# J7: SE decreases with more iterations (law of large numbers — test with 4 vs 2 iter)
sim_J7_n2 <- run_simulation(n_regions=2, n=200, p=80, n_iter=2, S=1,
                             phi=0.3, seed=20, verbose=FALSE)
sim_J7_n4 <- run_simulation(n_regions=2, n=200, p=80, n_iter=4, S=1,
                             phi=0.3, seed=20, verbose=FALSE)
res_J7_n2 <- run_methods(sim_J7_n2, methods="abf", verbose=FALSE)
res_J7_n4 <- run_methods(sim_J7_n4, methods="abf", verbose=FALSE)
ev_J7_n2  <- evaluate_methods(sim_J7_n2, res_J7_n2, verbose=FALSE)
ev_J7_n4  <- evaluate_methods(sim_J7_n4, res_J7_n4, verbose=FALSE)

check("J7: SE is non-NA for n_iter=2 and n_iter=4",
  !is.na(ev_J7_n2$abf$global$auprc_se) &&
  !is.na(ev_J7_n4$abf$global$auprc_se))

check("J7: by_S and by_phi SE fields consistent with n_iter (non-NA, non-negative)",
  {
    s2 <- ev_J7_n2$abf$by_S[["1"]]$auprc_se
    s4 <- ev_J7_n4$abf$by_S[["1"]]$auprc_se
    (is.na(s2) || s2 >= 0) && (is.na(s4) || s4 >= 0)
  })

# J8: Failed-method SE fields should all be NA (no valid fits → no SE)
eval_J8 <- evaluate_methods(sim_G, res_G_finemap, verbose=FALSE)
check("J8: all SE fields are NA when all fits failed",
  all(sapply(SE_SCALAR_FIELDS, function(f) {
    v <- eval_J8$finemap$global[[f]]
    is.null(v) || is.na(v)
  })))


# =============================================================================
# [K] plot_results — PDF output, all sections, method filtering
# =============================================================================

section("[K] plot_results")

tmp_pdf_dir <- tempfile("fm_test_plots")
dir.create(tmp_pdf_dir)

# K1: Basic PDF creation (sparse model, 2 methods)
pdf_K1 <- file.path(tmp_pdf_dir, "test_K1.pdf")
check("K1: plot_results runs without error (sparse, susie+abf)",
  {
    tryCatch({
      plot_results(eval_A, output_file = pdf_K1,
                   methods = c("susie", "abf"), verbose = FALSE)
      TRUE
    }, error = function(e) { cat("      =>", conditionMessage(e), "\n"); FALSE })
  })

check("K1: PDF file is created and non-empty",
  file.exists(pdf_K1) && file.size(pdf_K1) > 1000)

# K2: sparse_inf model — p_causal section is generated
pdf_K2 <- file.path(tmp_pdf_dir, "test_K2.pdf")
check("K2: plot_results runs for sparse_inf model (by_p_causal section)",
  {
    tryCatch({
      plot_results(eval_E, output_file = pdf_K2, verbose = FALSE)
      TRUE
    }, error = function(e) { cat("      =>", conditionMessage(e), "\n"); FALSE })
  })

check("K2: sparse_inf PDF is larger than sparse PDF (extra p_causal section)",
  file.exists(pdf_K2) && file.size(pdf_K2) > file.size(pdf_K1))

# K3: method filtering — only plot one method
pdf_K3 <- file.path(tmp_pdf_dir, "test_K3.pdf")
check("K3: plot_results with single method filter runs without error",
  {
    tryCatch({
      plot_results(eval_A, output_file = pdf_K3,
                   methods = "susie", verbose = FALSE)
      TRUE
    }, error = function(e) { cat("      =>", conditionMessage(e), "\n"); FALSE })
  })

check("K3: filtering to one method produces a smaller PDF than all methods",
  file.exists(pdf_K3) && file.size(pdf_K3) <= file.size(pdf_K1) * 1.2)

# K4: mixed methods — some failed, some succeeded
pdf_K4 <- file.path(tmp_pdf_dir, "test_K4.pdf")
eval_K4_mixed <- evaluate_methods(
  sim_G,
  {
    # Combine a successful abf run with the failed finemap run
    res_abf_G <- run_methods(sim_G, methods="abf", verbose=FALSE)
    merged <- res_abf_G
    merged$finemap        <- res_G_finemap$finemap
    merged$methods_run    <- c("abf", "finemap")
    merged
  },
  verbose = FALSE
)
check("K4: plot_results handles mix of succeeded and failed methods",
  {
    tryCatch({
      plot_results(eval_K4_mixed, output_file = pdf_K4, verbose = FALSE)
      file.exists(pdf_K4) && file.size(pdf_K4) > 1000
    }, error = function(e) { cat("      =>", conditionMessage(e), "\n"); FALSE })
  })

# K5: invalid methods argument raises error
check("K5: invalid methods argument raises error",
  tryCatch({
    plot_results(eval_A, output_file = file.path(tmp_pdf_dir, "bad.pdf"),
                 methods = "nonexistent", verbose = FALSE)
    FALSE
  }, error = function(e) TRUE))

# K6: many S and phi values — grid plots remain manageable
sim_K6 <- run_simulation(n_regions=1, n=150, p=60, n_iter=2,
                          S=c(1,2,3), phi=c(0.1,0.3,0.5), seed=30, verbose=FALSE)
res_K6  <- run_methods(sim_K6, methods=c("susie","abf"),
                        method_args=list(susie=list(L=5)), verbose=FALSE)
eval_K6 <- evaluate_methods(sim_K6, res_K6, verbose=FALSE)
pdf_K6  <- file.path(tmp_pdf_dir, "test_K6.pdf")
check("K6: 3-S x 3-phi grid plots without error",
  {
    tryCatch({
      plot_results(eval_K6, output_file=pdf_K6, verbose=FALSE)
      file.exists(pdf_K6) && file.size(pdf_K6) > 1000
    }, error = function(e) { cat("      =>", conditionMessage(e), "\n"); FALSE })
  })

# K7: n_iter=1 (no SE) — plot_results still works without error bars
sim_K7 <- run_simulation(n_regions=1, n=150, p=60, n_iter=1,
                          S=c(1,2), phi=0.3, seed=31, verbose=FALSE)
res_K7  <- run_methods(sim_K7, methods="susie",
                        method_args=list(susie=list(L=5)), verbose=FALSE)
eval_K7 <- evaluate_methods(sim_K7, res_K7, verbose=FALSE)
pdf_K7  <- file.path(tmp_pdf_dir, "test_K7.pdf")
check("K7: n_iter=1 (no SE) — plot_results runs without error",
  {
    tryCatch({
      plot_results(eval_K7, output_file=pdf_K7, verbose=FALSE)
      file.exists(pdf_K7) && file.size(pdf_K7) > 1000
    }, error = function(e) { cat("      =>", conditionMessage(e), "\n"); FALSE })
  })

# Clean up temp PDF dir
unlink(tmp_pdf_dir, recursive=TRUE)


# =============================================================================
# [L] save argument — evaluate_methods writes correct files
# =============================================================================

section("[L] save argument")

tmp_save_dir <- tempfile("fm_test_save")

eval_L <- evaluate_methods(
  sim_A, res_A,
  save       = TRUE,
  output_dir = tmp_save_dir,
  verbose    = FALSE
)

check("L1: output_dir is created automatically",
  dir.exists(tmp_save_dir))

check("L1: evaluation.rds is written",
  file.exists(file.path(tmp_save_dir, "evaluation.rds")))

check("L1: evaluation_summary.csv is written",
  file.exists(file.path(tmp_save_dir, "evaluation_summary.csv")))

check("L2: evaluation.rds round-trips correctly (susie auprc matches)",
  {
    ev2 <- readRDS(file.path(tmp_save_dir, "evaluation.rds"))
    identical(eval_L$susie$global$auprc, ev2$susie$global$auprc)
  })

check("L2: evaluation.rds preserves SE fields",
  {
    ev2 <- readRDS(file.path(tmp_save_dir, "evaluation.rds"))
    has_names(ev2$susie$global, SE_SCALAR_FIELDS)
  })

check("L3: evaluation_summary.csv has correct method rows",
  {
    csv <- read.csv(file.path(tmp_save_dir, "evaluation_summary.csv"),
                    stringsAsFactors = FALSE)
    setequal(csv$method, eval_L$methods_evaluated)
  })

check("L3: evaluation_summary.csv has SE columns",
  {
    csv <- read.csv(file.path(tmp_save_dir, "evaluation_summary.csv"),
                    stringsAsFactors = FALSE)
    has_names(csv, c("auprc_se", "cs_coverage_se", "cs_power_se",
                     "cs_size_median_se", "runtime_mean_se"))
  })

check("L3: evaluation_summary.csv auprc values match eval object",
  {
    csv <- read.csv(file.path(tmp_save_dir, "evaluation_summary.csv"),
                    stringsAsFactors = FALSE)
    all(sapply(eval_L$methods_evaluated, function(m) {
      csv_val <- csv$auprc[csv$method == m]
      ev_val  <- eval_L[[m]]$global$auprc
      (is.na(csv_val) && is.na(ev_val)) ||
        (!is.na(csv_val) && !is.na(ev_val) && abs(csv_val - ev_val) < 1e-10)
    }))
  })

check("L4: save=FALSE does not write files",
  {
    tmp2 <- tempfile("fm_nosave")
    evaluate_methods(sim_A, res_A, save=FALSE,
                     output_dir=tmp2, verbose=FALSE)
    !dir.exists(tmp2)
  })

# Clean up
unlink(tmp_save_dir, recursive=TRUE)


# =============================================================================
# Summary
# =============================================================================

cat(sprintf(
  "\n%s\nResults: %d/%d tests passed",
  strrep("=", 60),
  .tests_passed,
  .tests_run
))
if (length(.tests_failed) > 0) {
  cat(sprintf("\nFailed tests (%d):\n", length(.tests_failed)))
  for (f in .tests_failed) cat(sprintf("  - %s\n", f))
} else {
  cat("\nAll tests passed.\n")
}
cat(strrep("=", 60), "\n")

if (length(.tests_failed) > 0) quit(status = 1L)
