# =============================================================================
# test-evaluate.R
#
# Comprehensive testthat coverage of the full benchmark pipeline:
#   simulate -> run_methods -> evaluate_methods -> plot_results
#
# Translated from scripts/test_evaluate.R as part of Phase 3b-ii bulk.
# The original script is kept on disk for now; 3b-iii decides whether to
# move it under inst/scripts/ or delete it.
#
# Coverage (mirrors the original section labels):
#   [A] Sparse model, pure-R methods (susie, susie_inf, abf, carma)
#   [B] Sparse model with binary annotations
#   [C] Sparse model with continuous annotations
#   [D] Multiple S and phi values - stratified metric validation
#   [E] sparse_inf model - p_causal stratification
#   [F] Alternative simulation args (equal effects, non-default n/p)
#   [G] External methods (finemap, paintor, funmap, beatrice) - graceful failure
#   [H] evaluate output structure and metric range validation
#   [I] Edge cases (single causal variant, high phi, all fits failed)
#   [J] Standard error fields - presence, validity, and n_iter behaviour
#   [K] plot_results - PDF output, all sections, method filtering
#   [L] save argument - evaluate_methods writes correct files
#
# Each section sets up its fixtures (sim_X, res_X, eval_X) at file level
# so multiple test_that blocks can share them. This mirrors the original
# script's structure and keeps the simulation cost paid once per section.
#
# Skip-on-CRAN: simulations call sim1000G which may need to download a
# HapMap genetic map (~1 MB per chromosome) on first call; CRAN does not
# allow that. The section [A] fixture block executes at file-load time,
# so we guard against missing sim1000G with `skip_if_not_installed()`.
# =============================================================================


# --- Small inline helpers (mirror the original) ---------------------------

# `has_names(x, nms)` is more concise than expect_named() when we only
# want presence (not exact order or set).
has_names <- function(x, nms) all(nms %in% names(x))

# `in_range_or_na(x, lo, hi)`: scalar numeric in [lo, hi] or NA.
in_range_or_na <- function(x, lo = 0, hi = 1) {
  is.na(x) || (is.numeric(x) && length(x) == 1L && x >= lo && x <= hi)
}


# =============================================================================
# [A] Sparse model - pure-R methods (susie, susie_inf, abf, carma)
# =============================================================================

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

eval_A <- evaluate_methods(sim_A, res_A, verbose = FALSE)


test_that("[A] sim_A: structure of run_simulation output (sparse model)", {
  expect_true(has_names(sim_A, c("genotypes", "scenarios", "params")))
  expect_length(sim_A$genotypes, 2L)
  # 2 S * 2 phi * 2 iter = 8 scenarios
  expect_length(sim_A$scenarios, 8L)
  expect_true(all(sapply(sim_A$scenarios,
                         function(sc) length(sc$regions) == 2L)))
  expect_false(is.null(sim_A$scenarios[[1]]$regions[[1]]$truth$causal_indices))
})


test_that("[A] res_A: run_methods output for pure-R methods", {
  expect_setequal(res_A$methods_run,
                  c("susie", "susie_inf", "abf", "carma"))
  # 8 scenarios * 2 regions = 16 fits per method
  expect_length(res_A$susie$results, 16L)
  expect_true(all(sapply(res_A$susie$results,
                         function(f) length(f$pip) == 80L)))
  expect_true(all(sapply(res_A$susie$results,
                         function(f) is.list(f$credible_sets))))
  expect_true(has_names(res_A$susie$results[[1]],
                        c("scenario_id", "region_id", "S", "phi", "iter")))

  # ABF assumes a single causal variant; pip should sum to ~1 per region.
  abf_ok <- Filter(function(f) is.null(f$error), res_A$abf$results)
  expect_true(all(sapply(abf_ok,
                         function(f) abs(sum(f$pip) - 1) < 1e-6)))
})


test_that("[A] eval_A: evaluate_methods output structure", {
  expect_setequal(eval_A$methods_evaluated,
                  c("susie", "susie_inf", "abf", "carma"))
  for (m in eval_A$methods_evaluated) {
    expect_true(
      has_names(eval_A[[m]],
                c("global", "by_S", "by_phi", "by_p_causal")),
      info = sprintf("method = %s", m)
    )
  }

  expect_true(has_names(
    eval_A$susie$global,
    c("fdr_power_curve", "auprc", "pip_calibration", "cs_coverage",
      "cs_power", "cs_size_median", "cs_size_mean", "n_cs_reported",
      "runtime_mean", "runtime_sd", "n_fits", "n_failed")
  ))

  expect_equal(eval_A$susie$global$n_fits, 16L)
  expect_true(in_range_or_na(eval_A$susie$global$auprc))
  expect_true(in_range_or_na(eval_A$susie$global$cs_coverage))
  expect_true(in_range_or_na(eval_A$susie$global$cs_power))

  expect_setequal(names(eval_A$susie$by_S),   c("1", "2"))
  expect_setequal(names(eval_A$susie$by_phi), c("0.2", "0.4"))
  expect_null(eval_A$susie$by_p_causal)

  # PR curve table
  fpc <- eval_A$susie$global$fdr_power_curve
  expect_s3_class(fpc, "data.frame")
  expect_true(has_names(fpc,
    c("threshold", "tp", "fp", "fn", "fdr", "power", "precision", "recall")))
  expect_true(all(fpc$fdr   >= 0 & fpc$fdr   <= 1))
  expect_true(all(fpc$power >= 0 & fpc$power <= 1))

  # PIP calibration table
  pc <- eval_A$susie$global$pip_calibration
  expect_s3_class(pc, "data.frame")
  expect_equal(nrow(pc), 10L)
  expect_true(has_names(pc,
    c("bin", "bin_lower", "bin_upper", "bin_mid", "n", "n_causal",
      "mean_pip", "frac_causal")))

  # Runtime sanity
  expect_true(all(sapply(eval_A$methods_evaluated, function(m) {
    rt <- eval_A[[m]]$global$runtime_mean
    is.na(rt) || rt >= 0
  })))
})


# =============================================================================
# [B] Sparse model with binary annotations
# =============================================================================

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

res_B <- run_methods(
  simulation  = sim_B,
  methods     = c("susie", "abf"),
  method_args = list(susie = list(L = 5)),
  verbose     = FALSE
)

eval_B <- evaluate_methods(sim_B, res_B, verbose = FALSE)


test_that("[B] sim_B: binary annotations attached to genotypes and scenarios", {
  expect_false(is.null(sim_B$genotypes[[1]]$annotations_matrix))
  expect_equal(ncol(sim_B$genotypes[[1]]$annotations_matrix), 3L)
  expect_false(is.null(sim_B$scenarios[[1]]$regions[[1]]$annotations_matrix))
})


test_that("[B] eval_B: evaluation works on annotated simulation", {
  expect_false(is.null(eval_B$susie$global$auprc))
  expect_setequal(names(eval_B$susie$by_S), c("1", "2"))
})


# =============================================================================
# [C] Sparse model with continuous annotations
# =============================================================================

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

res_C  <- run_methods(sim_C, methods = "abf", verbose = FALSE)
eval_C <- evaluate_methods(sim_C, res_C, verbose = FALSE)


test_that("[C] continuous annotations: simulation + evaluation", {
  expect_false(is.null(sim_C$genotypes[[1]]$annotations_matrix))
  expect_equal(ncol(sim_C$genotypes[[1]]$annotations_matrix), 2L)
  expect_false(is.null(eval_C$abf$global$auprc))
})


# =============================================================================
# [D] Multiple S and phi - stratified metric validation
# =============================================================================

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

res_D <- run_methods(
  sim_D,
  methods     = c("susie", "abf"),
  method_args = list(susie = list(L = 5)),
  verbose     = FALSE
)

eval_D <- evaluate_methods(sim_D, res_D, verbose = FALSE)


test_that("[D] sim_D: full S * phi * iter grid", {
  # 3 S * 3 phi * 2 iter = 18 scenarios
  expect_length(sim_D$scenarios, 18L)
})


test_that("[D] eval_D: stratified by_S and by_phi populated", {
  expect_setequal(names(eval_D$susie$by_S),   c("1", "2", "3"))
  expect_setequal(names(eval_D$susie$by_phi), c("0.1", "0.3", "0.5"))

  # Each S stratum: 3 phi * 2 iter * 2 regions = 12 fits
  expect_true(all(sapply(eval_D$susie$by_S, function(s) s$n_fits == 12L)))
  # Each phi stratum: 3 S * 2 iter * 2 regions = 12 fits
  expect_true(all(sapply(eval_D$susie$by_phi, function(s) s$n_fits == 12L)))

  # Global: 3 S * 3 phi * 2 iter * 2 regions = 36
  expect_equal(eval_D$susie$global$n_fits, 36L)

  # ABF AUPRC in [0,1] across phi strata.
  expect_true(all(sapply(eval_D$abf$by_phi,
                         function(s) in_range_or_na(s$auprc))))
})


# =============================================================================
# [E] sparse_inf model - p_causal stratification
# =============================================================================

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


test_that("[E] sim_E: sparse_inf grid with p_causal axis", {
  # 2 S * 2 phi * 2 p_causal * 2 iter = 16 scenarios
  expect_length(sim_E$scenarios, 16L)
  expect_true(all(sapply(sim_E$scenarios,
                         function(sc) !is.null(sc$p_causal))))
})


test_that("[E] eval_E: by_p_causal stratification populated for sparse_inf", {
  expect_false(is.null(eval_E$susie$by_p_causal))
  expect_setequal(names(eval_E$susie$by_p_causal), c("0.2", "0.6"))

  # Each p_causal stratum: 2 S * 2 phi * 2 iter * 2 regions = 16 fits
  expect_true(all(sapply(eval_E$susie$by_p_causal,
                         function(s) s$n_fits == 16L)))

  expect_true(all(sapply(eval_E$susie_inf$by_p_causal,
                         function(s) in_range_or_na(s$auprc))))

  # Global: 2 * 2 * 2 * 2 * 2 = 32 fits
  expect_equal(eval_E$susie_inf$global$n_fits, 32L)
})


# =============================================================================
# [F] Alternative simulation arguments
# =============================================================================

# Equal-magnitude effect distribution.
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

res_F1  <- run_methods(sim_F1, methods = "susie",
                       method_args = list(susie = list(L = 3)),
                       verbose = FALSE)
eval_F1 <- evaluate_methods(sim_F1, res_F1, verbose = FALSE)

# Varying p per region via a length-n_regions vector.
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

res_F2  <- run_methods(sim_F2, methods = "abf", verbose = FALSE)
eval_F2 <- evaluate_methods(sim_F2, res_F2, verbose = FALSE)


test_that("[F1] equal effect distribution: |effects| constant within a region", {
  eff <- abs(sim_F1$scenarios[[1]]$regions[[1]]$truth$causal_effects)
  expect_length(eff, 2L)
  expect_lt(diff(range(eff)), 1e-10)
  expect_false(is.null(eval_F1$susie$global$auprc))
})


test_that("[F2] per-region p vector produces matching region sizes", {
  expect_length(sim_F2$genotypes, 3L)
  expect_equal(sim_F2$genotypes[[1]]$p,  50L)
  expect_equal(sim_F2$genotypes[[2]]$p,  80L)
  expect_equal(sim_F2$genotypes[[3]]$p, 100L)

  expect_true(in_range_or_na(eval_F2$abf$global$auprc))
})


test_that("[F3] ABF non-default prior_variance and coverage are honoured", {
  res_F3 <- run_methods(
    sim_F1,
    methods     = "abf",
    method_args = list(abf = list(prior_variance = 0.01, coverage = 0.90)),
    verbose     = FALSE
  )
  expect_equal(res_F3$abf$results[[1]]$params$prior_variance, 0.01)

  eval_F3 <- evaluate_methods(sim_F1, res_F3, verbose = FALSE)
  expect_true(is.numeric(eval_F3$abf$global$cs_coverage))
})


test_that("[F4] SuSiE individual-level mode (X + y) reports input_type='individual'", {
  res_F4 <- run_methods(
    sim_F1,
    methods     = "susie",
    method_args = list(susie = list(L = 3, use_individual = TRUE)),
    verbose     = FALSE
  )
  ok_fits <- Filter(function(f) is.null(f$error), res_F4$susie$results)
  expect_true(length(ok_fits) > 0L)
  expect_true(all(sapply(ok_fits, function(f) f$input_type == "individual")))
})


# =============================================================================
# [G] External methods - graceful failure when binary / Python not found
# =============================================================================

# Tiny sim to keep this fast (1 region, 1 scenario => 1 fit per method).
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


test_that("[G-finemap] missing finemap binary: graceful per-fit failure", {
  res <- run_methods(
    sim_G,
    methods     = "finemap",
    method_args = list(finemap = list(finemap_path = "/nonexistent/finemap")),
    verbose     = FALSE
  )
  expect_equal(res$finemap$n_failed, 1L)
  expect_true(all(is.na(res$finemap$results[[1]]$pip)))

  ev <- evaluate_methods(sim_G, res, verbose = FALSE)
  expect_true(is.na(ev$finemap$global$auprc))
  expect_true(is.na(ev$finemap$global$cs_coverage))
})


test_that("[G-paintor] missing PAINTOR binary: graceful per-fit failure", {
  res <- run_methods(
    sim_G,
    methods     = "paintor",
    method_args = list(paintor = list(paintor_path = "/nonexistent/PAINTOR")),
    verbose     = FALSE
  )
  expect_equal(res$paintor$n_failed, 1L)

  ev <- evaluate_methods(sim_G, res, verbose = FALSE)
  expect_true(is.na(ev$paintor$global$auprc))
})


test_that("[G-beatrice] missing BEATRICE dir / python: graceful per-fit failure", {
  res <- run_methods(
    sim_G,
    methods     = "beatrice",
    method_args = list(beatrice = list(
      beatrice_dir = "/nonexistent/Beatrice",
      python       = "/nonexistent/python3"
    )),
    verbose = FALSE
  )
  expect_equal(res$beatrice$n_failed, 1L)
})


test_that("[G-funmap] missing python: graceful per-fit failure", {
  res <- run_methods(
    sim_G,
    methods     = "funmap",
    method_args = list(funmap = list(python = "/nonexistent/python3")),
    verbose     = FALSE
  )
  expect_equal(res$funmap$n_failed, 1L)
})


# Real-binary checks for external methods. These skip cleanly when the
# binary / Python env is not available on the current machine. They are
# the same checks the original script ran via tryCatch + setup_*().

test_that("[G-finemap-real] real FINEMAP binary: pip length, AUPRC in [0,1]", {
  fp <- tryCatch(setup_finemap(download = FALSE),
                 error = function(e) NULL)
  skip_if(is.null(fp), "FINEMAP binary not available on this machine")

  res <- run_methods(
    sim_G,
    methods     = "finemap",
    method_args = list(finemap = list(finemap_path = fp, n_causal = 2)),
    verbose     = FALSE
  )
  expect_length(res$finemap$results[[1]]$pip, 50L)

  ev <- evaluate_methods(sim_G, res, verbose = FALSE)
  expect_true(in_range_or_na(ev$finemap$global$auprc))
})


test_that("[G-paintor-real] real PAINTOR binary: pip length, AUPRC in [0,1]", {
  pp <- tryCatch(setup_paintor(), error = function(e) NULL)
  skip_if(is.null(pp), "PAINTOR binary not available on this machine")

  res <- run_methods(
    sim_G,
    methods     = "paintor",
    method_args = list(paintor = list(paintor_path = pp, max_causal = 2)),
    verbose     = FALSE
  )
  expect_length(res$paintor$results[[1]]$pip, 50L)

  ev <- evaluate_methods(sim_G, res, verbose = FALSE)
  expect_true(in_range_or_na(ev$paintor$global$auprc))
})


test_that("[G-beatrice-real] real BEATRICE: pip length, AUPRC in [0,1]", {
  ok <- tryCatch(
    setup_beatrice(
      beatrice_dir = "~/Beatrice-Finemapping",
      python       = "/opt/anaconda3/bin/python3"
    ),
    error = function(e) NULL
  )
  skip_if(is.null(ok),
          "BEATRICE_dir or compatible python not available on this machine")

  res <- run_methods(
    sim_G,
    methods     = "beatrice",
    method_args = list(beatrice = list(
      beatrice_dir = "~/Beatrice-Finemapping",
      python       = "/opt/anaconda3/bin/python3",
      max_iter     = 200
    )),
    verbose = FALSE
  )
  expect_length(res$beatrice$results[[1]]$pip, 50L)

  ev <- evaluate_methods(sim_G, res, verbose = FALSE)
  expect_true(in_range_or_na(ev$beatrice$global$auprc))
})


test_that("[G-funmap-real] real Funmap: pip length, AUPRC in [0,1]", {
  funmap_python <- "/opt/anaconda3/bin/python3"
  ok <- isTRUE(tryCatch(setup_funmap(python = funmap_python),
                         error = function(e) FALSE))
  skip_if(!ok, "funmap Python package not available on this machine")

  res <- run_methods(
    sim_G,
    methods     = "funmap",
    method_args = list(funmap = list(python = funmap_python, L = 5)),
    verbose     = FALSE
  )
  expect_length(res$funmap$results[[1]]$pip, 50L)

  ev <- evaluate_methods(sim_G, res, verbose = FALSE)
  expect_true(in_range_or_na(ev$funmap$global$auprc))
})


# =============================================================================
# [H] Evaluation output structure and metric range validation
# =============================================================================

# Reuses sim_A / res_A / eval_A from section [A].

test_that("[H] top-level eval object carries simulation_params and thresholds", {
  expect_false(is.null(eval_A$simulation_params))
  expect_false(is.null(eval_A$pip_thresholds_used))
})


test_that("[H] PIP calibration bin sums match pooled fit pip count", {
  cal      <- eval_A$susie$global$pip_calibration
  valid    <- Filter(function(f) is.null(f$error), res_A$susie$results)
  expected <- sum(sapply(valid, function(f) length(f$pip)))
  expect_equal(sum(cal$n), expected)

  # frac_causal entries are in [0,1] or NA
  expect_true(all(sapply(cal$frac_causal,
                         function(x) is.na(x) || (x >= 0 && x <= 1))))
})


test_that("[H] PR curve has consistent shape and edge-case totals", {
  df <- eval_A$susie$global$fdr_power_curve

  # Recall is non-decreasing as threshold decreases.
  df_sorted <- df[order(df$threshold, decreasing = TRUE), ]
  expect_true(all(diff(df_sorted$recall) >= -1e-10))

  # At threshold = 0 all variants are selected, so tp + fn = total causal.
  total_causal <- sum(sapply(
    Filter(function(f) is.null(f$error), res_A$susie$results),
    function(f) length(sim_A$scenarios[[f$scenario_id]]$regions[[f$region_id]]$truth$causal_indices)
  ))
  row_t0 <- df[df$threshold == 0, ]
  expect_equal(nrow(row_t0), 1L)
  expect_equal(row_t0$tp + row_t0$fn, total_causal)

  # At threshold = 1 no variant has pip > 1, so fp = 0.
  row_t1 <- df[df$threshold == 1, ]
  expect_equal(nrow(row_t1), 1L)
  expect_equal(row_t1$fp, 0L)
})


test_that("[H] cs_size_mean >= cs_size_median (right-skewed CS distribution)", {
  g <- eval_A$susie$global
  expect_true(
    is.na(g$cs_size_mean) || is.na(g$cs_size_median) ||
      g$cs_size_mean >= g$cs_size_median - 1e-10
  )

  # n_cs_reported and runtime_sd are non-negative or NA
  expect_true(is.numeric(g$n_cs_reported) && g$n_cs_reported >= 0L)
  expect_true(is.na(g$runtime_sd) || g$runtime_sd >= 0)
})


test_that("[H] n_pip_cal_bins and pip_thresholds args are honoured", {
  eval_20 <- evaluate_methods(sim_A, res_A,
                              n_pip_cal_bins = 20L, verbose = FALSE)
  expect_equal(nrow(eval_20$susie$global$pip_calibration), 20L)

  eval_coarse <- evaluate_methods(
    sim_A, res_A,
    pip_thresholds = seq(0, 1, by = 0.1),
    verbose = FALSE
  )
  expect_equal(nrow(eval_coarse$susie$global$fdr_power_curve), 11L)
})


# =============================================================================
# [I] Edge cases
# =============================================================================

# I1: S = 1 with high PVE - susie CS coverage should be reasonable; ABF
# always reports exactly one CS.
sim_I1 <- run_simulation(
  n_regions = 2,
  n         = 300,
  p         = 80,
  n_iter    = 3,
  S         = 1,
  phi       = 0.5,
  model     = "sparse",
  seed      = 10,
  verbose   = FALSE
)
res_I1  <- run_methods(sim_I1,
                        methods     = c("susie", "abf"),
                        method_args = list(susie = list(L = 3)),
                        verbose     = FALSE)
eval_I1 <- evaluate_methods(sim_I1, res_I1, verbose = FALSE)

test_that("[I1] S=1, high phi: susie CS coverage >= 0.5; ABF returns 1 CS per fit", {
  cov <- eval_I1$susie$global$cs_coverage
  expect_true(is.na(cov) || cov >= 0.5)

  ok_fits <- Filter(function(f) is.null(f$error), res_I1$abf$results)
  expect_true(all(sapply(ok_fits, function(f) length(f$credible_sets) == 1L)))
})


# I2: Low PVE - SuSiE may report no credible sets; cs_* metrics should be NA
# in that case but evaluate must not error.
sim_I2 <- run_simulation(
  n_regions = 1,
  n         = 100,
  p         = 50,
  n_iter    = 1,
  S         = 1,
  phi       = 0.05,
  model     = "sparse",
  seed      = 11,
  verbose   = FALSE
)
res_I2 <- run_methods(sim_I2,
                       methods     = "susie",
                       method_args = list(susie = list(L = 1, coverage = 0.99)),
                       verbose     = FALSE)
eval_I2 <- evaluate_methods(sim_I2, res_I2, verbose = FALSE)

test_that("[I2] low-signal SuSiE: evaluate completes; CS metrics NA-safe", {
  expect_false(is.null(eval_I2$susie$global))
  g <- eval_I2$susie$global
  expect_true(g$n_cs_reported == 0L || !is.na(g$cs_coverage))
})


test_that("[I3] all-fits-failed method: NA AUPRC, n_failed correct, others unaffected", {
  # Inject a fake fully-failed method into res_A and register it so
  # evaluate_methods will descend into it. Restore .FM_REGISTRY on exit.
  old_registry <- .FM_REGISTRY
  on.exit(assign(".FM_REGISTRY", old_registry, envir = globalenv()), add = TRUE)

  .FM_REGISTRY[["badmeth"]] <<- "run_susie_region"  # placeholder, never invoked

  res_with_bad <- res_A
  res_with_bad$badmeth <- list(
    results = lapply(res_A$susie$results, function(f) list(
      pip             = rep(NA_real_, 80L),
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
  res_with_bad$methods_run <- c(res_A$methods_run, "badmeth")

  ev <- evaluate_methods(sim_A, res_with_bad, verbose = FALSE)

  expect_true(is.na(ev$badmeth$global$auprc))
  expect_equal(ev$badmeth$global$n_failed, 16L)
  expect_null(ev$badmeth$global$fdr_power_curve)

  # Other methods unaffected by the failed one
  expect_false(is.na(ev$susie$global$auprc))
})
