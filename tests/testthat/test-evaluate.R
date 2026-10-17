# =============================================================================
# test-evaluate.R
#
# Comprehensive testthat coverage of the full benchmark pipeline:
#   simulate -> run_methods -> evaluate_methods -> plot_results
#
# Translated from the former scripts/test_evaluate.R, which was removed
# once this testthat version superseded it.
#
# Coverage (mirrors the original section labels):
#   [A] Sparse model, pure-R methods (susie, susie_inf, abf, carma)
#   [B] Sparse model with binary annotations
#   [C] Sparse model with continuous annotations
#   [D] Multiple S and phi values
#   [E] sparse_inf model
#   [F] Alternative simulation args (equal effects, non-default n/p)
#   [G] External methods (finemap, paintor, funmap, beatrice) - graceful failure
#   [H] evaluate output structure and metric range validation
#   [I] Edge cases (single causal variant, high phi, all fits failed)
#   [L] save argument - evaluate_methods writes correct files
#   [M] Per-region metrics: nothing pooled across regions
#   [N] FDR and calibration are region-weighted, not call-weighted
#   [O] The unit is one region of one iteration
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
  expect_false(is.null(eval_C$abf$global$ap))
})


# =============================================================================
# [D] Multiple S and phi
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


# =============================================================================
# [E] sparse_inf model
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
  expect_false(is.null(eval_F1$susie$global$ap))
})


test_that("[F2] per-region p vector produces matching region sizes", {
  expect_length(sim_F2$genotypes, 3L)
  expect_equal(sim_F2$genotypes[[1]]$p,  50L)
  expect_equal(sim_F2$genotypes[[2]]$p,  80L)
  expect_equal(sim_F2$genotypes[[3]]$p, 100L)

  expect_true(in_range_or_na(eval_F2$abf$global$ap))
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
  expect_true(is.na(ev$finemap$global$ap))
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
  expect_true(is.na(ev$paintor$global$ap))
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

test_that("[G-finemap-real] real FINEMAP binary: pip length, AP in [0,1]", {
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
  expect_true(in_range_or_na(ev$finemap$global$ap))
})


test_that("[G-paintor-real] real PAINTOR binary: pip length, AP in [0,1]", {
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
  expect_true(in_range_or_na(ev$paintor$global$ap))
})


test_that("[G-beatrice-real] real BEATRICE: pip length, AP in [0,1]", {
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
  expect_true(in_range_or_na(ev$beatrice$global$ap))
})


test_that("[G-funmap-real] real Funmap: pip length, AP in [0,1]", {
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
  expect_true(in_range_or_na(ev$funmap$global$ap))
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


test_that("[H] cs_size summaries are well-formed (non-negative or NA)", {
  g <- eval_A$susie$global
  # cs_size_mean / cs_size_median are non-negative when defined. We do NOT
  # assert mean >= median ("right-skewed"): for the small number of credible
  # sets in this fixture that ordering is not a stable invariant and flips
  # across BLAS/LAPACK implementations (e.g. CI Linux vs. local macOS).
  expect_true(is.na(g$cs_size_mean)   || g$cs_size_mean   >= 0)
  expect_true(is.na(g$cs_size_median) || g$cs_size_median >= 0)

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


test_that("[I3] all-fits-failed method: NA AP, n_failed correct, others unaffected", {
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

  expect_true(is.na(ev$badmeth$global$ap))
  expect_equal(ev$badmeth$global$n_failed, 16L)
  expect_null(ev$badmeth$global$fdr_power_curve)

  # Other methods unaffected by the failed one
  expect_false(is.na(ev$susie$global$ap))
})


SE_SCALAR_FIELDS <- c("ap_se", "cs_coverage_se", "cs_power_se",
                      "cs_size_median_se", "cs_size_mean_se", "runtime_mean_se")


# =============================================================================
# [L] save argument - evaluate_methods writes correct files
# =============================================================================

test_that("[L1] save=TRUE creates output_dir and writes evaluation.rds + .csv", {
  out_dir <- tempfile("test_L1_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  ev <- evaluate_methods(sim_A, res_A,
                          save = TRUE, output_dir = out_dir, verbose = FALSE)

  expect_true(dir.exists(out_dir))
  expect_true(file.exists(file.path(out_dir, "evaluation.rds")))
  expect_true(file.exists(file.path(out_dir, "evaluation_summary.csv")))
})


test_that("[L4] save=FALSE does not create the output_dir", {
  tmp <- tempfile("test_L4_nosave_")
  evaluate_methods(sim_A, res_A,
                   save = FALSE, output_dir = tmp, verbose = FALSE)
  expect_false(dir.exists(tmp))
})


# =============================================================================
# [M] Per-region metrics: nothing pooled across regions
#
# The rule being tested: per-region SCALARS are averaged, COUNTS are summed and
# the rate formed afterwards. Pooling a ranking across regions lets one
# region's false positives outrank another region's true positive, which is a
# different quantity - not a noisier estimate of the same one.
# =============================================================================

test_that("[M] per_region has exactly one row per (method, scenario, region)", {
  pr <- per_region_metrics(sim_A, res_A, verbose = FALSE)
  expect_s3_class(pr, "data.frame")
  key <- pr[c("method", "scenario_id", "region_id")]
  expect_equal(nrow(pr), nrow(unique(key)))
  expect_setequal(unique(pr$method), res_A$methods_run)
  # 2 regions x (2 S x 2 phi x 2 iter) scenarios x 4 methods
  expect_equal(nrow(pr), 2L * length(sim_A$scenarios) * length(res_A$methods_run))
})

test_that("[M] per_region carries the design factors and the counts", {
  pr <- per_region_metrics(sim_A, res_A, verbose = FALSE)
  for (col in c("method", "scenario_id", "region_id", "S", "phi",
                "region_size", "n_variants", "n_causal", "ap",
                "failed", "bad_pip", "runtime_seconds",
                "n_cs", "n_cs_with_causal", "n_causal_captured",
                "nsel_at_095", "tp_at_095", "n_band_top", "c_band_top"))
    expect_true(col %in% names(pr), info = col)
  ok <- !pr$failed & !pr$bad_pip
  expect_true(all(pr$ap[ok] >= 0 & pr$ap[ok] <= 1, na.rm = TRUE))
  expect_true(all(pr$tp_at_095[ok] <= pr$nsel_at_095[ok], na.rm = TRUE))
  expect_true(all(pr$n_causal_captured[ok] <= pr$n_causal[ok], na.rm = TRUE))
})

test_that("[M] a failed fit becomes an all-NA metric row, not a dropped row", {
  res_bad <- res_A
  res_bad$susie$results[[1]]$error <- "forced failure"
  pr <- per_region_metrics(sim_A, res_bad, verbose = FALSE)
  expect_equal(nrow(pr),
               nrow(per_region_metrics(sim_A, res_A, verbose = FALSE)))
  row <- pr[pr$method == "susie" &
            pr$scenario_id == res_A$susie$results[[1]]$scenario_id &
            pr$region_id   == res_A$susie$results[[1]]$region_id, ]
  expect_true(row$failed)
  expect_true(is.na(row$ap))
})

test_that("[M] non-finite PIPs are marked failed, never repaired", {
  res_bad <- res_A
  res_bad$susie$results[[1]]$pip[1] <- NaN
  pr <- per_region_metrics(sim_A, res_bad, verbose = FALSE)
  row <- pr[pr$method == "susie" &
            pr$scenario_id == res_A$susie$results[[1]]$scenario_id &
            pr$region_id   == res_A$susie$results[[1]]$region_id, ]
  expect_true(row$bad_pip)
  expect_true(is.na(row$ap))
})

test_that("[M] evaluate_methods()$ap is the mean of the per-region APs", {
  pr <- per_region_metrics(sim_A, res_A, verbose = FALSE)
  for (m in res_A$methods_run) {
    v <- pr$ap[pr$method == m]
    v <- v[is.finite(v)]
    if (!length(v)) next
    expect_equal(eval_A[[m]]$global$ap, mean(v), tolerance = 1e-12,
                 info = sprintf("method=%s", m))
    expect_equal(eval_A[[m]]$global$ap_n_fits, length(v),
                 info = sprintf("method=%s", m))
  }
})


test_that("[M] compute_ap_exact scores ties as a block and guards non-finites", {
  # A perfect ranking.
  expect_equal(compute_ap_exact(c(0.9, 0.1, 0.05), c(TRUE, FALSE, FALSE)), 1)
  # All tied: the causal shares one block of 4, so precision is 1/4.
  expect_equal(compute_ap_exact(rep(0.25, 4), c(TRUE, FALSE, FALSE, FALSE)),
               0.25)
  # Two causals, ranked 1st and 3rd of 4: (1/1 + 2/3) / 2.
  expect_equal(compute_ap_exact(c(0.4, 0.3, 0.2, 0.1),
                                c(TRUE, FALSE, TRUE, FALSE)),
               (1 + 2 / 3) / 2)
  # No causal variants, and non-finite input.
  expect_true(is.na(compute_ap_exact(runif(10), rep(FALSE, 10))))
  expect_true(is.na(compute_ap_exact(c(1, NaN, 0.5), c(TRUE, FALSE, FALSE))))
  expect_true(is.na(compute_ap_exact(c(1, Inf, 0.5), c(TRUE, FALSE, FALSE))))
})


test_that("[M] evaluate_methods(per_region = FALSE) omits the table", {
  ev <- evaluate_methods(sim_A, res_A, per_region = FALSE, verbose = FALSE)
  expect_null(ev$per_region)
  ev2 <- evaluate_methods(sim_A, res_A, verbose = FALSE)
  expect_s3_class(ev2$per_region, "data.frame")
})


# =============================================================================
# [N] FDR and calibration are region-weighted, not call-weighted
#
# Fine-mapping is applied to a region, so the question is "how trustworthy are
# the calls in a region", not "what fraction of all calls made anywhere were
# wrong". The two differ whenever regions make unequal numbers of calls.
# =============================================================================

test_that("[N] the curve carries both weightings plus the region count", {
  d <- eval_A$susie$global$fdr_power_curve
  for (col in c("threshold", "fdr", "precision", "power", "recall",
                "n_regions", "n_regions_with_causal",
                "fdr_pooled", "precision_pooled", "power_pooled",
                "recall_pooled", "tp", "fp", "fn"))
    expect_true(col %in% names(d), info = col)
  expect_true(all(d$fdr   >= 0 & d$fdr   <= 1, na.rm = TRUE))
  expect_true(all(d$power >= 0 & d$power <= 1, na.rm = TRUE))
  # n_regions can only fall as the threshold rises.
  d2 <- d[order(d$threshold), ]
  expect_true(all(diff(d2$n_regions) <= 0))
})

test_that("[N] a region that selects nothing is excluded, not scored FDR 0", {
  d <- eval_A$susie$global$fdr_power_curve
  # Wherever no region selected anything, the region-weighted FDR is NA while
  # the pooled convention reports 0. That convention is why averaging it would
  # reward abstention.
  none <- d$n_regions == 0L
  if (any(none)) {
    expect_true(all(is.na(d$fdr[none])))
    expect_true(all(d$fdr_pooled[none] == 0))
  }
  expect_true(all(!is.na(d$fdr[d$n_regions > 0L])))
})


test_that("[N] the two weightings genuinely differ on unequal call counts", {
  # Two regions: one makes 1 correct call, the other 3 calls of which 1 is
  # correct. Region-weighted FDR = mean(0, 2/3) = 1/3. Call-weighted FDR =
  # 2 false of 4 calls = 1/2. Neither is a noisy version of the other.
  mk <- function(pip, causal) list(
    pip = pip, credible_sets = list(), method = "t", scenario_id = 1L,
    region_id = 1L, n_variants = length(pip), causal_indices = causal,
    runtime_seconds = 0)
  f1 <- mk(c(0.99, 0.01, 0.01, 0.01), 1L)            # 1 call, correct
  f2 <- mk(c(0.99, 0.99, 0.99, 0.01), 1L)            # 3 calls, 1 correct
  per_fdr <- vapply(list(f1, f2), function(f) {
    sel <- f$pip >= 0.95
    tp  <- sum(seq_along(f$pip)[sel] %in% f$causal_indices)
    (sum(sel) - tp) / sum(sel)
  }, numeric(1))
  expect_equal(mean(per_fdr), 1 / 3)
  pooled <- (sum(c(1, 3)) - sum(c(1, 1))) / sum(c(1, 3))
  expect_equal(pooled, 1 / 2)
  expect_false(isTRUE(all.equal(mean(per_fdr), pooled)))
})

test_that("[N] calibration is region-weighted and reports its region count", {
  pc <- eval_A$susie$global$pip_calibration
  for (col in c("bin", "bin_lower", "bin_upper", "mean_pip", "frac_causal",
                "n_regions", "mean_pip_pooled", "frac_causal_pooled",
                "n", "n_causal"))
    expect_true(col %in% names(pc), info = col)
  expect_true(all(pc$frac_causal >= 0 & pc$frac_causal <= 1, na.rm = TRUE))
  # A bin populated in no region is NA, never 0.
  expect_true(all(is.na(pc$frac_causal[pc$n_regions == 0L])))
  # The counts still total every variant of every valid fit.
  valid <- Filter(function(f) is.null(f$error), res_A$susie$results)
  expect_equal(sum(pc$n), sum(vapply(valid, function(f) length(f$pip),
                                     numeric(1))))
})

test_that("[N] per-region rates are NA where undefined, never 0", {
  pr <- per_region_metrics(sim_A, res_A, verbose = FALSE)
  ok <- !pr$failed & !pr$bad_pip
  # No selections at a threshold -> FDR and precision undefined for that region.
  empty <- ok & pr$nsel_at_099 == 0L
  if (any(empty)) {
    expect_true(all(is.na(pr$fdr_at_099[empty])))
    expect_true(all(is.na(pr$precision_at_099[empty])))
    # Power is still defined: zero of the causal variants were found.
    expect_true(all(pr$power_at_099[empty] == 0))
  }
  # A band with no variant in this region is NA, not 0.
  empty_band <- ok & pr$n_band_top == 0L
  if (any(empty_band))
    expect_true(all(is.na(pr$frac_causal_band_top[empty_band])))
})


# =============================================================================
# [O] The unit is one region of one iteration, and SEs need the iteration
# =============================================================================

test_that("[O] one per_region row is one (method, iteration, region)", {
  pr <- per_region_metrics(sim_A, res_A, verbose = FALSE)
  # The scenario id and the (S, phi, iter) triple identify the same thing.
  expect_equal(nrow(pr),
               nrow(unique(pr[c("method", "scenario_id", "region_id")])))
  expect_equal(nrow(pr),
               nrow(unique(pr[c("method", "S", "phi", "iter", "region_id")])))
  # Every scenario carries exactly one iteration.
  expect_true(all(tapply(pr$iter, pr$scenario_id,
                         function(v) length(unique(v))) == 1L))
  expect_false(anyNA(pr$iter))
  expect_setequal(unique(pr$iter), seq_len(2L))   # sim_A has n_iter = 2
})


test_that("[O] LD and annotations are shared across a region's iterations", {
  # This is why the iteration, not the region-fit, is the independent unit:
  # genotypes are stored per REGION and reused by every scenario, so only the
  # phenotype draw changes between iterations.
  expect_equal(length(sim_A$genotypes), 2L)          # per region, not per fit
  expect_true(length(sim_A$scenarios) > length(sim_A$genotypes))
  # The same region's truth differs between iterations, so something does vary.
  tr <- lapply(sim_A$scenarios, function(sc) sc$regions[[1]]$truth$causal_indices)
  expect_true(length(unique(vapply(tr, function(x) paste(x, collapse = ","),
                                   character(1)))) > 1L)
})


# =============================================================================
# [P] The output shape after the post-analysis layer was removed
#
# evaluate_methods() returns `global` per method plus one `per_region` table,
# and nothing else: no pre-chosen strata, no standard errors, no plots. Any
# other grouping is the caller's, from $per_region.
# =============================================================================

test_that("[P] each method carries global and nothing else", {
  for (m in eval_A$methods_evaluated) {
    expect_named(eval_A[[m]], "global", info = m)
    g <- eval_A[[m]]$global
    expect_true(has_names(g,
      c("ap", "ap_n_fits", "fdr_power_curve", "pip_calibration",
        "cs_coverage", "cs_power", "cs_size_median", "cs_size_mean",
        "cs_coverage_n", "cs_power_n", "n_cs_reported",
        "runtime_mean", "runtime_sd", "n_fits", "n_failed")), info = m)
    # The removed layers.
    for (gone in c("by_S", "by_phi", "by_p_causal", "by_region_size",
                   "by_causal_maf", "by_true_annotation_type",
                   "auprc", "ap_pooled_grid",
                   "ap_se", "cs_coverage_se", "runtime_mean_se"))
      expect_null(g[[gone]], info = paste(m, gone))
  }
  expect_null(eval_A$susie$by_S)
})

test_that("[P] the top level carries per_region and the run metadata", {
  expect_s3_class(eval_A$per_region, "data.frame")
  expect_true(all(c("methods_evaluated", "simulation_params",
                    "pip_thresholds_used") %in% names(eval_A)))
  # Every method evaluated appears in the one shared table.
  expect_setequal(unique(eval_A$per_region$method), eval_A$methods_evaluated)
})

test_that("[P] global is the per-region metrics over all of a method's fits", {
  pr <- eval_A$per_region
  for (m in eval_A$methods_evaluated) {
    v <- pr$ap[pr$method == m]
    v <- v[is.finite(v)]
    if (!length(v)) next
    expect_equal(eval_A[[m]]$global$ap, mean(v), tolerance = 1e-12, info = m)
    expect_equal(eval_A[[m]]$global$ap_n_fits, length(v), info = m)
  }
})

test_that("[P] annotated simulations evaluate and carry the annotation factor", {
  for (m in eval_B$methods_evaluated) {
    g <- eval_B[[m]]$global
    expect_true(is.numeric(g$ap))
    expect_true(is.na(g$ap) || (g$ap >= 0 && g$ap <= 1), info = m)
  }
  pr <- eval_B$per_region
  expect_true("true_annotation_type" %in% names(pr))
  expect_setequal(unique(pr$true_annotation_type), "binary")
  # And a sparse_inf run carries p_causal, so it can be grouped on.
  expect_false(all(is.na(eval_E$per_region$p_causal)))
})

test_that("[P] the pooled-ranking AP estimators are gone", {
  # auprc (trapezoid on the pooled ranking) and ap_pooled_grid (step integral
  # on the same) estimated the same quantity as `ap` but on the wrong
  # population, and disagreed with it about which method was best.
  for (m in res_A$methods_run) {
    g <- eval_A[[m]]$global
    expect_null(g$auprc)
    expect_null(g$ap_pooled_grid)
    expect_true(is.numeric(g$ap))
  }
})

test_that("[P] the region-weighted FDR is the mean of the per-region FDRs", {
  pr <- eval_A$per_region
  for (m in res_A$methods_run) {
    curve <- eval_A[[m]]$global$fdr_power_curve
    v <- pr$fdr_at_095[pr$method == m]
    v <- v[is.finite(v)]
    if (is.null(curve)) {
      # No usable fit, so there is no curve at all - and correspondingly no
      # per-region FDR to average. A method that produced nothing must not be
      # handed an FDR of 0.
      expect_length(v, 0L)
      next
    }
    i <- which(abs(curve$threshold - 0.95) < 1e-9)
    got <- curve$fdr[i]
    if (!length(v)) expect_true(is.na(got), info = m)
    else expect_equal(got, mean(v), tolerance = 1e-9, info = m)
  }
})

test_that("[P] a method whose every fit failed gets no curve, not a zero FDR", {
  sim_P <- run_simulation(n_regions = 1, n = 150, p = 50, n_iter = 1,
                          S = 1, phi = 0.3, seed = 11, verbose = FALSE)
  res_P <- run_methods(
    sim_P, methods = "finemap",
    method_args = list(finemap = list(finemap_path = "/nonexistent/finemap")),
    verbose = FALSE)
  ev <- evaluate_methods(sim_P, res_P, verbose = FALSE)
  g  <- ev$finemap$global
  expect_equal(g$n_failed, g$n_fits)
  expect_true(is.na(g$ap))
  expect_equal(g$ap_n_fits, 0L)
  expect_null(g$fdr_power_curve)
  expect_null(g$pip_calibration)
  # The row survives in the table, all-NA, so the absence stays visible.
  expect_equal(nrow(ev$per_region), 1L)
  expect_true(ev$per_region$failed[1])
  expect_true(is.na(ev$per_region$ap[1]))
})

test_that("[P] evaluation.rds and the summary CSV round-trip", {
  out_dir <- tempfile("test_P_save_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  ev <- evaluate_methods(sim_A, res_A, save = TRUE, output_dir = out_dir,
                         verbose = FALSE)
  ev2 <- readRDS(file.path(out_dir, "evaluation.rds"))
  expect_identical(ev$susie$global$ap, ev2$susie$global$ap)
  expect_s3_class(ev2$per_region, "data.frame")
  expect_equal(nrow(ev2$per_region), nrow(ev$per_region))

  csv <- read.csv(file.path(out_dir, "evaluation_summary.csv"))
  expect_setequal(csv$method, ev$methods_evaluated)
  expect_false(any(grepl("_se$", names(csv))))
  for (m in ev$methods_evaluated)
    expect_equal(csv$ap[csv$method == m], ev[[m]]$global$ap,
                 tolerance = 1e-10, info = m)
})

test_that("[P] plot_results and the aggregation helpers are gone", {
  for (fn in c("plot_results", "aggregate_region_metrics", "add_rate_metrics"))
    expect_false(exists(fn, where = asNamespace("fmbenchmark"),
                        inherits = FALSE), info = fn)
  # The metric functions that evaluate the fine-mapping output stay.
  for (fn in c("evaluate_methods", "per_region_metrics", "compute_ap_exact"))
    expect_true(exists(fn, where = asNamespace("fmbenchmark"),
                       inherits = FALSE), info = fn)
})
