# =============================================================================
# test-pipeline.R
#
# End-to-end smoke test of the full benchmark pipeline:
#   simulate -> run_methods -> evaluate -> plot
#
# Translated from scripts/test_pipeline.R as part of Phase 3b-ii. The
# original script (still in scripts/) is kept for now; 3b-iii decides
# whether to move it under inst/scripts/ or delete it.
#
# Differences from the original script:
#
# - Uses a per-test tempdir for all save = TRUE output, so tests do not
#   write to the source tree.
# - Runs only the pure-R methods (susie, abf, marginal_z) end-to-end so
#   the test passes on any machine without external binaries / Python.
# - Skips on CRAN: the simulator's first call may need to download a
#   HapMap genetic map from GitHub (~1 MB per chromosome), which CRAN
#   does not permit.
# - Drops the diagnostic pretty-printed summary tables. testthat reports
#   pass/fail itself; the diagnostic output was useful only when run
#   interactively as a script.
# - Replaces `stopifnot(...)` sanity checks with `expect_*()` calls.
# =============================================================================

test_that("end-to-end pipeline runs and produces the expected structure", {
  skip_on_cran()  # sim1000G may download a genetic map on first run

  out_dir <- tempfile("pipeline_test_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # --- 1. Simulate -----------------------------------------------------------

  sim <- run_simulation(
    n_regions     = 3,
    n             = 300,
    p             = 100,
    n_iter        = 2,
    S             = c(1, 2),
    phi           = c(0.2, 0.4),
    model         = "sparse",
    annotations   = "binary",
    n_annotations = 3,
    seed          = 42,
    save          = TRUE,
    output_dir    = out_dir,
    verbose       = FALSE
  )

  expect_length(sim$genotypes, 3L)
  expect_length(sim$scenarios, 2L * 2L * 2L)  # |S| * |phi| * n_iter
  expect_false(is.null(sim$genotypes[[1]]$annotations_matrix))
  expect_false(is.null(sim$scenarios[[1]]$regions[[1]]$truth))

  # --- 2. Run methods --------------------------------------------------------
  # Only the Tier 1 pure-R methods that are always available. External
  # methods (FINEMAP, PAINTOR, BEATRICE, Funmap) get their own
  # graceful-failure test below so we do not gate this on installs.

  pure_r_methods <- c("susie", "abf", "marginal_z")

  results <- run_methods(
    simulation  = sim,
    methods     = pure_r_methods,
    method_args = list(
      susie      = list(L = 10, coverage = 0.95),
      abf        = list(prior_variance = 0.04),
      marginal_z = list(coverage = 0.95)
    ),
    save       = TRUE,
    output_dir = out_dir,
    verbose    = FALSE
  )

  for (m in pure_r_methods) {
    expect_true(m %in% names(results),
                info = sprintf("%s missing from run_methods() output", m))
    expect_equal(results[[m]]$n_failed, 0L,
                 info = sprintf("%s had unexpected fit failures", m))
  }

  # --- 3. Evaluate -----------------------------------------------------------

  eval_out <- evaluate_methods(
    simulation = sim,
    results    = results,
    save       = TRUE,
    output_dir = out_dir,
    verbose    = FALSE
  )

  for (m in pure_r_methods) {
    g <- eval_out[[m]]$global
    expect_false(is.na(g$auprc),
                 info = sprintf("%s global AUPRC is NA", m))
    expect_gte(g$auprc, 0)
    expect_lte(g$auprc, 1)
  }

  # Sanity: marginal_z (the model-free floor) should not beat susie or abf
  # on a setting where they have signal to work with.
  expect_lte(eval_out$marginal_z$global$auprc, eval_out$susie$global$auprc)
  expect_lte(eval_out$marginal_z$global$auprc, eval_out$abf$global$auprc)

  # --- 4. Plot ---------------------------------------------------------------

  pdf_path <- file.path(out_dir, "results.pdf")
  plot_results(eval_out, output_file = pdf_path, save = TRUE, verbose = FALSE)

  expect_true(file.exists(pdf_path))
  expect_gt(file.info(pdf_path)$size, 0L)
})


test_that("run_methods returns graceful-error results for missing external methods", {
  # No need to skip on CRAN here — this test never touches the network
  # or any external binary. FINEMAP is asked to run with a path that
  # cannot exist; the wrapper should return per-fit error results
  # without crashing the batch.

  sim <- run_simulation(
    n_regions   = 2,
    n           = 200,
    p           = 50,
    n_iter      = 1,
    S           = 1,
    phi         = 0.2,
    model       = "sparse",
    annotations = "none",
    vcf_files   = NULL,    # bundled sim1000G example VCF
    seed        = 1,
    save        = FALSE,
    verbose     = FALSE
  )

  results <- run_methods(
    simulation  = sim,
    methods     = "finemap",
    method_args = list(finemap = list(finemap_path = "/definitely/not/a/binary")),
    save        = FALSE,
    verbose     = FALSE
  )

  # 2 regions × 1 scenario = 2 fits, all expected to fail gracefully
  expect_equal(results$finemap$n_total, 2L)
  expect_equal(results$finemap$n_failed, 2L)
  for (fit in results$finemap$results) {
    expect_false(is.null(fit$error))
    expect_equal(length(fit$pip), 50L)
    expect_true(all(is.na(fit$pip)))
    expect_equal(fit$method, "finemap")
  }
})
