#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/smoke_test.R
#
# Quick sanity check before submitting the full array. Runs the same code path
# as run_benchmark_job.R but with tiny parameters - finishes in a minute or two
# on a laptop. Use this to catch dependency / path issues before burning
# cluster compute.
#
# Usage (from project root):
#   Rscript scripts/hpc/smoke_test.R
# =============================================================================

# Load the package
if (requireNamespace("fmbenchmark", quietly = TRUE)) {
  library(fmbenchmark)
} else if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(quiet = TRUE)
} else {
  stop("fmbenchmark is not installed. Run from project root after `renv::restore()`.")
}

cat("====================================================================\n")
cat("  fmbenchmark HPC smoke test\n")
cat("====================================================================\n\n")

t0 <- Sys.time()

cat("[1/3] Simulating (small) ... ")
sim <- run_simulation(
  n_regions       = 3,
  n               = 200,
  p               = 100,
  n_iter          = 2,
  S               = c(1, 2),
  phi             = c(0.1, 0.4),
  model           = "sparse",
  annotations     = "binary",
  n_annotations   = 2,
  enrichment      = 5,
  genetic_map_dir = "data/genetic_maps",
  seed            = 42,
  verbose         = FALSE
)
cat(sprintf("done (%.1fs)\n", as.numeric(Sys.time() - t0, units = "secs")))
cat(sprintf("       %d regions x %d scenarios\n",
            length(sim$genotypes), length(sim$scenarios)))

cat("[2/3] Running methods ... ")
t1 <- Sys.time()
methods <- c("susie", "susie_inf", "abf", "carma",
             "marginal_z", "polyfun_oracle", "polyfun_est",
             "polyfun_ldsc", "sbayesrc")
results <- run_methods(
  sim,
  methods = methods,
  method_args = list(
    susie          = list(L = 5, coverage = 0.95),
    susie_inf      = list(L = 5),
    abf            = list(prior_variance = 0.04),
    carma          = list(num.causal = 3),
    marginal_z     = list(coverage = 0.95),
    polyfun_oracle = list(L = 5),
    polyfun_est    = list(L = 5),
    polyfun_ldsc   = list(L = 5),
    sbayesrc       = list(n_iter = 200, burn_in = 100)
  ),
  verbose = FALSE
)
cat(sprintf("done (%.1fs)\n", as.numeric(Sys.time() - t1, units = "secs")))
for (m in methods) {
  cat(sprintf("       %-16s n_fits=%d  failed=%d\n",
              m, results[[m]]$n_total, results[[m]]$n_failed))
}

cat("[3/3] Evaluating ... ")
t2 <- Sys.time()
evaluation <- evaluate_methods(sim, results, verbose = FALSE)
cat(sprintf("done (%.1fs)\n", as.numeric(Sys.time() - t2, units = "secs")))

cat("\nGlobal AUPRC by method (sanity check):\n")
for (m in methods) {
  a <- evaluation[[m]]$global$auprc
  cat(sprintf("  %-16s %s\n", m,
              if (is.null(a) || is.na(a)) "  NA" else sprintf("%.3f", a)))
}

# --- §0.3 metric-availability checks -----------------------------------------
# Confirm the FDR curve and PIP calibration curve compute correctly for
# every method, not just the AUPRC scalar. If any of these are absent the
# Phase 2 calibration gate and separation-based dataset ranking cannot run.
cat("\nMetric availability (calibration + FDR + n_pip bins):\n")
metrics_ok <- TRUE
for (m in methods) {
  g <- evaluation[[m]]$global
  cal_ok <- !is.null(g$pip_calibration) && nrow(g$pip_calibration) > 0
  fdr_ok <- !is.null(g$fdr_power_curve) && nrow(g$fdr_power_curve) > 0
  n_bins <- if (cal_ok) sum(g$pip_calibration$n > 0, na.rm = TRUE) else 0L
  cat(sprintf("  %-16s calibration=%s  fdr_curve=%s  non_empty_cal_bins=%d\n",
              m, ifelse(cal_ok, "YES", "NO "),
              ifelse(fdr_ok, "YES", "NO "), n_bins))
  metrics_ok <- metrics_ok && cal_ok && fdr_ok
}

# --- plot_results end-to-end (catches plotting regressions) -----------------
tmp_pdf <- tempfile(fileext = ".pdf")
plot_ok <- tryCatch({
  plot_results(evaluation, output_file = tmp_pdf, verbose = FALSE)
  file.exists(tmp_pdf) && file.info(tmp_pdf)$size > 1000L
}, error = function(e) { cat("plot_results error:", conditionMessage(e), "\n"); FALSE })
cat(sprintf("plot_results:      %s\n", ifelse(plot_ok, "YES", "NO ")))

# --- AUPRC ordering sanity ---------------------------------------------------
ord_ok <- all(c(
  evaluation$marginal_z$global$auprc <= evaluation$susie$global$auprc,
  evaluation$marginal_z$global$auprc <= evaluation$polyfun_oracle$global$auprc,
  # polyfun_ldsc should be at least as good as marginal_z on a signal-rich
  # setting - if not, something is broken in the corrected regressor
  evaluation$marginal_z$global$auprc <= evaluation$polyfun_ldsc$global$auprc
))
cat(sprintf("AUPRC ordering:    %s\n", ifelse(ord_ok, "OK ", "OFF")))

cat(sprintf("\nTotal runtime: %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))
if (metrics_ok && plot_ok && ord_ok) {
  cat("Smoke test PASSED. Submit with: sbatch scripts/hpc/submit_benchmark.sh\n")
} else {
  cat("WARNING: smoke test issues detected - investigate before submitting.\n")
}
