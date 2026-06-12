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
methods <- c("susie", "abf", "marginal_z", "polyfun_oracle", "polyfun_est")
results <- run_methods(
  sim,
  methods = methods,
  method_args = list(
    susie          = list(L = 5, coverage = 0.95),
    abf            = list(prior_variance = 0.04),
    marginal_z     = list(coverage = 0.95),
    polyfun_oracle = list(L = 5),
    polyfun_est    = list(L = 5)
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

# Sanity: marginal_z should not beat susie or polyfun on a setting with signal.
ok <- all(c(
  evaluation$marginal_z$global$auprc <= evaluation$susie$global$auprc,
  evaluation$marginal_z$global$auprc <= evaluation$polyfun_oracle$global$auprc
))

cat(sprintf("\nTotal runtime: %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))
if (ok) {
  cat("Smoke test PASSED. Submit with: sbatch scripts/hpc/submit_benchmark.sh\n")
} else {
  cat("WARNING: AUPRC ordering looks off (marginal_z >= susie/polyfun_oracle).\n",
      "         Investigate before submitting the full array.\n")
}
