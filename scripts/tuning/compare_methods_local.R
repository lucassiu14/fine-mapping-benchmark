#!/usr/bin/env Rscript
# =============================================================================
# scripts/tuning/compare_methods_local.R
#
# Run several methods on ONE scenario of a local sim and print the metrics
# side by side. Built to answer a specific question: when Functional BEATRICE
# does badly under `sparse_inf`, is that the generic penalty every sparse method
# pays under a polygenic background, or is it specific to FB?
#
# The answer on a real-1000G-LD bed (sparse_inf rho=0.2, S=3, phi=0.2, 2 regions
# of 500/1000, n=5000) was that it is partly specific: susie_inf 0.595,
# susie 0.418, fb_xregion 0.368. SuSiE is also a sparse model, but it runs with
# estimate_residual_variance = TRUE and can absorb the background into its
# residual variance, whereas FB takes sigma_sq and n as given and has no
# residual-variance adaptation at all.
#
# Usage:
#   Rscript scripts/tuning/compare_methods_local.R <sim.rds> <m1,m2,...> [prior_mass]
#
# Build a bed first with make_tuning_sim.R. Note that run_simulation samples one
# VCF per region (R/run_simulation.R:264), so n_regions cannot exceed the number
# of .vcf.gz files in data/vcf.
# =============================================================================
suppressWarnings(suppressMessages({
  if (requireNamespace("fmbenchmark", quietly = TRUE)) {
    library(fmbenchmark)
  } else if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else stop("fmbenchmark not installed and pkgload/DESCRIPTION unavailable")
}))

# `%||%` is an internal package helper and is not exported; define it locally so
# this runs whether fmbenchmark is installed or loaded via pkgload (same as
# fb_objective.R does).
`%||%` <- function(x, y) if (is.null(x)) y else x

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2L)
  stop("usage: compare_methods_local.R <sim.rds> <m1,m2,...> [prior_mass]",
       call. = FALSE)
sim <- readRDS(a[1]); meths <- strsplit(a[2], ",")[[1]]
pm <- if (length(a) >= 3) as.numeric(a[3]) else 0
margs <- list()
for (m in meths) {
  if (m %in% c("fb_xregion", "functional_beatrice", "beatrice")) {
    margs[[m]] <- list(beatrice_dir = normalizePath("BEATRICE_annot_sparse"),
                       python = "python3", max_iter = 2000L, n_caus = 10L,
                       sigma_sq = 0.0899, sparse_concrete = 51L,
                       prior_regularisation = 0.1255, lambda_l1 = 0.1223,
                       hierarchy_M = 10.15, gamma = 0.1, gamma_key = 0.2,
                       prior_mass = pm)
  }
}
res <- run_methods(sim, methods = meths, method_args = margs,
                   save = FALSE, verbose = FALSE)
ev <- evaluate_methods(sim, res, save = FALSE, verbose = FALSE)
# mass, FDR@.95 and reliability come from global$pip_calibration /
# global$fdr_power_curve, exactly as fb_objective.R derives them - `global` has
# no scalar fdr_at_95 or mass field, and reading one gives NA for every method.
metrics <- function(g) {
  fdr95 <- NA_real_; mass <- NA_real_; reliab <- NA_real_
  fpc <- g$fdr_power_curve
  if (!is.null(fpc) && nrow(fpc) > 0) {
    i <- which.min(abs(fpc$threshold - 0.95))
    if (length(i) && (fpc$tp[i] + fpc$fp[i]) > 0)
      fdr95 <- fpc$fp[i] / (fpc$tp[i] + fpc$fp[i])
  }
  pc <- g$pip_calibration
  if (!is.null(pc) && nrow(pc) > 0) {
    mass <- sum(pc$n * pc$mean_pip, na.rm = TRUE) /
            max(sum(pc$n_causal, na.rm = TRUE), 1)
    # Reliability is the TOP calibration bin ([0.9, 1.0]), taken exactly as
    # fb_objective.R takes it: the last row of pip_calibration. Selecting rows
    # by mean_pip >= 0.9 instead picks up NA rows for empty bins.
    top <- nrow(pc)
    if (isTRUE(pc$n[top] > 0)) reliab <- pc$n_causal[top] / pc$n[top]
  }
  list(fdr95 = fdr95, mass = mass, reliab = reliab)
}
cat(sprintf("%-22s %8s %8s %9s %9s\n", "method", "AP", "mass", "fdr@.95", "reliab"))
for (m in meths) {
  g <- tryCatch(ev[[m]]$global, error = function(e) NULL)
  if (is.null(g)) { cat(sprintf("%-22s %8s\n", m, "FAILED")); next }
  x <- metrics(g)
  cat(sprintf("%-22s %8.4f %8.3f %9s %9s\n", m, g$ap %||% NA_real_, x$mass,
              format(round(x$fdr95, 3)), format(round(x$reliab, 3))))
}
