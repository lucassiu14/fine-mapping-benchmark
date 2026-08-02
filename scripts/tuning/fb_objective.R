#!/usr/bin/env Rscript
# =============================================================================
# scripts/tuning/fb_objective.R
#
# Evaluate ONE hyperparameter configuration of Functional BEATRICE on ONE
# scenario of a fixed tuning simulation, and print the metrics.
#
# Called once per scenario by scripts/tuning/optuna_fb.py, which accumulates the
# per-scenario values, reports them to Optuna as intermediate results (enabling
# pruning of bad configurations after the first scenario or two), and averages
# them into the trial's objective.
#
# Metrics are computed with the SAME formulas as scripts/hpc/collect_results.R
# (which were audited in Iteration 002) so tuning optimises exactly the
# quantities the benchmark reports - not a re-implementation that could drift:
#   ap                    step-integral average precision (macro over the
#                         scenario's fits, as evaluate_methods computes it)
#   max_fdr_violation_n20 max over thresholds of (FDR(t) - (1-t)), restricted to
#                         thresholds selecting >= 20 variants (tail-noise guard)
#   total_mass_ratio      sum(PIP) / #causal        (1.0 = calibrated)
#   hi_pip_reliab         fraction of PIP>=0.9 calls that are truly causal
#
# I/O is deliberately dependency-free (no jsonlite): hyperparameters arrive as
# --key value flags, results are printed as `key=value` lines.
#
# Usage:
#   Rscript scripts/tuning/fb_objective.R --sim <tuning_sim.rds> --scenario 1 \
#       --regions 4 --beatrice_dir BEATRICE_annot_sparse --python ~/tools/py-venv-runner.sh \
#       --max_iter 2000 --lambda_l1 0.01 --prior_regularisation 1.0 --n_caus 5 \
#       --sigma_sq 0.05 --sparse_concrete 50 --hierarchy_M 10
# =============================================================================

suppressWarnings(suppressMessages({
  if (requireNamespace("fmbenchmark", quietly = TRUE)) {
    library(fmbenchmark)
  } else if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    stop("fmbenchmark not installed and pkgload/DESCRIPTION unavailable")
  }
}))

# `%||%` is an internal helper in the package; define it locally so this script
# works whether fmbenchmark is installed (not exported) or loaded via pkgload.
`%||%` <- function(x, y) if (is.null(x)) y else x

# --- parse --key value flags -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (startsWith(a, "--")) {
    key <- sub("^--", "", a)
    val <- if (i + 1L <= length(args) && !startsWith(args[i + 1L], "--")) args[i + 1L] else "TRUE"
    opt[[key]] <- val
    i <- i + 2L
  } else i <- i + 1L
}
getn <- function(k, d) if (!is.null(opt[[k]])) as.numeric(opt[[k]]) else d
geti <- function(k, d) if (!is.null(opt[[k]])) as.integer(as.numeric(opt[[k]])) else d
gets <- function(k, d) if (!is.null(opt[[k]])) opt[[k]] else d

SIM_PATH  <- gets("sim", "")
SCENARIO  <- geti("scenario", 1L)
N_REGIONS <- geti("regions", 0L)          # 0 = use all regions in the sim
if (!nzchar(SIM_PATH) || !file.exists(SIM_PATH)) {
  stop("--sim must point at an existing tuning sim .rds (see make_tuning_sim.R)")
}

sim <- readRDS(SIM_PATH)
n_sc <- length(sim$scenarios)
if (SCENARIO < 1L || SCENARIO > n_sc) {
  stop(sprintf("--scenario %d out of range (sim has %d scenarios)", SCENARIO, n_sc))
}

# --- subset to one scenario (and optionally fewer regions) -------------------
# evaluate_methods() indexes simulation$scenarios[[fit$scenario_id]] as a LIST
# INDEX, so after subsetting the single scenario must be renumbered to 1.
mini <- sim
mini$scenarios <- sim$scenarios[SCENARIO]
mini$scenarios[[1]]$scenario_id <- 1L
if (N_REGIONS > 0L && N_REGIONS < length(mini$genotypes)) {
  keep <- seq_len(N_REGIONS)
  mini$genotypes <- mini$genotypes[keep]
  mini$scenarios[[1]]$regions <- mini$scenarios[[1]]$regions[keep]
}

# --- hyperparameters under test ---------------------------------------------
fb_args <- list(
  beatrice_dir         = gets("beatrice_dir", normalizePath("BEATRICE_annot_sparse", mustWork = FALSE)),
  python               = gets("python", "python"),
  max_iter             = geti("max_iter", 2000L),
  n_caus               = geti("n_caus", 5L),
  sigma_sq             = getn("sigma_sq", 0.05),
  gamma_coverage       = getn("gamma_coverage", 0.95),
  sparse_concrete      = geti("sparse_concrete", 50L),
  prior_regularisation = getn("prior_regularisation", 1.0),
  lambda_l1            = getn("lambda_l1", 0.01),
  hierarchy_M          = getn("hierarchy_M", 10.0)
)

t0 <- Sys.time()
res <- tryCatch(
  run_methods(mini, methods = "functional_beatrice",
              method_args = list(functional_beatrice = fb_args),
              save = FALSE, verbose = FALSE),
  error = function(e) structure(conditionMessage(e), class = "fb_fail")
)
if (inherits(res, "fb_fail")) {
  cat(sprintf("status=error\nmessage=%s\n", gsub("[\r\n]+", " ", as.character(res))))
  quit(status = 0)
}

ev <- tryCatch(
  evaluate_methods(mini, res, save = FALSE, verbose = FALSE),
  error = function(e) NULL
)
g <- if (!is.null(ev)) ev[["functional_beatrice"]]$global else NULL
if (is.null(g)) {
  cat("status=error\nmessage=evaluate_methods returned no global metrics\n")
  quit(status = 0)
}

# --- metrics: identical formulas to collect_results.R ------------------------
ap <- g$ap %||% NA_real_

max_fdr_violation_n20 <- NA_real_
fpc <- g$fdr_power_curve
if (!is.null(fpc) && nrow(fpc) > 0) {
  viol <- pmax(0, fpc$fdr - (1 - fpc$threshold))
  nsel <- fpc$tp + fpc$fp
  ok <- is.finite(viol) & nsel >= 20L
  if (any(ok)) max_fdr_violation_n20 <- max(viol[ok])
}

total_mass_ratio <- NA_real_
hi_pip_reliab    <- NA_real_
ece_hi           <- NA_real_
pc <- g$pip_calibration
if (!is.null(pc) && nrow(pc) > 0) {
  sp <- pc$n * pc$mean_pip                       # sum of PIPs in the bin
  denom <- max(sum(pc$n_causal, na.rm = TRUE), 1)
  total_mass_ratio <- sum(sp, na.rm = TRUE) / denom
  top <- nrow(pc)                                # highest PIP bin ([0.9, 1.0])
  if (isTRUE(pc$n[top] > 0)) hi_pip_reliab <- pc$n_causal[top] / pc$n[top]

  # ece_hi: expected calibration error over the INFORMATIVE range (PIP >= 0.1).
  # Formula copied from scripts/hpc/collect_results.R so the tuning objective
  # and the benchmark report the same number. Plain ECE is n-weighted and ~94%
  # of variant-observations sit in bin 1, so it mostly measures how well
  # near-zero PIPs are calibrated and can rank a badly-calibrated method well.
  keep <- pc$n > 0
  if (any(keep)) {
    mp  <- sp[keep] / pc$n[keep]                 # mean predicted PIP per bin
    fc  <- pc$n_causal[keep] / pc$n[keep]        # observed frequency per bin
    idx <- which(keep)
    hi  <- idx > 1L                              # bin 1 is [0, 0.1)
    if (any(hi) && sum(pc$n[idx][hi]) > 0) {
      wh <- pc$n[idx][hi] / sum(pc$n[idx][hi])
      ece_hi <- sum(wh * abs(mp[hi] - fc[hi]))
    }
  }
}

n_failed <- res[["functional_beatrice"]]$n_failed %||% NA_integer_
elapsed  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf(
  "status=ok\nap=%.6f\nmax_fdr_violation_n20=%s\ntotal_mass_ratio=%s\nhi_pip_reliab=%s\nece_hi=%s\nn_failed=%s\nseconds=%.1f\n",
  ap,
  ifelse(is.na(max_fdr_violation_n20), "NA", sprintf("%.6f", max_fdr_violation_n20)),
  ifelse(is.na(total_mass_ratio),      "NA", sprintf("%.6f", total_mass_ratio)),
  ifelse(is.na(hi_pip_reliab),         "NA", sprintf("%.6f", hi_pip_reliab)),
  ifelse(is.na(ece_hi),                "NA", sprintf("%.6f", ece_hi)),
  ifelse(is.na(n_failed),              "NA", as.character(n_failed)),
  elapsed))
