#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/phase1_analysis.R
#
# Iteration 001 (Phase 1) analysis. Reads the collect_results.R aggregates and
# produces:
#   1. AP vs trapezoid-AUPRC reordering diff
#   2. PIP calibration suite   (ECE, MCE, signed bias, slope, total-mass ratio,
#                               high-PIP reliability)
#   3. FDR control suite       (observed FDR vs the <=1-t bound, max violation,
#                               power at fixed FDR, pAUC)
#   4. Best-in-stratum leaderboards (never a single aggregate ranking)
#
# All pooling sums COUNTS then computes rates - averaging per-scenario rates
# would be badly biased where denominators are tiny.
#
# Usage:  Rscript scripts/analysis/phase1_analysis.R [results_dir]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
RES  <- if (length(args) >= 1) args[1] else "results"

ev  <- readRDS(file.path(RES, "combined_evaluation.rds"))
cal <- readRDS(file.path(RES, "combined_pip_calibration.rds"))
fdr <- readRDS(file.path(RES, "combined_fdr_curves.rds"))

# Annotation arm from the job label (anNone vs anBinary).
arm_of <- function(job_dir) ifelse(grepl("anNone", job_dir), "none", "binary")
cal$arm <- arm_of(cal$job_dir); fdr$arm <- arm_of(fdr$job_dir)

# =============================================================================
# 1. Calibration suite  (pool bin COUNTS over `by`, then compute)
# =============================================================================
calib_suite <- function(cal, by = NULL) {
  keys <- c("method", by)
  g <- split(cal, lapply(keys, function(k) cal[[k]]), drop = TRUE)
  do.call(rbind, lapply(g, function(d) {
    b <- aggregate(cbind(n, n_causal, sum_pip) ~ bin, data = d, FUN = sum)
    b <- b[b$n > 0, ]
    if (!nrow(b)) return(NULL)
    mp <- b$sum_pip / b$n          # pooled mean predicted PIP per bin
    fc <- b$n_causal / b$n         # pooled observed frequency per bin
    w  <- b$n / sum(b$n)
    # weighted least squares of observed on predicted
    fit <- tryCatch(stats::lm(fc ~ mp, weights = b$n), error = function(e) NULL)
    top <- b[b$bin == max(b$bin), ]
    data.frame(
      as.list(setNames(lapply(keys, function(k) d[[k]][1]), keys)),
      ece            = sum(w * abs(mp - fc)),
      mce            = max(abs(mp - fc)),
      signed_bias    = sum(w * (mp - fc)),   # >0 = OVER-confident
      slope          = if (is.null(fit)) NA_real_ else unname(coef(fit)[2]),
      intercept      = if (is.null(fit)) NA_real_ else unname(coef(fit)[1]),
      total_mass_ratio = sum(b$n * mp) / max(sum(b$n_causal), 1),
      hi_pip_reliab  = if (nrow(top)) top$n_causal / top$n else NA_real_,
      hi_pip_n       = if (nrow(top)) top$n else 0L,
      bins_occupied  = nrow(b),
      stringsAsFactors = FALSE, row.names = NULL)
  }))
}

# =============================================================================
# 2. FDR suite  (pool tp/fp/fn over `by`, then compute)
# =============================================================================
fdr_suite <- function(fdr, by = NULL, report_t = c(0.5, 0.7, 0.9, 0.95, 0.99)) {
  keys <- c("method", by)
  g <- split(fdr, lapply(keys, function(k) fdr[[k]]), drop = TRUE)
  do.call(rbind, lapply(g, function(d) {
    a <- aggregate(cbind(tp, fp, fn) ~ threshold, data = d, FUN = sum)
    a <- a[order(a$threshold), ]
    nsel <- a$tp + a$fp
    a$fdr    <- ifelse(nsel > 0, a$fp / nsel, 0)
    a$power  <- ifelse((a$tp + a$fn) > 0, a$tp / (a$tp + a$fn), 0)
    # A calibrated method selecting PIP >= t should have FDR <= 1 - t.
    a$bound  <- 1 - a$threshold
    viol     <- pmax(0, a$fdr - a$bound)
    at <- function(t) { i <- which.min(abs(a$threshold - t)); a[i, ] }
    pw_at_fdr <- function(target) {
      ok <- a[a$fdr <= target & nsel[match(a$threshold, a$threshold)] > 0, ]
      if (!nrow(ok)) return(NA_real_)
      max(ok$power)
    }
    # partial AUC of power vs FDR, up to FDR = 0.10
    sub <- a[a$fdr <= 0.10, ]
    pauc <- if (nrow(sub) > 1) {
      o <- order(sub$fdr); sum(diff(sub$fdr[o]) * head(sub$power[o], -1))
    } else NA_real_
    row <- data.frame(
      as.list(setNames(lapply(keys, function(k) d[[k]][1]), keys)),
      max_fdr_violation = max(viol),
      pauc_fdr10        = pauc,
      pw_at_fdr05       = pw_at_fdr(0.05),
      pw_at_fdr10       = pw_at_fdr(0.10),
      pw_at_fdr20       = pw_at_fdr(0.20),
      stringsAsFactors = FALSE, row.names = NULL)
    for (t in report_t) {
      r <- at(t)
      row[[paste0("fdr_at_", t)]]   <- r$fdr
      row[[paste0("power_at_", t)]] <- r$power
    }
    row
  }))
}

cat("### 1. PIP CALIBRATION (pooled over everything, per method)\n")
cs <- calib_suite(cal)
cs <- cs[order(cs$ece), ]
print(format(cs[, c("method","ece","mce","signed_bias","slope",
                    "total_mass_ratio","hi_pip_reliab","hi_pip_n")],
             digits = 3), row.names = FALSE)

cat("\n### 2. FDR CONTROL (pooled, per method)\n")
fs <- fdr_suite(fdr)
fs <- fs[order(fs$max_fdr_violation), ]
print(format(fs[, c("method","max_fdr_violation","fdr_at_0.95","power_at_0.95",
                    "pw_at_fdr05","pw_at_fdr10","pauc_fdr10")],
             digits = 3), row.names = FALSE)

saveRDS(list(calibration = cs, fdr = fs), file.path(RES, "phase1_gates.rds"))
cat("\nwrote", file.path(RES, "phase1_gates.rds"), "\n")
