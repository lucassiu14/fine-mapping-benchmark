#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_check_piptail.R
#
# Validate the rescued PIP tails against the already-collected L3 counts BEFORE
# anything is rebuilt from them.
#
#   Rscript scripts/analysis/iter004_check_piptail.R <piptail_dir> <l3_rds>
#
# WHY THIS EXISTS. extract_pip_tail.R reads truth as
#     sc$regions[[region_id]]$truth$causal_indices
# and falls back to integer(0) if that path is empty. A silent failure there
# does not error - it produces a tail in which NOTHING is causal, and the
# calibration figure rebuilt from it would show every method sitting on the
# floor. That reads as a dramatic finding rather than as a broken join.
#
# The cross-check is exact. The frozen bands in iter004_collect.R were computed
# from the same PIPs, so for any stratum:
#
#     sum(n_band_top) over cells  ==  # tail entries with pip in [0.9, 1]
#     sum(c_band_top) over cells  ==  # of those that are causal
#
# and likewise for the hi band [0.5, 0.9). Both live above the 0.01 floor, so
# truncation cannot explain a difference. Anything more than a rounding
# discrepancy means fits are missing or the truth join failed.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
piptail_dir <- args[1]
l3_file     <- if (length(args) >= 2) args[2] else "results/iter004/combined_scenario_metrics.rds"
stopifnot(dir.exists(piptail_dir), file.exists(l3_file))

files <- sort(list.files(piptail_dir, pattern = "^piptail_.*\\.rds$", full.names = TRUE))
message("piptail files: ", length(files), " of 45 expected")
if (length(files) < 45L) message("  NOTE: incomplete; the comparison below is scaled accordingly")

L3     <- readRDS(l3_file)
design <- unique(L3[, c("job_dir", "model", "annotation_type")])

`%||%` <- function(a, b) if (is.null(a)) b else a
BANDS <- list(hi = c(0.5, 0.9), top = c(0.9, 1 + 1e-9))
obs <- new.env(parent = emptyenv())
bump <- function(k, v) obs[[k]] <- (obs[[k]] %||% c(0, 0)) + v

n_fits <- 0L; floors <- numeric(0); seen_jobs <- character(0)
for (f in files) {
  o <- readRDS(f); floors <- c(floors, o$floor); seen_jobs <- c(seen_jobs, o$job_dir)
  d <- design[match(o$job_dir, design$job_dir), ]
  for (fit in o$tail) {
    n_fits <- n_fits + 1L
    for (b in names(BANDS)) {
      e <- BANDS[[b]]; sel <- fit$pip >= e[1] & fit$pip < e[2]
      if (any(sel)) bump(paste(d$model[1], d$annotation_type[1], fit$method, b, sep = "\r"),
                         c(sum(sel), sum(fit$is_causal[sel])))
    }
  }
}
message("fits in tails: ", n_fits)
if (length(unique(floors)) != 1L) message("  WARNING: mixed floors: ", paste(unique(floors), collapse = ", "))

# Expected, from the frozen band counts, restricted to the rows actually rescued.
E <- L3[L3$job_dir %in% seen_jobs & L3$n_failed < L3$n_fits, ]
key <- function(b) paste(E$model, E$annotation_type, E$method, b, sep = "\r")
exp_tab <- do.call(rbind, lapply(names(BANDS), function(b) data.frame(
  k = key(b), n = E[[paste0("n_band_", b)]], c = E[[paste0("c_band_", b)]],
  stringsAsFactors = FALSE)))
agg <- aggregate(cbind(n, c) ~ k, exp_tab, sum, na.rm = TRUE)

agg$n_obs <- vapply(agg$k, function(k) (obs[[k]] %||% c(0, 0))[1], 0)
agg$c_obs <- vapply(agg$k, function(k) (obs[[k]] %||% c(0, 0))[2], 0)
p          <- do.call(rbind, strsplit(agg$k, "\r", fixed = TRUE))
agg$method <- p[, 3]; agg$band <- p[, 4]

n_ratio <- sum(agg$n_obs) / max(sum(agg$n), 1)
c_ratio <- sum(agg$c_obs) / max(sum(agg$c), 1)
message(sprintf("\nvariants above 0.5:  expected %d  rescued %d  (ratio %.4f)",
                sum(agg$n), sum(agg$n_obs), n_ratio))
message(sprintf("of those causal:     expected %d  rescued %d  (ratio %.4f)",
                sum(agg$c), sum(agg$c_obs), c_ratio))

fail <- character(0)
dead <- agg[agg$c > 0 & agg$c_obs == 0, ]
if (nrow(dead))
  fail <- c(fail, sprintf("%d method/band groups have ZERO rescued causal variants where L3 has some (truth join failed): %s",
                          nrow(dead), paste(unique(dead$method), collapse = ", ")))
if (!is.finite(n_ratio) || n_ratio < 0.97 || n_ratio > 1.03)
  fail <- c(fail, sprintf("variant-count ratio %.4f is outside [0.97, 1.03] - fits are missing", n_ratio))
if (!is.finite(c_ratio) || c_ratio < 0.97 || c_ratio > 1.03)
  fail <- c(fail, sprintf("causal-count ratio %.4f is outside [0.97, 1.03]", c_ratio))

off <- agg[agg$n > 200 & abs(agg$n_obs / pmax(agg$n, 1) - 1) > 0.05, ]
if (nrow(off)) {
  message("\nper-group discrepancies over 5% (n > 200):")
  print(head(off[order(-abs(off$n_obs / pmax(off$n, 1) - 1)),
                 c("method", "band", "n", "n_obs", "c", "c_obs")], 12), row.names = FALSE)
}

if (length(fail)) {
  message("\nFAILED:"); for (x in fail) message("  - ", x)
  quit(save = "no", status = 1L)
}
message("\nPASS - the rescued tails reproduce the frozen band counts. Safe to rebuild.")
