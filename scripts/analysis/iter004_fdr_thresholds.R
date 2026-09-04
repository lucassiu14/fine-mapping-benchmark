#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_fdr_thresholds.R
#
# Recompute threshold-based selection counts on a FINER grid than the five
# frozen into iter004_collect.R.
#
#   Rscript scripts/analysis/iter004_fdr_thresholds.R <piptail_dir> <l1_rds> <out_rds>
#
# THRESH in the collector is c(0.50, 0.80, 0.90, 0.95, 0.99), fixed at
# collection time. The rescued PIP tails make any threshold at or above the
# 0.01 floor exactly recoverable, since every PIP >= 0.01 is retained.
#
# The grid used here is the set of calibration band edges at or above 0.05, so
# the calibration and false discovery figures read against each other
# band-for-threshold - which was the reason those band edges were chosen.
#
# Counts are accumulated to the CELL (job_dir x S x phi x region_size x method)
# and the rate is formed there, matching how the collector forms fdr_at_*.
# Cells are keyed by joining each tail row to L1 on
# (job_dir, scenario_id, region_id, method) rather than by reconstructing the
# scenario expansion order, which is not recorded and would be a guess.
# =============================================================================

THRESH <- c(0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.99)
FLOOR  <- 0.01

args <- commandArgs(trailingOnly = TRUE)
piptail_dir <- args[1]
l1_file     <- if (length(args) >= 2) args[2] else "results/iter004/combined_fit_metrics.rds"
out_file    <- if (length(args) >= 3) args[3] else "results/iter004/fdr_thresholds.rds"
stopifnot(dir.exists(piptail_dir), file.exists(l1_file), min(THRESH) >= FLOOR)

files <- sort(list.files(piptail_dir, pattern = "^piptail_.*\\.rds$", full.names = TRUE))
message("piptail files: ", length(files))

message("reading L1 for the cell keys ...")
L1 <- readRDS(l1_file)
# See iter004_calibration_bands.R: the stratum is whatever the design varies.
STRATUM <- intersect(c("model", "annotation_type", "relationship"), names(L1))
message("stratum keyed on: ", paste(STRATUM, collapse = ", "))
KEEP <- c("job_dir", "scenario_id", "region_id", "method",
          "S", "phi", "region_size", STRATUM)
L1 <- L1[, KEEP]
L1$fitkey <- do.call(paste, c(L1[c("job_dir", "scenario_id", "region_id", "method")], sep = "\r"))

acc <- new.env(parent = emptyenv())
n_fits <- 0L; n_bad <- 0L; n_unmatched <- 0L

for (f in files) {
  o <- readRDS(f)
  sub <- L1[L1$job_dir == o$job_dir, ]
  if (!nrow(sub)) { n_unmatched <- n_unmatched + length(o$tail); next }
  idx <- setNames(seq_len(nrow(sub)), sub$fitkey)

  for (fit in o$tail) {
    n_fits <- n_fits + 1L
    ok <- is.finite(fit$pip)
    n_bad <- n_bad + sum(!ok)
    if (!any(ok)) next
    p <- fit$pip[ok]; y <- fit$is_causal[ok]

    i <- idx[[paste(o$job_dir, fit$scenario_id, fit$region_id, fit$method, sep = "\r")]]
    if (is.null(i)) { n_unmatched <- n_unmatched + 1L; next }
    r <- sub[i, ]
    cell <- paste(paste(unlist(r[STRATUM]), collapse = "\r"), r$method, r$job_dir,
                  r$S, r$phi, r$region_size, sep = "\r")

    v <- acc[[cell]]
    if (is.null(v)) v <- matrix(0, nrow = length(THRESH), ncol = 3)
    for (k in seq_along(THRESH)) {
      sel <- p >= THRESH[k]
      if (!any(sel)) next
      v[k, ] <- v[k, ] + c(sum(sel), sum(y[sel]), sum(p[sel]))
    }
    acc[[cell]] <- v
  }
}
message("fits: ", n_fits,
        if (n_bad) paste0("   (dropped ", n_bad, " non-finite PIPs)") else "",
        if (n_unmatched) paste0("   (", n_unmatched, " unmatched to L1)") else "")
if (!length(ls(acc))) stop("nothing accumulated", call. = FALSE)

keys  <- ls(acc)
parts <- do.call(rbind, strsplit(keys, "\r", fixed = TRUE))
nS <- length(STRATUM)
out <- do.call(rbind, lapply(seq_along(keys), function(j) {
  v <- acc[[keys[j]]]
  data.frame(setNames(as.list(parts[j, seq_len(nS)]), STRATUM),
             method = parts[j, nS + 1L],
             job_dir = parts[j, nS + 2L], S = as.numeric(parts[j, nS + 3L]),
             phi = as.numeric(parts[j, nS + 4L]),
             region_size = as.numeric(parts[j, nS + 5L]),
             t = THRESH, nsel = v[, 1], tp = v[, 2], sum_pip_sel = v[, 3],
             stringsAsFactors = FALSE)
}))
out$fdr <- ifelse(out$nsel > 0, (out$nsel - out$tp) / out$nsel, NA_real_)
saveRDS(out, out_file)
message("wrote ", out_file, "  (", nrow(out), " cell x threshold rows, ",
        length(unique(paste(out$model, out$annotation_type))), " strata)")
