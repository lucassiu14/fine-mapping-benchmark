#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_calibration_bands.R
#
# Rebuild the PIP calibration table at NINE bands instead of the four frozen
# into iter004_collect.R (lo/mid/hi/top). Consumes the output of
# scripts/hpc/extract_pip_tail.R, which is the only surviving source of raw
# PIPs once $EPHEMERAL purges results.rds.
#
#   Rscript scripts/analysis/iter004_calibration_bands.R <piptail_dir> <l3_rds> <out_rds>
#
# Band edges were chosen 2026-09-02 to align with the five FDR decision
# thresholds, so the calibration and FDR figures can be read band-for-threshold.
#
# The sub-floor band [0, 0.01) is recoverable and IS included: extract_pip_tail
# stores n_below and sum_pip_below per fit, and retains every causal variant
# regardless of its PIP, so n, c and sum_pip are all exact down there. Set
# INCLUDE_SUB_FLOOR <- FALSE to start the curve at the floor instead.
#
# Counts are pooled to the stratum (model x annotation_type x method) before
# any rate is formed - a reliability estimate inside one cell rests on as few
# as ten variants in the top band, which is not a measurement. Strata are never
# pooled across model or annotation type.
# =============================================================================

INCLUDE_SUB_FLOOR <- TRUE
EDGES <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 1 + 1e-9)
MIN_N <- 50L                       # bands thinner than this are not reported

args <- commandArgs(trailingOnly = TRUE)
piptail_dir <- args[1]
l3_file     <- args[2]
out_file    <- if (length(args) >= 3) args[3] else "results/iter004/calibration_bands9.rds"
stopifnot(dir.exists(piptail_dir), file.exists(l3_file))

files <- sort(list.files(piptail_dir, pattern = "^piptail_.*\\.rds$", full.names = TRUE))
if (!length(files)) stop("no piptail_*.rds in ", piptail_dir, call. = FALSE)
message("piptail files: ", length(files))

# job_dir -> (model, annotation_type). 45 rows, one per parameter row.
L3     <- readRDS(l3_file)
design <- unique(L3[, c("job_dir", "model", "annotation_type")])
stopifnot(!anyDuplicated(design$job_dir))

lab <- function(i) {
  if (i == 0L) return("[0,0.01)")
  close <- if (i == length(EDGES) - 1L) "]" else ")"
  sprintf("[%g,%g%s", EDGES[i], min(EDGES[i + 1L], 1), close)
}
BAND_LABELS <- c(if (INCLUDE_SUB_FLOOR) lab(0L), vapply(seq_len(length(EDGES) - 1L), lab, ""))

acc <- new.env(parent = emptyenv())      # key -> c(n, c, sum_pip)
bump <- function(key, n, cc, sp) {
  v <- acc[[key]]
  acc[[key]] <- if (is.null(v)) c(n, cc, sp) else v + c(n, cc, sp)
}

n_fits <- 0L; n_skipped <- 0L; n_bad <- 0L
for (f in files) {
  obj <- readRDS(f)
  if (!isTRUE(all.equal(obj$floor, EDGES[1]))) {
    warning("floor in ", basename(f), " is ", obj$floor, ", not ", EDGES[1],
            " - the lowest band is not exact", call. = FALSE)
  }
  d <- design[match(obj$job_dir, design$job_dir), ]
  if (is.na(d$model[1])) { n_skipped <- n_skipped + length(obj$tail); next }
  stratum <- paste(d$model[1], d$annotation_type[1], sep = "\r")

  for (fit in obj$tail) {
    n_fits <- n_fits + 1L
    base <- paste(stratum, fit$method, sep = "\r")
    # Non-finite PIPs (L1's bad_pip fits) would poison every comparison below.
    ok <- is.finite(fit$pip)
    n_bad <- n_bad + sum(!ok)
    if (!any(ok)) next
    p <- fit$pip[ok]; y <- fit$is_causal[ok]

    if (INCLUDE_SUB_FLOOR) {
      below <- p < EDGES[1]                       # retained causal variants only
      bump(paste(base, "[0,0.01)", sep = "\r"),
           fit$n_below, sum(y[below]), fit$sum_pip_below)
    }
    for (i in seq_len(length(EDGES) - 1L)) {
      sel <- p >= EDGES[i] & p < EDGES[i + 1L]
      if (!any(sel)) next
      bump(paste(base, lab(i), sep = "\r"), sum(sel), sum(y[sel]), sum(p[sel]))
    }
  }
}
if (n_bad) message("dropped ", n_bad, " non-finite PIPs")
message("fits binned: ", n_fits, if (n_skipped) paste0("  (skipped ", n_skipped, " with no design row)") else "")

keys <- ls(acc)
parts <- do.call(rbind, strsplit(keys, "\r", fixed = TRUE))
vals  <- do.call(rbind, lapply(keys, function(k) acc[[k]]))
tab <- data.frame(
  model = parts[, 1], annotation_type = parts[, 2], method = parts[, 3],
  band  = factor(parts[, 4], levels = BAND_LABELS),
  n = vals[, 1], c = vals[, 2], sum_pip = vals[, 3],
  stringsAsFactors = FALSE)

tab$x  <- tab$sum_pip / tab$n                       # mean PIP assigned in band
tab$y  <- tab$c / tab$n                             # proportion actually causal
tab$lo <- qbeta(0.025, tab$c + 0.5, tab$n - tab$c + 0.5)   # Jeffreys
tab$hi <- qbeta(0.975, tab$c + 0.5, tab$n - tab$c + 0.5)
tab$reportable <- tab$n >= MIN_N

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(tab, out_file)
message("wrote ", out_file, "  (", nrow(tab), " rows, ",
        sum(!tab$reportable), " below the ", MIN_N, "-variant floor)")
