#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_calibration_bands.R
#
# Rebuild PIP calibration at ten bands, at CELL level, from the rescued PIP
# tails.
#
#   Rscript scripts/analysis/iter004_calibration_bands.R <piptail_dir> <l1_rds> <out_rds>
#
# ESTIMATOR. The reliability of a band is formed within a cell (one parameter
# row x S x phi x region size x method, 20 fits) and then averaged over cells,
# with a standard error across cells. This matches how the false discovery rate
# is formed in the results section, and the denominators support it: the top
# band holds a median of 27 variants per cell, which is the same denominator the
# FDR has at t = 0.9 - necessarily, since {PIP >= 0.9} IS the top band.
#
# The pooled estimator, sum(c)/sum(n) over all cells in the stratum, is emitted
# alongside it rather than discarded. The two diverge for exactly the methods
# that are worst calibrated, because pooling weights a cell by how many variants
# the method put in the band and those are the cells where it is wrong; that
# divergence is reported rather than hidden.
#
# Cells holding no variants in a band have no rate and are excluded from the
# cell mean. The count is carried in cells_dropped: a method is not penalised
# for declining to populate a band, exactly as an abstainer is not penalised on
# FDR, and the reader is owed the number.
#
# Band edges are the FDR thresholds, so the two figures read band for threshold.
# The sub-floor band [0, 0.01) is exact: extract_pip_tail stores n_below and
# sum_pip_below per fit and retains every causal variant whatever its PIP.
# =============================================================================

INCLUDE_SUB_FLOOR <- TRUE
EDGES <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 1 + 1e-9)
MIN_N <- 50L        # a stratum-level band thinner than this is not reported
# ...nor one resting on fewer cells than MIN_CELLS. This has to track the
# design: Iteration 004 has 500 cells per stratum and 20 is a sensible floor,
# but Iteration 005 has only 4, where a floor of 20 rejects everything. Passed
# as the fourth argument so the caller states it rather than inheriting it.
MIN_CELLS <- 20L

args <- commandArgs(trailingOnly = TRUE)
piptail_dir <- args[1]
l1_file     <- if (length(args) >= 2) args[2] else "results/iter004/combined_fit_metrics.rds"
out_file    <- if (length(args) >= 3) args[3] else "results/iter004/calibration_bands9.rds"
if (length(args) >= 4) MIN_CELLS <- as.integer(args[4])
# Band edges are an argument because the grid a design can support is a
# property of the design. Iteration 004 supports the ten-band split at 97% or
# better in every band; Iteration 005, with 25 iterations over four cells, puts
# a median of 49 variants in [0.9,0.95) and can report it in only 41% of
# method-by-relationship combinations. Merging the sparse upper bands there
# keeps the edges a SUBSET of the ten-band set, so the two remain readable
# against each other.
if (length(args) >= 5 && nzchar(args[5]))
  EDGES <- c(as.numeric(strsplit(args[5], ",", fixed = TRUE)[[1]]), 1 + 1e-9)
stopifnot(dir.exists(piptail_dir), file.exists(l1_file), is.finite(MIN_CELLS))
message("reporting floor: >= ", MIN_N, " variants and >= ", MIN_CELLS, " cells")

files <- sort(list.files(piptail_dir, pattern = "^piptail_.*\\.rds$", full.names = TRUE))
if (!length(files)) stop("no piptail_*.rds in ", piptail_dir, call. = FALSE)
message("piptail files: ", length(files))

lab <- function(i) {
  if (i == 0L) return("[0,0.01)")
  close <- if (i == length(EDGES) - 1L) "]" else ")"
  sprintf("[%g,%g%s", EDGES[i], min(EDGES[i + 1L], 1), close)
}
BAND_LABELS <- c(if (INCLUDE_SUB_FLOOR) lab(0L), vapply(seq_len(length(EDGES) - 1L), lab, ""))

message("reading L1 for the cell keys ...")
L1 <- readRDS(l1_file)
# The stratum is whatever the design actually varies. Iteration 004 varies model
# and annotation type; Iteration 005 is entirely sparse and varies the
# annotation-to-causality relationship instead, so `relationship` joins the key
# wherever the column exists. Omitting it there would collapse every
# relationship into one stratum and silently average them together.
STRATUM <- intersect(c("model", "annotation_type", "relationship"), names(L1))
message("stratum keyed on: ", paste(STRATUM, collapse = ", "))
L1 <- L1[, c("job_dir", "scenario_id", "region_id", "method",
             "S", "phi", "region_size", STRATUM)]
L1$fitkey <- do.call(paste, c(L1[c("job_dir","scenario_id","region_id","method")], sep = "\r"))

acc <- new.env(parent = emptyenv())          # cell x band -> c(n, c, sum_pip)
n_fits <- 0L; n_bad <- 0L; n_unmatched <- 0L

for (f in files) {
  o <- readRDS(f)
  if (!isTRUE(all.equal(o$floor, EDGES[1])))
    warning("floor in ", basename(f), " is ", o$floor, ", not ", EDGES[1], call. = FALSE)
  sub <- L1[L1$job_dir == o$job_dir, ]
  if (!nrow(sub)) { n_unmatched <- n_unmatched + length(o$tail); next }
  idx <- setNames(seq_len(nrow(sub)), sub$fitkey)

  for (fit in o$tail) {
    n_fits <- n_fits + 1L
    i <- idx[[paste(o$job_dir, fit$scenario_id, fit$region_id, fit$method, sep = "\r")]]
    if (is.null(i)) { n_unmatched <- n_unmatched + 1L; next }
    r <- sub[i, ]
    cell <- paste(paste(unlist(r[STRATUM]), collapse = "\r"), r$method,
                  r$job_dir, r$S, r$phi, r$region_size, sep = "\r")

    ok <- is.finite(fit$pip)                 # the bad_pip fits
    n_bad <- n_bad + sum(!ok)
    p <- fit$pip[ok]; y <- fit$is_causal[ok]

    v <- acc[[cell]]
    if (is.null(v)) v <- matrix(0, nrow = length(BAND_LABELS), ncol = 3)
    off <- 0L
    if (INCLUDE_SUB_FLOOR) {
      below <- p < EDGES[1]                  # retained causal variants only
      v[1, ] <- v[1, ] + c(fit$n_below, sum(y[below]), fit$sum_pip_below)
      off <- 1L
    }
    for (k in seq_len(length(EDGES) - 1L)) {
      sel <- p >= EDGES[k] & p < EDGES[k + 1L]
      if (any(sel)) v[k + off, ] <- v[k + off, ] + c(sum(sel), sum(y[sel]), sum(p[sel]))
    }
    acc[[cell]] <- v
  }
}
message("fits: ", n_fits,
        if (n_bad) paste0("   (dropped ", n_bad, " non-finite PIPs)") else "",
        if (n_unmatched) paste0("   (", n_unmatched, " unmatched to L1)") else "")

keys  <- ls(acc)
parts <- do.call(rbind, strsplit(keys, "\r", fixed = TRUE))
nS <- length(STRATUM)
stratum <- apply(parts[, seq_len(nS + 1L), drop = FALSE], 1, paste, collapse = "\r")

out <- do.call(rbind, lapply(split(seq_along(keys), stratum), function(ix) {
  pr <- parts[ix[1], ]
  do.call(rbind, lapply(seq_along(BAND_LABELS), function(b) {
    n  <- vapply(ix, function(j) acc[[keys[j]]][b, 1], 0)
    c_ <- vapply(ix, function(j) acc[[keys[j]]][b, 2], 0)
    sp <- vapply(ix, function(j) acc[[keys[j]]][b, 3], 0)
    k  <- n > 0
    if (!any(k)) return(NULL)
    y_cell <- c_[k] / n[k]; x_cell <- sp[k] / n[k]
    data.frame(
      setNames(as.list(pr[seq_len(nS)]), STRATUM),
      method = pr[nS + 1L],
      band = BAND_LABELS[b],
      # cell-mean estimator, and its SE across cells
      x = mean(x_cell), y = mean(y_cell),
      se = if (sum(k) > 1L) sd(y_cell) / sqrt(sum(k)) else NA_real_,
      # pooled estimator, retained for the appendix comparison
      x_pooled = sum(sp) / sum(n), y_pooled = sum(c_) / sum(n),
      n_total = sum(n), cells_used = sum(k), cells_dropped = sum(!k),
      stringsAsFactors = FALSE)
  }))
}))
out$band <- factor(out$band, levels = BAND_LABELS)
out$lo <- out$y - 2 * out$se
out$hi <- out$y + 2 * out$se
out$reportable <- out$n_total >= MIN_N & out$cells_used >= MIN_CELLS

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, out_file)
message("wrote ", out_file, "  (", nrow(out), " rows, ",
        sum(!out$reportable), " below the reporting floor)")
message(sprintf("max |cell-mean - pooled| = %.3f", max(abs(out$y - out$y_pooled), na.rm = TRUE)))
