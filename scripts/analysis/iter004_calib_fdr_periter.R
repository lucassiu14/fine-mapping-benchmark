#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_calib_fdr_periter.R
#
# PIP calibration and false discovery rate under three estimators, for the
# seven methods displayed in the LSR results section, from the rescued PIP
# tails. No cluster job is needed: every PIP >= 0.01 is retained per fit with
# its causal flag, every causal variant is retained whatever its PIP, and the
# sub-floor count is stored exactly, so all ten bands and all eight thresholds
# are recomputable at any grouping.
#
#   Rscript scripts/analysis/iter004_calib_fdr_periter.R <piptail_dir> <l1_rds> <out_rds>
#
# ESTIMATORS. Each forms a rate within a unit, then averages the unit rates
# within a stratum (model x annotation type x method), with a standard error
# across units. They differ only in the unit:
#
#   cell   job_dir x S x phi x region_size    20 fits  (the current estimator)
#   draw   job_dir x scenario_id              10 fits  one iteration of one
#          scenario: a single (S, phi, iter) draw holding its ten regions
#   fit    job_dir x scenario_id x region_id   1 fit
#
# "draw" is the per-iteration estimator requested in supervisor feedback on the
# LSR: compute the rate for each iteration of each scenario, then average.
# Before it is used, scenario_id is verified to index (S, phi, iter)
# one-to-one within every job_dir.
#
# Never pooled across model or annotation type.
#
# A unit with no variants in a band, or no selections at a threshold, has no
# rate and is excluded rather than scored; units_dropped carries the count.
# Smaller units drop more, so draw and fit condition more heavily than cell on
# the method having produced a number at all.
#
# VALIDATION. The cell estimator is rebuilt here from scratch and checked
# against results/iter004/calibration_bands9.rds and fdr_thresholds.rds. If it
# does not reproduce them the accumulation is wrong, the draw and fit numbers
# cannot be trusted either, and the script stops without writing.
# =============================================================================

M7 <- c("susie", "polyfun_ldsc", "finemap_inf", "polyfun_oracle",
        "sbayesrc", "beatrice", "fb_xregion")
EDGES  <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 1 + 1e-9)
NB     <- length(EDGES) - 1L                    # nine bands at or above the floor
THRESH <- EDGES[2:NB]                           # 0.05 ... 0.99, as fdr_thresholds
BAND_LABELS <- c("[0,0.01)",
                 sprintf("[%g,%g%s", EDGES[1:NB], pmin(EDGES[2:(NB + 1L)], 1),
                         c(rep(")", NB - 1L), "]")))
MIN_N <- 50L; MIN_UNITS <- 20L                  # same reporting rule as the cell figures

args <- commandArgs(trailingOnly = TRUE)
piptail_dir <- if (length(args) >= 1) args[1] else "results/iter004/piptail"
l1_file     <- if (length(args) >= 2) args[2] else "results/iter004/combined_fit_metrics.rds"
out_file    <- if (length(args) >= 3) args[3] else "results/iter004/calib_fdr_periter.rds"
REF_CAL <- "results/iter004/calibration_bands9.rds"
REF_FDR <- "results/iter004/fdr_thresholds.rds"
stopifnot(dir.exists(piptail_dir), file.exists(l1_file),
          file.exists(REF_CAL), file.exists(REF_FDR))

message("reading L1 keys ...")
L1 <- readRDS(l1_file)
L1 <- L1[L1$method %in% M7,
         c("job_dir", "scenario_id", "region_id", "method", "S", "phi", "iter",
           "region_size", "model", "annotation_type")]
u <- unique(L1[, c("job_dir", "scenario_id", "S", "phi", "iter")])
if (any(duplicated(u[, c("job_dir", "scenario_id")])))
  stop("scenario_id does not index a single (S, phi, iter) within job_dir", call. = FALSE)
if (any(duplicated(u[, c("job_dir", "S", "phi", "iter")])))
  stop("an (S, phi, iter) maps to more than one scenario_id", call. = FALSE)
message("verified: scenario_id <-> (S, phi, iter) is one-to-one in all ",
        length(unique(u$job_dir)), " rows")
L1$fitkey <- paste(L1$job_dir, L1$scenario_id, L1$region_id, L1$method, sep = "\r")

files <- sort(list.files(piptail_dir, pattern = "^piptail_.*\\.rds$", full.names = TRUE))
message("piptail files: ", length(files))
revcum <- function(M) { for (k in (ncol(M) - 1L):1L) M[, k] <- M[, k] + M[, k + 1L]; M }

cellL <- drawL <- fitL <- list()
n_fit <- 0L; n_bad <- 0L; n_unm <- 0L

for (f in files) {
  o  <- readRDS(f)
  tl <- o$tail
  tl <- tl[vapply(tl, function(z) z$method %in% M7, NA)]
  if (!length(tl)) next
  meth <- vapply(tl, function(z) z$method, "")
  sc   <- vapply(tl, function(z) as.integer(z$scenario_id), 0L)
  rg   <- vapply(tl, function(z) as.integer(z$region_id), 0L)
  mi   <- match(paste(o$job_dir, sc, rg, meth, sep = "\r"), L1$fitkey)
  n_unm <- n_unm + sum(is.na(mi))
  keep <- !is.na(mi); tl <- tl[keep]; mi <- mi[keep]
  nf <- length(tl); n_fit <- n_fit + nf
  if (!nf) next

  nbw <- vapply(tl, function(z) if (is.null(z$n_below)) NA_real_ else as.numeric(z$n_below), 0)
  spb <- vapply(tl, function(z) if (is.null(z$sum_pip_below)) NA_real_ else as.numeric(z$sum_pip_below), 0)
  if (anyNA(nbw) || anyNA(spb)) stop("missing sub-floor counts in ", basename(f), call. = FALSE)
  pl  <- lapply(tl, function(z) z$pip)
  yl  <- lapply(tl, function(z) z$is_causal)
  if (!identical(lengths(pl), lengths(yl))) stop("pip/is_causal length mismatch in ", basename(f), call. = FALSE)
  fid <- rep.int(seq_len(nf), lengths(pl))
  p <- unlist(pl, use.names = FALSE); y <- unlist(yl, use.names = FALSE)
  if (anyNA(y)) stop("missing causal flags in ", basename(f), call. = FALSE)
  ok <- is.finite(p); n_bad <- n_bad + sum(!ok)
  p <- p[ok]; y <- y[ok]; fid <- fid[ok]

  # sub-floor band, exactly as iter004_calibration_bands.R forms it
  c0 <- tabulate(fid[p < EDGES[1] & y == 1], nbins = nf)
  # nine above-floor bands
  b  <- findInterval(p, EDGES); a <- b >= 1L
  Nv <- Cv <- Pv <- numeric(nf * NB)
  if (any(a)) {
    ix <- (fid[a] - 1L) * NB + b[a]
    rs <- rowsum(cbind(1, y[a], p[a]), ix)
    pos <- as.integer(rownames(rs))
    Nv[pos] <- rs[, 1]; Cv[pos] <- rs[, 2]; Pv[pos] <- rs[, 3]
  }
  Nn <- matrix(Nv, nf, NB, byrow = TRUE)
  Cn <- matrix(Cv, nf, NB, byrow = TRUE)
  Pn <- matrix(Pv, nf, NB, byrow = TRUE)
  N <- cbind(nbw, Nn); C <- cbind(c0, Cn); P <- cbind(spb, Pn)       # 10 bands
  # {PIP >= EDGES[j]} is the union of bands j..NB, so selection counts at every
  # threshold follow from the band counts exactly.
  NS <- revcum(Nn)[, 2:NB, drop = FALSE]
  TP <- revcum(Cn)[, 2:NB, drop = FALSE]
  SS <- revcum(Pn)[, 2:NB, drop = FALSE]

  r <- L1[mi, ]
  strat <- paste(r$model, r$annotation_type, r$method, sep = "\r")
  X <- cbind(N, C, P, NS, TP, SS)                                    # 30 + 24 columns
  cellL[[f]] <- rowsum(X, paste(strat, r$job_dir, r$S, r$phi, r$region_size, sep = "\r"))
  drawL[[f]] <- rowsum(X, paste(strat, r$job_dir, r$scenario_id, sep = "\r"))

  # per-fit estimator: rates formed per fit, summarised by stratum here
  kc <- N > 0; yr <- ifelse(kc, C / N, 0); xr <- ifelse(kc, P / N, 0)
  kf <- NS > 0; fr <- ifelse(kf, (NS - TP) / NS, 0)
  fitL[[f]] <- rowsum(cbind(kc, yr, yr^2, xr, N, kf, fr, fr^2), strat)
}
message(sprintf("fits read: %d   non-finite PIPs dropped: %d   unmatched to L1: %d",
                n_fit, n_bad, n_unm))

split_key <- function(k) do.call(rbind, strsplit(k, "\r", fixed = TRUE))

unit_estimates <- function(U, unit) {
  kp <- split_key(rownames(U))
  s  <- paste(kp[, 1], kp[, 2], kp[, 3], sep = "\r")
  N <- U[, 1:10, drop = FALSE]; C <- U[, 11:20, drop = FALSE]; P <- U[, 21:30, drop = FALSE]
  NS <- U[, 31:38, drop = FALSE]; TP <- U[, 39:46, drop = FALSE]
  cal <- fdr <- list()
  for (g in unique(s)) {
    ii <- which(s == g); sp <- strsplit(g, "\r", fixed = TRUE)[[1]]
    for (bb in 1:10) {
      n <- N[ii, bb]; k <- n > 0
      if (!any(k)) next
      yv <- C[ii, bb][k] / n[k]; xv <- P[ii, bb][k] / n[k]
      cal[[length(cal) + 1L]] <- data.frame(
        model = sp[1], annotation_type = sp[2], method = sp[3], band = BAND_LABELS[bb],
        x = mean(xv), y = mean(yv), se = if (sum(k) > 1L) sd(yv) / sqrt(sum(k)) else NA_real_,
        x_pooled = sum(P[ii, bb]) / sum(n), y_pooled = sum(C[ii, bb]) / sum(n),
        n_total = sum(n), units_used = sum(k), units_dropped = sum(!k),
        estimator = unit, stringsAsFactors = FALSE)
    }
    for (j in seq_along(THRESH)) {
      ns <- NS[ii, j]; k <- ns > 0
      if (!any(k)) next
      fv <- (ns[k] - TP[ii, j][k]) / ns[k]
      fdr[[length(fdr) + 1L]] <- data.frame(
        model = sp[1], annotation_type = sp[2], method = sp[3], t = THRESH[j],
        m = mean(fv), se = if (sum(k) > 1L) sd(fv) / sqrt(sum(k)) else NA_real_,
        units_used = sum(k), units_dropped = sum(!k), estimator = unit,
        stringsAsFactors = FALSE)
    }
  }
  cal <- do.call(rbind, cal); fdr <- do.call(rbind, fdr)
  cal$band <- factor(cal$band, levels = BAND_LABELS)
  cal$lo <- cal$y - 2 * cal$se; cal$hi <- cal$y + 2 * cal$se
  cal$reportable <- cal$n_total >= MIN_N & cal$units_used >= MIN_UNITS
  fdr$reportable <- fdr$units_used >= MIN_UNITS
  list(cal = cal, fdr = fdr)
}

CELL <- do.call(rbind, cellL); DRAW <- do.call(rbind, drawL)
FITS <- do.call(rbind, fitL); FITS <- rowsum(FITS, rownames(FITS))
est_cell <- unit_estimates(CELL, "cell")
est_draw <- unit_estimates(DRAW, "draw")

# per-fit estimator from the stratum sums
kp <- split_key(rownames(FITS))
fit_cal <- fit_fdr <- list()
for (i in seq_len(nrow(FITS))) {
  v <- FITS[i, ]
  kc <- v[1:10]; sy <- v[11:20]; sy2 <- v[21:30]; sx <- v[31:40]; nn <- v[41:50]
  kf <- v[51:58]; sf <- v[59:66]; sf2 <- v[67:74]
  for (bb in which(kc > 0)) {
    mu <- sy[bb] / kc[bb]
    se <- if (kc[bb] > 1) sqrt(max(sy2[bb] / kc[bb] - mu^2, 0) * kc[bb] / (kc[bb] - 1) / kc[bb]) else NA_real_
    fit_cal[[length(fit_cal) + 1L]] <- data.frame(
      model = kp[i, 1], annotation_type = kp[i, 2], method = kp[i, 3], band = BAND_LABELS[bb],
      x = sx[bb] / kc[bb], y = mu, se = se, n_total = nn[bb], units_used = kc[bb],
      estimator = "fit", stringsAsFactors = FALSE)
  }
  for (j in which(kf > 0)) {
    mu <- sf[j] / kf[j]
    se <- if (kf[j] > 1) sqrt(max(sf2[j] / kf[j] - mu^2, 0) * kf[j] / (kf[j] - 1) / kf[j]) else NA_real_
    fit_fdr[[length(fit_fdr) + 1L]] <- data.frame(
      model = kp[i, 1], annotation_type = kp[i, 2], method = kp[i, 3], t = THRESH[j],
      m = mu, se = se, units_used = kf[j], estimator = "fit", stringsAsFactors = FALSE)
  }
}
fit_cal <- do.call(rbind, fit_cal); fit_fdr <- do.call(rbind, fit_fdr)
fit_cal$band <- factor(fit_cal$band, levels = BAND_LABELS)
fit_cal$lo <- fit_cal$y - 2 * fit_cal$se; fit_cal$hi <- fit_cal$y + 2 * fit_cal$se
fit_cal$reportable <- fit_cal$n_total >= MIN_N & fit_cal$units_used >= MIN_UNITS
fit_fdr$reportable <- fit_fdr$units_used >= MIN_UNITS

## ---- validation: the rebuilt cell estimator must reproduce the frozen outputs
message("validating the rebuilt cell estimator ...")
rc <- readRDS(REF_CAL); rc <- rc[rc$method %in% M7, ]
mc <- merge(est_cell$cal, rc, by = c("model", "annotation_type", "method", "band"),
            suffixes = c("", ".ref"))
dy <- max(abs(mc$y - mc$y.ref)); dx <- max(abs(mc$x - mc$x.ref))
dn <- max(abs(mc$n_total - mc$n_total.ref)); du <- max(abs(mc$units_used - mc$cells_used))
message(sprintf("  calibration: %d/%d rows matched   max|dy| %.2e  max|dx| %.2e  max|dn| %g  max|dcells| %g",
                nrow(mc), nrow(rc), dy, dx, dn, du))

rf <- readRDS(REF_FDR); rf <- rf[rf$method %in% M7, ]
ck <- split_key(rownames(CELL))
cell_long <- do.call(rbind, lapply(seq_along(THRESH), function(j) data.frame(
  model = ck[, 1], annotation_type = ck[, 2], method = ck[, 3], job_dir = ck[, 4],
  S = as.numeric(ck[, 5]), phi = as.numeric(ck[, 6]), region_size = as.numeric(ck[, 7]),
  t = THRESH[j], nsel = CELL[, 30 + j], tp = CELL[, 38 + j], stringsAsFactors = FALSE)))
mf <- merge(cell_long, rf, by = c("model", "annotation_type", "method", "job_dir",
                                  "S", "phi", "region_size", "t"), suffixes = c("", ".ref"))
dsel <- max(abs(mf$nsel - mf$nsel.ref)); dtp <- max(abs(mf$tp - mf$tp.ref))
extra <- sum(!paste(cell_long$job_dir, cell_long$method, cell_long$S, cell_long$phi,
                    cell_long$region_size, cell_long$t)[cell_long$nsel > 0] %in%
               paste(rf$job_dir, rf$method, rf$S, rf$phi, rf$region_size, rf$t))
message(sprintf("  FDR counts : %d/%d rows matched   max|dnsel| %g  max|dtp| %g   selecting cells absent from ref: %d",
                nrow(mf), nrow(rf), dsel, dtp, extra))

ok_val <- nrow(mc) == nrow(rc) && dy < 1e-9 && dx < 1e-9 && dn == 0 && du == 0 &&
          nrow(mf) == nrow(rf) && dsel == 0 && dtp == 0 && extra == 0
if (!ok_val) stop("rebuilt cell estimator does NOT reproduce the frozen outputs; not writing", call. = FALSE)
message("  PASS - cell estimator reproduced exactly")

saveRDS(list(
  cal = rbind(est_cell$cal[, names(fit_cal)], est_draw$cal[, names(fit_cal)], fit_cal),
  fdr = rbind(est_cell$fdr[, names(fit_fdr)], est_draw$fdr[, names(fit_fdr)], fit_fdr),
  cal_draw_full = est_draw$cal, fdr_draw_full = est_draw$fdr,
  # unit-level counts (rows = units, rownames = keys; columns: N[10] C[10]
  # P[10] then NS[8] TP[8] SS[8] at THRESH), so any per-unit quantity - e.g.
  # the typical number of selections an iteration makes - is recoverable
  # without another pass over the tails
  units = list(cell = CELL, draw = DRAW),
  meta = list(methods = M7, edges = EDGES, thresholds = THRESH, min_n = MIN_N,
              min_units = MIN_UNITS, fits = n_fit, nonfinite_dropped = n_bad,
              unmatched = n_unm, validated = TRUE, created = Sys.time())),
  out_file)
message("wrote ", out_file)
