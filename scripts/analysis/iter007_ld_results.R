#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter007_ld_results.R
#
# >>> ITERATION 007 ONLY - TEMPORARY. See docs/autoresearch/iteration-007-REVERT.md
#
# First results for Iteration 007 (LD misspecification). Every row is the sparse
# model, so a stratum is LD level x annotation type - ten cells, one grid row
# each - and nothing is pooled across them.
#
#   Rscript scripts/analysis/iter007_ld_results.R [iter007_dir]
#
# Inputs (read only): L1/, piptail/ (floor 0), aux/, params_grid_iter007.csv.
# Output: <dir>/analysis/iter007_report.txt and iter007_results.rds.
#
#   0. checks: L1 joins the grid; every fit is present; band and threshold counts
#      rebuilt from the full PIPs reproduce the ones the collector froze in L1
#   1. failures: fits that errored, and fits that returned non-finite PIPs, per
#      method and stratum. Both are excluded from every metric below.
#   2. average precision. A draw (scenario: S x phi x iteration, ten regions) has
#      the mean AP of its valid regions; a stratum has the mean over its 60 draws,
#      +- 2 SE. Each panel level is compared with in-sample LD at the same
#      annotation type.
#   3. false discovery rate by the per-iteration estimator the LSR adopted, and
#      every crossing of 1 - t by more than 2 SE
#   4. calibration: the largest overconfidence in a reportable band at or above
#      PIP 0.5, and the two top bands
#   5. runtime: median seconds per scenario, and its ratio to in-sample LD
#
# Estimators follow iter006_lsr_figures.R: ten bands with edges 0, .01, .05, .1,
# .2, .5, .8, .9, .95, .99, 1; eight thresholds .05-.99; a draw with no variant in
# a band, or no selection at a threshold, is excluded rather than scored; a point
# is reportable with at least 50 variants (calibration) and 20 draws.
#
# CAVEATS, repeated in the report:
#   * UNPAIRED: each LD level is its own simulation (seed 1000 + row), so a
#     difference between levels carries draw-to-draw noise; its SE is
#     sqrt(se_a^2 + se_b^2).
#   * n = 1000, while the BEATRICE-family hyperparameters were tuned at n = 5000.
#   * fb_xregion ran the two-logit head (submitted before its removal).
# =============================================================================
options(width = 220, stringsAsFactors = FALSE)
`%||%` <- function(x, y) if (is.null(x)) y else x
args <- commandArgs(trailingOnly = TRUE)
I7   <- if (length(args) >= 1L) args[1] else "results/iter007"
OUT  <- file.path(I7, "analysis")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "iter007_report.txt"), split = TRUE)
hdr <- function(s) cat("\n", strrep("=", 110), "\n", s, "\n", strrep("=", 110), "\n", sep = "")

LD_LEV <- c("in-sample", "2000", "1500", "750", "500")    # least to most mismatched
ANN    <- c("binary", "continuous")
wide <- function(d, row, col, val, rows, cols = LD_LEV) {
  M <- matrix("", length(rows), length(cols), dimnames = list(rows, cols))
  for (i in seq_len(nrow(d))) if (d[[row]][i] %in% rows && d[[col]][i] %in% cols)
    M[d[[row]][i], d[[col]][i]] <- d[[val]][i]
  print(noquote(M))
}

# ---- 0. data and checks --------------------------------------------------------
g <- read.csv(file.path(I7, "params_grid_iter007.csv"))
g$job_dir <- sprintf("job_%03d_%s", g$job_id, g$label)
g$ld <- ifelse(is.na(g$n_ref), "in-sample", as.character(g$n_ref))
stopifnot(nrow(g) == 10L, all(g$model == "sparse"), setequal(g$ld, LD_LEV),
          setequal(g$annotation_type, ANN))

L1 <- do.call(rbind, lapply(sort(list.files(file.path(I7, "L1"), "^L1_.*[.]rds$", full.names = TRUE)), readRDS))
L1 <- merge(L1, g[, c("job_dir", "ld", "annotation_type")], by = "job_dir", all.x = TRUE)
if (anyNA(L1$ld)) stop("L1 rows that do not join the grid", call. = FALSE)
L1$failed <- as.logical(L1$failed); L1$bad_pip <- as.logical(L1$bad_pip)
L1$valid  <- !L1$failed & !L1$bad_pip
METHODS <- sort(unique(L1$method))
u <- unique(L1[, c("job_dir", "scenario_id", "S", "phi", "iter")])
stopifnot(length(METHODS) == 12L, nrow(L1) == 10L * 60L * 10L * 12L, nrow(u) == 600L,
          !anyDuplicated(u[, c("job_dir", "scenario_id")]), !anyDuplicated(u[, c("job_dir", "S", "phi", "iter")]))

hdr("0. Design and checks")
cat("Iteration 007: LD misspecification. Sparse model; GWAS n = 1000; regions of 100-400 variants;\n",
    "S {1,3,5} x phi {0.1,0.4} x 10 iterations = 60 draws per stratum; 10 annotations, 5 at fold 5.4.\n",
    "Strata: LD level (in-sample; reference panel of 2000, 1500, 750, 500) x annotation type (binary, continuous).\n\n",
    "CAVEATS\n",
    "  * UNPAIRED: each LD level is a separate simulation, so level-to-level differences carry draw-to-draw noise.\n",
    "  * n = 1000 here; the BEATRICE-family hyperparameters were tuned at n = 5000.\n",
    "  * fb_xregion ran the two-logit LassoNet head.\n\n", sep = "")
cat(sprintf("L1: %d fits = 10 rows x 60 draws x 10 regions x %d methods; all joined to the grid.\n", nrow(L1), length(METHODS)))

# ---- 1. failures -----------------------------------------------------------------
hdr("1. Failures per stratum: errored fits + fits with non-finite PIPs (of 600 fits per method per cell)")
fl <- aggregate(cbind(errored = failed, nonfinite = bad_pip) ~ method + annotation_type + ld, L1, sum)
fl$cell <- ifelse(fl$errored + fl$nonfinite == 0, "", sprintf("%d + %d", fl$errored, fl$nonfinite))
bad_methods <- sort(unique(fl$method[fl$errored + fl$nonfinite > 0]))
for (at in ANN) {
  cat(sprintf("\n-- %s annotations\n", at))
  if (length(bad_methods)) wide(fl[fl$annotation_type == at, ], "method", "ld", "cell", rows = bad_methods)
}
cat(sprintf("\nMethods with no failure anywhere: %s\n", paste(setdiff(METHODS, bad_methods), collapse = ", ")))

# ---- 2. average precision ---------------------------------------------------------
hdr("2. Average precision: mean over draws +- 2 SE; panel minus in-sample [SE of the difference], * if |z| > 2")
v  <- L1[L1$valid & is.finite(L1$ap), ]
dr <- aggregate(ap ~ job_dir + scenario_id + method + annotation_type + ld, v, mean)
apS <- do.call(rbind, lapply(split(dr, list(dr$annotation_type, dr$ld, dr$method), drop = TRUE), function(r)
  data.frame(annotation_type = r$annotation_type[1], ld = r$ld[1], method = r$method[1],
             m = mean(r$ap), se = sd(r$ap) / sqrt(nrow(r)), draws = nrow(r))))
ORD <- list(); apD <- list()
for (at in ANN) {
  a <- apS[apS$annotation_type == at, ]
  o <- a[a$ld == "in-sample", ]; ORD[[at]] <- o$method[order(-o$m)]
  a$cell <- sprintf("%.3f ± %.3f", a$m, 2 * a$se)
  cat(sprintf("\n-- %s annotations (draws per cell: %s); methods ordered by in-sample AP\n", at,
              paste(unique(range(a$draws)), collapse = "-")))
  wide(a, "method", "ld", "cell", rows = ORD[[at]])
  base <- setNames(a[a$ld == "in-sample", c("method", "m", "se")], c("method", "m0", "se0"))
  d <- merge(a[a$ld != "in-sample", ], base, by = "method")
  d$delta <- d$m - d$m0; d$se_d <- sqrt(d$se^2 + d$se0^2); d$z <- d$delta / d$se_d
  d$cell <- sprintf("%+.3f [%.3f]%s", d$delta, d$se_d, ifelse(abs(d$z) > 2, "*", " "))
  cat(sprintf("\n-- %s annotations: panel minus in-sample\n", at))
  wide(d, "method", "ld", "cell", rows = ORD[[at]], cols = LD_LEV[-1])
  a$rank <- ave(-a$m, a$ld, FUN = function(x) rank(x, ties.method = "min"))
  a$rcell <- as.character(a$rank)
  cat(sprintf("\n-- %s annotations: rank by AP within each LD level (1 = highest)\n", at))
  wide(a, "method", "ld", "rcell", rows = ORD[[at]])
  apD[[at]] <- d
}
apD <- do.call(rbind, apD)

# ---- per-draw band and threshold counts from the full PIPs ----------------------
EDGES  <- c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 1 + 1e-9)
NB     <- length(EDGES) - 1L
THRESH <- EDGES[3:NB]
BAND_LABELS <- sprintf("[%g,%g%s", EDGES[1:NB], pmin(EDGES[2:(NB + 1L)], 1), c(rep(")", NB - 1L), "]"))
MIN_N <- 50L; MIN_UNITS <- 20L
L1key <- paste(L1$job_dir, L1$scenario_id, L1$region_id, L1$method, sep = "\r")
stopifnot(!anyDuplicated(L1key))
REFCOLS <- c(paste0("n_band_", c("lo", "mid", "hi", "top")), paste0("c_band_", c("lo", "mid", "hi", "top")),
             paste0("nsel_at_", c(50, 80, 90, 95, 99)), paste0("tp_at_", c(50, 80, 90, 95, 99)))
revcum <- function(M) { for (k in (ncol(M) - 1L):1L) M[, k] <- M[, k] + M[, k + 1L]; M }
drawL <- list(); n_fit <- 0L; n_used <- 0L; max_dev <- 0
for (f in sort(list.files(file.path(I7, "piptail"), "^piptail_.*[.]rds$", full.names = TRUE))) {
  o <- readRDS(f)
  if (!identical(as.numeric(o$floor), 0)) stop(basename(f), " was not extracted at floor 0", call. = FALSE)
  tl <- o$tail; n_fit <- n_fit + length(tl)
  mi <- match(paste(o$job_dir, vapply(tl, function(z) as.integer(z$scenario_id), 0L),
                    vapply(tl, function(z) as.integer(z$region_id), 0L), vapply(tl, `[[`, "", "method"), sep = "\r"), L1key)
  if (anyNA(mi)) stop("fits in ", basename(f), " not found in L1", call. = FALSE)
  # extract_pip_tail.R keeps is.finite(pip) & pip >= floor, so a failed fit's NA
  # PIPs are counted as "below the floor". Failed fits are excluded here, so the
  # completeness check applies to the valid fits only.
  keep <- L1$valid[mi]; tl <- tl[keep]; mi <- mi[keep]; nf <- length(tl); n_used <- n_used + nf
  if (any(vapply(tl, function(z) z$n_below != 0, TRUE)))
    stop("valid fits with PIPs missing below the floor in ", basename(f), call. = FALSE)
  pl <- lapply(tl, `[[`, "pip"); yl <- lapply(tl, `[[`, "is_causal")
  fid <- rep.int(seq_len(nf), lengths(pl)); p <- unlist(pl, use.names = FALSE); y <- unlist(yl, use.names = FALSE)
  if (!all(is.finite(p))) stop("non-finite PIPs in a fit L1 calls valid: ", basename(f), call. = FALSE)
  bi <- findInterval(p, EDGES)
  if (!all(bi >= 1L & bi <= NB)) stop("PIPs outside [0, 1] in ", basename(f), call. = FALSE)
  rs <- rowsum(cbind(1, y, p), (fid - 1L) * NB + bi); pos <- as.integer(rownames(rs))
  Nv <- Cv <- Pv <- numeric(nf * NB); Nv[pos] <- rs[, 1]; Cv[pos] <- rs[, 2]; Pv[pos] <- rs[, 3]
  N <- matrix(Nv, nf, NB, byrow = TRUE); C <- matrix(Cv, nf, NB, byrow = TRUE); P <- matrix(Pv, nf, NB, byrow = TRUE)
  RN <- revcum(N); RC <- revcum(C)                       # {pip >= EDGES[k]} = bands k..NB
  mine <- cbind(rowSums(N[, 1:3]), rowSums(N[, 4:5]), rowSums(N[, 6:7]), rowSums(N[, 8:10]),
                rowSums(C[, 1:3]), rowSums(C[, 4:5]), rowSums(C[, 6:7]), rowSums(C[, 8:10]),
                RN[, 6:10], RC[, 6:10])
  max_dev <- max(max_dev, max(abs(mine - as.matrix(L1[mi, REFCOLS]))))
  r <- L1[mi, ]
  drawL[[f]] <- rowsum(cbind(N, C, P, RN[, 3:NB], RC[, 3:NB]),
                       paste(r$annotation_type, r$ld, r$method, r$job_dir, r$scenario_id, sep = "\r"))
}
cat(sprintf("\nPIPs: %d fits read, %d valid fits used; rebuilt band and threshold counts against L1: max |diff| %g\n",
            n_fit, n_used, max_dev))
if (n_fit != nrow(L1) || n_used != sum(L1$valid) || max_dev != 0)
  stop("rebuilt counts do not reproduce the pipeline's; stopping", call. = FALSE)

DRAW <- do.call(rbind, drawL)
kp   <- do.call(rbind, strsplit(rownames(DRAW), "\r", fixed = TRUE))
grp  <- paste(kp[, 1], kp[, 2], kp[, 3], sep = "\r")
cal <- fdr <- list()
for (gg in unique(grp)) {
  ii <- which(grp == gg); sp <- strsplit(gg, "\r", fixed = TRUE)[[1]]
  for (bb in seq_len(NB)) {
    n <- DRAW[ii, bb]; k <- n > 0
    if (!any(k)) next
    yv <- DRAW[ii, NB + bb][k] / n[k]; xv <- DRAW[ii, 2L * NB + bb][k] / n[k]
    cal[[length(cal) + 1L]] <- data.frame(annotation_type = sp[1], ld = sp[2], method = sp[3], band = BAND_LABELS[bb],
      x = mean(xv), y = mean(yv), se = if (sum(k) > 1L) sd(yv) / sqrt(sum(k)) else NA_real_,
      n_total = sum(n), units_used = sum(k), units_dropped = sum(!k))
  }
  for (j in seq_along(THRESH)) {
    ns <- DRAW[ii, 3L * NB + j]; k <- ns > 0
    if (!any(k)) next
    fv <- (ns[k] - DRAW[ii, 3L * NB + length(THRESH) + j][k]) / ns[k]
    fdr[[length(fdr) + 1L]] <- data.frame(annotation_type = sp[1], ld = sp[2], method = sp[3], t = THRESH[j],
      m = mean(fv), se = if (sum(k) > 1L) sd(fv) / sqrt(sum(k)) else NA_real_,
      units_used = sum(k), units_dropped = sum(!k))
  }
}
cal <- do.call(rbind, cal); fdr <- do.call(rbind, fdr)
cal$reportable <- cal$n_total >= MIN_N & cal$units_used >= MIN_UNITS
fdr$reportable <- fdr$units_used >= MIN_UNITS
fdr$z <- (fdr$m - (1 - fdr$t)) / fdr$se

# ---- 3. false discovery rate ------------------------------------------------------
hdr("3. False discovery rate, per iteration: mean [z against 1 - t]; blank = fewer than 20 draws selecting anything")
for (at in ANN) for (tt in c(0.5, 0.8, 0.9, 0.95)) {
  x <- fdr[fdr$annotation_type == at & abs(fdr$t - tt) < 1e-9 & fdr$reportable, ]
  x$cell <- sprintf("%.3f[%+.1f]", x$m, x$z)
  cat(sprintf("\n-- %s annotations, t = %g (bound %.2f)\n", at, tt, 1 - tt))
  wide(x, "method", "ld", "cell", rows = ORD[[at]])
}
cr <- fdr[fdr$reportable & is.finite(fdr$z) & fdr$z > 2, ]
cr <- cr[order(cr$annotation_type, cr$method, match(cr$ld, LD_LEV), cr$t), ]
cat("\n-- every crossing of 1 - t by more than 2 SE, all eight thresholds\n")
if (nrow(cr)) {
  cr$crossings <- sprintf("%s: t=%g %.3f(z %+.1f)", cr$ld, cr$t, cr$m, cr$z)
  for (k in unique(paste(cr$annotation_type, cr$method)))
    cat(sprintf("   %-26s %s\n", k, paste(cr$crossings[paste(cr$annotation_type, cr$method) == k], collapse = "; ")))
} else cat("   none\n")

# ---- 4. calibration ------------------------------------------------------------------
hdr("4. Calibration, per iteration: largest overconfidence (mean PIP - proportion causal) in a reportable band at PIP >= 0.5")
ca <- cal[cal$reportable & cal$x >= 0.5, ]
ca$gap <- ca$x - ca$y
worst <- do.call(rbind, lapply(split(ca, list(ca$annotation_type, ca$ld, ca$method), drop = TRUE),
                               function(r) r[which.max(r$gap), ]))
worst$cell <- sprintf("%+.2f %s", worst$gap, worst$band)
top <- cal[cal$reportable & cal$band %in% BAND_LABELS[9:10], ]
top <- do.call(rbind, lapply(split(top, list(top$annotation_type, top$ld, top$method), drop = TRUE), function(r)
  data.frame(annotation_type = r$annotation_type[1], ld = r$ld[1], method = r$method[1],
             cell = paste(sprintf("%.2f", r$y[match(BAND_LABELS[9:10], r$band)]), collapse = "/"))))
for (at in ANN) {
  cat(sprintf("\n-- %s annotations: largest gap and its band\n", at))
  wide(worst[worst$annotation_type == at, ], "method", "ld", "cell", rows = ORD[[at]])
  cat(sprintf("\n-- %s annotations: proportion causal in [0.95,0.99) / [0.99,1] (NA = band not reportable)\n", at))
  wide(top[top$annotation_type == at, ], "method", "ld", "cell", rows = ORD[[at]])
}

# ---- 5. runtime -----------------------------------------------------------------------
hdr("5. Runtime: median seconds per scenario (ten regions) [ratio to in-sample median]")
rt <- do.call(rbind, lapply(sort(list.files(file.path(I7, "aux"), "^aux_.*[.]rds$", full.names = TRUE)), function(f) {
  o <- readRDS(f); r <- Filter(function(x) identical(x$scope, "scenario_total"), o$runtime)
  data.frame(job_dir = o$job_dir, scenario_id = vapply(r, function(x) as.integer(x$scenario_id), 0L),
             method = vapply(r, `[[`, "", "method"), s = vapply(r, function(x) as.numeric(x$runtime_seconds), 0))
}))
rt <- merge(rt, g[, c("job_dir", "ld", "annotation_type")], by = "job_dir")
stopifnot(nrow(rt) == 600L * 12L)
rS <- aggregate(s ~ annotation_type + ld + method, rt, median)
b0 <- setNames(rS[rS$ld == "in-sample", c("annotation_type", "method", "s")], c("annotation_type", "method", "s0"))
rS <- merge(rS, b0, by = c("annotation_type", "method"))
rS$cell <- ifelse(rS$ld == "in-sample", sprintf("%.1f", rS$s), sprintf("%.1f [%.2f]", rS$s, rS$s / rS$s0))
for (at in ANN) {
  cat(sprintf("\n-- %s annotations\n", at))
  x <- rS[rS$annotation_type == at, ]
  wide(x, "method", "ld", "cell", rows = x$method[x$ld == "in-sample"][order(-x$s[x$ld == "in-sample"])])
}

saveRDS(list(failures = fl, ap_draws = dr, ap = apS, ap_delta = apD, fdr = fdr, cal = cal,
             cal_worst = worst, runtime = rt, runtime_median = rS, draw_counts = DRAW,
             meta = list(edges = EDGES, thresholds = THRESH, min_n = MIN_N, min_units = MIN_UNITS,
                         fits = n_fit, valid_fits = n_used, validated = TRUE, created = Sys.time())),
        file.path(OUT, "iter007_results.rds"))
cat(sprintf("\nwrote %s and iter007_results.rds\n", file.path(OUT, "iter007_report.txt")))
sink()
