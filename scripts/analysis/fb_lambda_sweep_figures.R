#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/fb_lambda_sweep_figures.R
#
# Figures and tables for the supplementary section on the lambda_l1 sweep of the
# cross-region model: AP, FDR and PIP calibration for five of the eight arms in the
# LSR's format, annotation pass-through for all eight, and every number the text
# quotes.
#
#   Rscript scripts/analysis/fb_lambda_sweep_figures.R [sweep_dir] [iter004_dir]
#
# Defaults: results/fb_lambda_sweep and results/iter004. Figures go to
# <sweep_dir>/figures and tables to <sweep_dir>/analysis; inputs are only read.
# RECOMPUTE=1 rebuilds the per-draw counts, which otherwise load from
# analysis/lambda_figures_data.rds.
#
# The five arms shown:
#   - the lowest lambda_l1 (0.03)
#   - the three with the highest average precision averaged over the four strata:
#     the identifiable head at 0.1223, the original two-logit head at 0.1223, and
#     the identifiable head at 0.25
#   - the highest (2.0)
#
# Format as iter004_lsr_figures.R: generative model down the rows, annotation type
# across the columns, never pooled; no captions baked in.
#   AP                average of the pipeline's AP over cells, ± 1.96 SE across cells,
#                     as in the LSR's AP figure. beatrice and polyfun_ldsc from the
#                     main study are shown in grey for reference.
#   FDR, calibration  the per-iteration estimator of iter004_calib_fdr_periter.R:
#                     a rate within each draw, averaged over draws, ± 2 SE, with ten
#                     bands and eight thresholds, shown with >= 50 variants and
#                     >= 20 draws. Built from the floor-0 PIPs, after the per-fit band
#                     and threshold counts are rebuilt and checked against the
#                     pipeline's frozen L1 counts.
#
# Annotation pass-through is measured per scenario from the shared head's
# importance vector (one per scenario; tracks 1-5 informative). An annotation is
# "let through" when its importance is at least 1e-8, the zeroing threshold used
# throughout the sweep.
#   precision    = informative let through / let through, over scenarios that let
#                  any through
#   FPR          = non-informative let through / 5
#   sensitivity  = informative let through / 5
# =============================================================================
suppressMessages({ library(ggplot2); library(grid) })
args  <- commandArgs(trailingOnly = TRUE)
SWEEP <- if (length(args) >= 1L) args[1] else "results/fb_lambda_sweep"
I004  <- if (length(args) >= 2L) args[2] else "results/iter004"
FIG <- file.path(SWEEP, "figures"); OUT <- file.path(SWEEP, "analysis")
for (d in c(FIG, OUT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
options(width = 220, stringsAsFactors = FALSE)
source("scripts/analysis/lsr_palette.R")            # GREY, REF

ARMS <- data.frame(
  method = c("fb_xregion_id_l03", "fb_xregion_id_l06", "fb_xregion_id", "fb_xregion_rerun",
             "fb_xregion_id_l25", "fb_xregion_id_l50", "fb_xregion_id_l100", "fb_xregion_id_l200"),
  lambda = c(0.03, 0.06, 0.1223, 0.1223, 0.25, 0.5, 1.0, 2.0),
  head   = c(rep("identifiable", 3), "original", rep("identifiable", 4)))
SHOW <- c("fb_xregion_id_l03", "fb_xregion_id", "fb_xregion_rerun", "fb_xregion_id_l25", "fb_xregion_id_l200")
LAB  <- c(fb_xregion_id_l03 = "λ = 0.03", fb_xregion_id = "λ = 0.1223",
          fb_xregion_rerun = "λ = 0.1223, original head", fb_xregion_id_l25 = "λ = 0.25",
          fb_xregion_id_l200 = "λ = 2")
COL  <- c(fb_xregion_id_l03 = "#56B4E9", fb_xregion_id = "#009E73", fb_xregion_rerun = "#D55E00",
          fb_xregion_id_l25 = "#0072B2", fb_xregion_id_l200 = "#CC79A7")
SHP  <- c(fb_xregion_id_l03 = 17, fb_xregion_id = 15, fb_xregion_rerun = 16,
          fb_xregion_id_l25 = 18, fb_xregion_id_l200 = 4)
REFS <- c("polyfun_ldsc", "beatrice")
MODEL_LAB <- c(sparse = "Sparse", sparse_inf = "Sparse + infinitesimal")
ARM_LAB   <- c(binary = "Binary annotations", continuous = "Continuous annotations")
strata <- function(d) {
  d$model <- factor(MODEL_LAB[as.character(d$model)], MODEL_LAB)
  d$arm   <- factor(ARM_LAB[as.character(d$annotation_type)], ARM_LAB)
  d
}
rename_sweep <- function(m) replace(m, m == "fb_xregion", "fb_xregion_rerun")
EDGES  <- c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 1 + 1e-9)
NB     <- length(EDGES) - 1L
THRESH <- EDGES[3:NB]
BAND_LABELS <- sprintf("[%g,%g%s", EDGES[1:NB], pmin(EDGES[2:(NB + 1L)], 1), c(rep(")", NB - 1L), "]"))
MIN_N <- 50L; MIN_UNITS <- 20L
TRUTH <- rep(c(TRUE, FALSE), each = 5L); ZERO_TOL <- 1e-8

th <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linewidth = .25, colour = "grey90"),
        panel.border = element_rect(colour = "grey70", linewidth = .35),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 7.4, face = "bold"),
        axis.title = element_text(size = 7.8),
        axis.text  = element_text(size = 6.9),
        legend.text  = element_text(size = 6.9),
        legend.key.size = unit(9, "pt"),
        legend.position = "right")
save_pdf <- function(p, file, w = 7.6, h = 5.3) ggsave(file.path(FIG, file), p, width = w, height = h, device = cairo_pdf)
arm_scales <- function() list(
  scale_colour_manual(values = COL[SHOW], limits = SHOW, labels = unname(LAB[SHOW]), name = NULL),
  scale_shape_manual(values = SHP[SHOW], limits = SHOW, labels = unname(LAB[SHOW]), name = NULL))

# ---- average precision -----------------------------------------------------------
L3s <- readRDS(file.path(SWEEP, "combined_scenario_metrics.rds")); L3s$method <- rename_sweep(L3s$method)
JOBS <- unique(L3s$job_dir[L3s$annotation_type != "none"])
L3s <- L3s[L3s$job_dir %in% JOBS & L3s$n_failed < L3s$n_fits & L3s$method %in% ARMS$method, ]
L3r <- readRDS(file.path(I004, "combined_scenario_metrics.rds"))
L3r <- L3r[L3r$job_dir %in% JOBS & L3r$method %in% REFS & L3r$n_failed < L3r$n_fits, ]
cols <- c("model", "annotation_type", "method", "ap")
apd <- rbind(L3s[, cols], L3r[, cols])
AP <- do.call(rbind, lapply(split(apd, list(apd$model, apd$annotation_type, apd$method), drop = TRUE), function(r) {
  v <- r$ap[is.finite(r$ap)]
  data.frame(model = r$model[1], annotation_type = r$annotation_type[1], method = r$method[1],
             m = mean(v), se = sd(v) / sqrt(length(v)), cells = length(v))
}))
a1 <- AP[AP$method %in% c(SHOW, REFS), ]
ylab <- c(LAB, polyfun_ldsc = "polyfun_ldsc", beatrice = "beatrice")
a1$row <- factor(ylab[a1$method], ylab[c(REFS, rev(SHOW))])
a1$grp <- factor(ifelse(a1$method %in% SHOW, a1$method, "reference"), c(SHOW, "reference"))
a1 <- strata(a1)
p1 <- ggplot(a1, aes(m, row, colour = grp, shape = grp)) +
  geom_errorbar(aes(xmin = m - 1.96 * se, xmax = m + 1.96 * se), orientation = "y", width = .42, linewidth = .45) +
  geom_point(size = 1.9, stroke = .6) +
  facet_grid(model ~ arm) +
  scale_colour_manual(values = c(COL[SHOW], reference = GREY), limits = c(SHOW, "reference"),
                      labels = c(unname(LAB[SHOW]), "main-study method"), name = NULL) +
  scale_shape_manual(values = c(SHP[SHOW], reference = 16), limits = c(SHOW, "reference"),
                     labels = c(unname(LAB[SHOW]), "main-study method"), name = NULL) +
  scale_x_continuous(breaks = seq(.5, .9, .1)) +
  labs(x = "Mean average precision over cells (95% CI)", y = NULL) +
  th + theme(panel.grid.major.y = element_line(linewidth = .2, colour = "grey93"))
save_pdf(p1, "fig_xregion_lambda_auprc.pdf")
message("AP figure written")

# ---- per-draw band and threshold counts, and annotation pass-through ------------
DATA <- file.path(OUT, "lambda_figures_data.rds")
if (file.exists(DATA) && Sys.getenv("RECOMPUTE") != "1") {
  dd <- readRDS(DATA); cal <- dd$cal; fdr <- dd$fdr; IMP <- dd$importance_units
  message("loaded ", DATA)
} else {
  L1 <- readRDS(file.path(SWEEP, "combined_fit_metrics.rds")); L1$method <- rename_sweep(L1$method)
  L1 <- L1[L1$job_dir %in% JOBS & L1$method %in% SHOW, ]
  u <- unique(L1[, c("job_dir", "scenario_id", "S", "phi", "iter")])
  if (any(duplicated(u[, c("job_dir", "scenario_id")]))) stop("scenario_id does not index (S, phi, iter)", call. = FALSE)
  L1key <- paste(L1$job_dir, L1$scenario_id, L1$region_id, L1$method, sep = "\r")
  REFCOLS <- c(paste0("n_band_", c("lo", "mid", "hi", "top")), paste0("c_band_", c("lo", "mid", "hi", "top")),
               paste0("nsel_at_", c(50, 80, 90, 95, 99)), paste0("tp_at_", c(50, 80, 90, 95, 99)))
  revcum <- function(M) { for (k in (ncol(M) - 1L):1L) M[, k] <- M[, k] + M[, k + 1L]; M }
  drawL <- list(); n_fit <- 0L; max_dev <- 0
  for (f in sort(list.files(file.path(SWEEP, "piptail"), "^piptail_.*[.]rds$", full.names = TRUE))) {
    o <- readRDS(f)
    if (!o$job_dir %in% JOBS) next
    if (!identical(as.numeric(o$floor), 0)) stop(basename(f), " was not extracted at floor 0", call. = FALSE)
    tl <- Filter(function(z) rename_sweep(z$method) %in% SHOW && length(z$pip) > 0L, o$tail); rm(o)
    nf <- length(tl)
    if (any(vapply(tl, function(z) z$n_below != 0, TRUE))) stop("PIPs missing below the floor in ", basename(f), call. = FALSE)
    meth <- rename_sweep(vapply(tl, `[[`, "", "method"))
    mi <- match(paste(tl[[1]]$job_dir, vapply(tl, function(z) as.integer(z$scenario_id), 0L),
                      vapply(tl, function(z) as.integer(z$region_id), 0L), meth, sep = "\r"), L1key)
    if (anyNA(mi)) stop("fits in ", basename(f), " not found in L1", call. = FALSE)
    pl <- lapply(tl, `[[`, "pip"); yl <- lapply(tl, `[[`, "is_causal"); rm(tl)
    fid <- rep.int(seq_len(nf), lengths(pl)); p <- unlist(pl, use.names = FALSE); y <- unlist(yl, use.names = FALSE); rm(pl, yl)
    ok <- is.finite(p); p <- p[ok]; y <- y[ok]; fid <- fid[ok]
    ix <- (fid - 1L) * NB + findInterval(p, EDGES)
    rs <- rowsum(cbind(1, y, p), ix); pos <- as.integer(rownames(rs))
    Nv <- Cv <- Pv <- numeric(nf * NB); Nv[pos] <- rs[, 1]; Cv[pos] <- rs[, 2]; Pv[pos] <- rs[, 3]
    N <- matrix(Nv, nf, NB, byrow = TRUE); C <- matrix(Cv, nf, NB, byrow = TRUE); P <- matrix(Pv, nf, NB, byrow = TRUE)
    RN <- revcum(N); RC <- revcum(C)
    mine <- cbind(rowSums(N[, 1:3]), rowSums(N[, 4:5]), rowSums(N[, 6:7]), rowSums(N[, 8:10]),
                  rowSums(C[, 1:3]), rowSums(C[, 4:5]), rowSums(C[, 6:7]), rowSums(C[, 8:10]), RN[, 6:10], RC[, 6:10])
    max_dev <- max(max_dev, max(abs(mine - as.matrix(L1[mi, REFCOLS]))))
    r <- L1[mi, ]
    drawL[[f]] <- rowsum(cbind(N, C, P, RN[, 3:NB], RC[, 3:NB]),
                         paste(r$model, r$annotation_type, r$method, r$job_dir, r$scenario_id, sep = "\r"))
    n_fit <- n_fit + nf
    message(sprintf("  %s: %d fits", basename(f), nf))
  }
  message(sprintf("fits read: %d of %d; rebuilt counts against L1: max |diff| %g", n_fit, nrow(L1), max_dev))
  if (n_fit != nrow(L1) || max_dev != 0) stop("rebuilt counts do not reproduce the pipeline's; not plotting", call. = FALSE)

  DRAW <- do.call(rbind, drawL)
  kp <- do.call(rbind, strsplit(rownames(DRAW), "\r", fixed = TRUE))
  grp <- paste(kp[, 1], kp[, 2], kp[, 3], sep = "\r")
  cal <- fdr <- list()
  for (g in unique(grp)) {
    ii <- which(grp == g); sp <- strsplit(g, "\r", fixed = TRUE)[[1]]
    for (bb in seq_len(NB)) {
      n <- DRAW[ii, bb]; k <- n > 0
      if (!any(k)) next
      yv <- DRAW[ii, NB + bb][k] / n[k]; xv <- DRAW[ii, 2L * NB + bb][k] / n[k]
      cal[[length(cal) + 1L]] <- data.frame(model = sp[1], annotation_type = sp[2], method = sp[3], band = BAND_LABELS[bb],
        x = mean(xv), y = mean(yv), se = if (sum(k) > 1L) sd(yv) / sqrt(sum(k)) else NA_real_,
        n_total = sum(n), units_used = sum(k), units_dropped = sum(!k))
    }
    for (j in seq_along(THRESH)) {
      ns <- DRAW[ii, 3L * NB + j]; k <- ns > 0
      if (!any(k)) next
      fv <- (ns[k] - DRAW[ii, 3L * NB + length(THRESH) + j][k]) / ns[k]
      fdr[[length(fdr) + 1L]] <- data.frame(model = sp[1], annotation_type = sp[2], method = sp[3], t = THRESH[j],
        m = mean(fv), se = if (sum(k) > 1L) sd(fv) / sqrt(sum(k)) else NA_real_, units_used = sum(k), units_dropped = sum(!k))
    }
  }
  cal <- do.call(rbind, cal); fdr <- do.call(rbind, fdr)
  cal$lo <- cal$y - 2 * cal$se; cal$hi <- cal$y + 2 * cal$se
  cal$reportable <- cal$n_total >= MIN_N & cal$units_used >= MIN_UNITS
  fdr$reportable <- fdr$units_used >= MIN_UNITS

  # annotation pass-through, all eight arms
  design <- unique(L3s[, c("job_dir", "model", "annotation_type")])
  impL <- lapply(sort(list.files(file.path(SWEEP, "aux"), "^aux_.*[.]rds$", full.names = TRUE)), function(f) {
    o <- readRDS(f)
    if (!o$job_dir %in% JOBS) return(NULL)
    im <- Filter(function(z) rename_sweep(z$method) %in% ARMS$method && !isTRUE(z$joint_fallback), o$importance)
    if (!length(im)) return(NULL)
    meta <- data.frame(job_dir = o$job_dir, scenario_id = vapply(im, function(z) as.integer(z$scenario_id), 0L),
                       method = rename_sweep(vapply(im, `[[`, "", "method")))
    keep <- !duplicated(meta[, c("scenario_id", "method")])           # the shared head: one vector per scenario
    V <- t(vapply(im[keep], function(z) { v <- rep(NA_real_, 10L); ix <- as.integer(z$importance$annotation_index)
      if (min(ix) == 0L) ix <- ix + 1L; v[ix] <- as.numeric(z$importance$importance); v }, numeric(10)))
    through <- V >= ZERO_TOL
    cbind(meta[keep, ], let_through = rowSums(through), informative = rowSums(through[, TRUTH, drop = FALSE]),
          inert = rowSums(through[, !TRUTH, drop = FALSE]))
  })
  IMP <- merge(do.call(rbind, impL), design, by = "job_dir")
  saveRDS(list(cal = cal, fdr = fdr, importance_units = IMP,
               meta = list(arms = ARMS, shown = SHOW, edges = EDGES, thresholds = THRESH, min_n = MIN_N,
                           min_units = MIN_UNITS, fits = n_fit, validated = TRUE, created = Sys.time())), DATA)
}

# ---- FDR and calibration figures ---------------------------------------------------
fd <- strata(fdr[fdr$reportable & fdr$method %in% SHOW, ])
fd$method <- factor(fd$method, SHOW); fd$pos <- match(fd$t, THRESH)
bound <- data.frame(pos = seq_along(THRESH), m = 1 - THRESH)
p2 <- ggplot(fd, aes(pos, m, group = method, colour = method, shape = method)) +
  geom_line(data = bound, aes(pos, m), inherit.aes = FALSE, colour = REF, linetype = "22", linewidth = .45) +
  geom_ribbon(aes(ymin = m - 2 * se, ymax = m + 2 * se, fill = method), colour = NA, alpha = .12, show.legend = FALSE, na.rm = TRUE) +
  geom_line(linewidth = .55) + geom_point(size = 1.4, stroke = .6) +
  facet_grid(model ~ arm) + arm_scales() +
  scale_fill_manual(values = COL[SHOW], guide = "none") +
  scale_x_continuous(breaks = seq_along(THRESH), labels = sub("^0", "", sprintf("%.2f", THRESH))) +
  labs(x = "PIP threshold t", y = "False discovery rate (mean over iterations, ±2 SE)") + th
save_pdf(p2, "fig_xregion_lambda_fdr.pdf")
ca <- strata(cal[cal$reportable & cal$method %in% SHOW, ])
ca$method <- factor(ca$method, SHOW); ca <- ca[order(ca$method, ca$x), ]
p3 <- ggplot(ca, aes(x, y, group = method, colour = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = REF, linewidth = .45) +
  geom_line(linewidth = .5) + geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .35, na.rm = TRUE) +
  geom_point(size = 1.5, stroke = .6) +
  facet_grid(model ~ arm) + arm_scales() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean PIP assigned within band", y = "Proportion of those variants causal") + th
save_pdf(p3, "fig_xregion_lambda_calibration.pdf")
message("FDR and calibration figures written")

# ---- numbers for the text -----------------------------------------------------------
IMPS <- do.call(rbind, lapply(split(IMP, list(IMP$method, IMP$model, IMP$annotation_type), drop = TRUE), function(d) {
  k <- d$let_through > 0; prec <- d$informative[k] / d$let_through[k]
  data.frame(method = d$method[1], model = d$model[1], annotation_type = d$annotation_type[1], scenarios = nrow(d),
             let_through = mean(d$let_through), informative = mean(d$informative), inert = mean(d$inert),
             precision = if (any(k)) mean(prec) else NA_real_, precision_se = if (sum(k) > 1) sd(prec) / sqrt(sum(k)) else NA_real_,
             fpr = mean(d$inert / 5), fpr_se = sd(d$inert / 5) / sqrt(nrow(d)),
             sensitivity = mean(d$informative / 5), none_through = mean(!k))
}))
IMPS <- IMPS[order(IMPS$model, IMPS$annotation_type, match(IMPS$method, ARMS$method)), ]
sw <- readRDS(file.path(OUT, "sweep_analysis.rds"))
sink(file.path(OUT, "lambda_supplement_numbers.txt"), split = TRUE)
cat("== Mean AP over the four strata, all arms (choice of the five shown) ==\n")
a8 <- AP[AP$method %in% ARMS$method, ]
print(sort(tapply(a8$m, a8$method, mean), decreasing = TRUE), digits = 4)
cat("\n== AP per stratum, mean ± 2 SE over cells ==\n")
print(within(AP[order(AP$model, AP$annotation_type, match(AP$method, c(ARMS$method, REFS))), ],
             { m <- round(m, 3); se2 <- round(2 * se, 3); se <- NULL }), row.names = FALSE)
cat("\n== AP paired by cell against Iteration 004's fb_xregion, beatrice, polyfun_ldsc (sweep_analysis.rds) ==\n")
pp <- sw$ap_paired; pp$diff <- round(pp$diff, 3); pp$diff_se <- round(pp$diff_se, 3)
print(pp[order(pp$versus, pp$model, pp$annotation_type, match(pp$method, ARMS$method)), ], row.names = FALSE)
cat("\n== FDR (per iteration) at t = .5 .8 .9 .95 .99 for the five arms: mean [z against 1 - t] ==\n")
fdr$z <- (fdr$m - (1 - fdr$t)) / fdr$se
tab <- do.call(rbind, lapply(split(fdr[fdr$reportable & fdr$t %in% c(.5, .8, .9, .95, .99), ],
  list(fdr$model[fdr$reportable & fdr$t %in% c(.5, .8, .9, .95, .99)], fdr$annotation_type[fdr$reportable & fdr$t %in% c(.5, .8, .9, .95, .99)],
       fdr$method[fdr$reportable & fdr$t %in% c(.5, .8, .9, .95, .99)]), drop = TRUE), function(d) {
  d <- d[order(d$t), ]
  data.frame(model = d$model[1], annotation_type = d$annotation_type[1], method = d$method[1],
             fdr = paste(sprintf("%.3f[%+.1f]", d$m, d$z), collapse = " "))
}))
print(tab[order(tab$model, tab$annotation_type, match(tab$method, SHOW)), ], row.names = FALSE, right = FALSE)
cat("\ncrossings of 1 - t by more than 2 SE:\n")
cr <- fdr[fdr$reportable & fdr$z > 2, c("model", "annotation_type", "method", "t", "m", "se", "z")]
if (nrow(cr)) print(cr, row.names = FALSE, digits = 3) else cat("  none\n")
cat("\n== Calibration: worst overconfidence (x - y) in bands with mean PIP >= 0.5, and the top two bands ==\n")
cal$gap <- cal$x - cal$y
hi <- cal[cal$reportable & cal$x >= 0.5, ]
w <- do.call(rbind, lapply(split(hi, list(hi$model, hi$annotation_type, hi$method), drop = TRUE), function(d) {
  i <- which.max(d$gap); top <- cal[cal$reportable & cal$model == d$model[1] & cal$annotation_type == d$annotation_type[1] &
                                    cal$method == d$method[1] & cal$band %in% c("[0.95,0.99)", "[0.99,1]"), ]
  data.frame(model = d$model[1], annotation_type = d$annotation_type[1], method = d$method[1], worst_band = d$band[i],
             x = round(d$x[i], 3), y = round(d$y[i], 3), gap = round(d$gap[i], 3), z = round(d$gap[i] / d$se[i], 1),
             top_bands_y = paste(sprintf("%.3f", top$y[order(top$x)]), collapse = "/"))
}))
print(w[order(w$model, w$annotation_type, match(w$method, SHOW)), ], row.names = FALSE)
cat("\n== Annotation pass-through, all eight arms (per scenario) ==\n")
print(within(IMPS, { let_through <- round(let_through, 2); informative <- round(informative, 2); inert <- round(inert, 2)
  precision <- round(precision, 3); precision_se <- round(precision_se, 3); fpr <- round(fpr, 3); fpr_se <- round(fpr_se, 3)
  sensitivity <- round(sensitivity, 3); none_through <- round(none_through, 3) }), row.names = FALSE)
cat("\n== Runtime: identifiable head at lambda 0.1223, and Iteration 004's fb_xregion ==\n")
rt <- sw$runtime
print(rt[rt$method %in% c("fb_xregion_id", "fb_xregion"), c("model", "annotation_type", "method", "scenarios", "median_s", "q25_s",
      "q75_s", "mean_s", "cpu_hours", "ratio_median", "ratio_q25", "ratio_q75", "ratio_of_totals")], row.names = FALSE, digits = 3)
sink()
saveRDS(list(ap = AP, importance = IMPS, calibration_worst = w, fdr = fdr, cal = cal), file.path(OUT, "lambda_supplement_numbers.rds"))
message("numbers written to ", file.path(OUT, "lambda_supplement_numbers.txt"))
