#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter006_lsr_figures.R
#
# >>> ITERATION 006 ONLY - TEMPORARY. See docs/autoresearch/iteration-006-REVERT.md
#
# Three Iteration 006 figures in the LSR's format - average precision, false
# discovery rate and PIP calibration - with the eight annotation-to-causality
# relationships down the rows. Annotations are continuous throughout, so each row
# is one stratum; nothing is pooled across relationships.
#
#   Rscript scripts/analysis/iter006_lsr_figures.R [iter006_dir] [fig_dir]
#
# Format as iter004_lsr_figures.R draws the LSR figures (branch lsr/feedback1-figures):
#   * colour and shape from lsr_palette.R. The displayed methods present here are
#     coloured (finemap_inf was not run); every other method is grey.
#   * average precision shows every method with the displayed ones highlighted;
#     calibration and FDR show the displayed ones only
#   * no title or caption in the PDF
#
# Estimators:
#   * AP. The pipeline's per-fit AP (.compute_ap_exact), averaged over the ten
#     regions of a draw, then over the 100 draws (S x phi x iteration) of a
#     relationship. 95% CI = ± 1.96 SE across draws. The LSR's AP figure takes its SE
#     across cells, but each relationship here has only four S x phi cells.
#   * FDR and calibration. The per-iteration ("draw") estimator of
#     iter004_calib_fdr_periter.R, which the LSR adopted: a rate within each draw,
#     averaged over draws, ± 2 SE.
#       - ten bands with edges 0, .01, .05, .1, .2, .5, .8, .9, .95, .99, 1
#       - eight thresholds, .05 to .99
#       - a draw with no variant in a band, or no selection at a threshold, is
#         excluded rather than scored
#       - a point is shown only with at least 50 variants and 20 draws (FDR: 20 draws)
#   * Inputs are the floor-0 PIPs in results/iter006/piptail, so the lowest band is
#     read directly rather than from sub-floor counts. Before anything is plotted,
#     the per-fit counts are rebuilt and checked against the band and threshold
#     counts the pipeline froze in L1; the script stops if any differ.
# =============================================================================
suppressMessages({ library(ggplot2); library(grid) })
args <- commandArgs(trailingOnly = TRUE)
I6  <- if (length(args) >= 1L) args[1] else "results/iter006"
FIG <- if (length(args) >= 2L) args[2] else file.path(I6, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(I6, "analysis"), recursive = TRUE, showWarnings = FALSE)
source("scripts/analysis/lsr_palette.R")

REL <- c("null", "additive", "wavg2", "valthresh", "valthresh_wavg2", "square", "cubic", "cosine10")
REL_LAB <- c(null            = "Null: annotations carry no information",
             additive        = "Additive (the log-linear form the methods assume)",
             wavg2           = "Weighted average over ±2 variants",
             valthresh       = "Threshold: a value counts only above 1",
             valthresh_wavg2 = "Threshold, then weighted average over ±2 variants",
             square          = "Square",
             cubic           = "Cubic",
             cosine10        = "Damped-cosine weighted average over ±10 variants")
EDGES  <- c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 1 + 1e-9)
NB     <- length(EDGES) - 1L                         # ten bands
THRESH <- EDGES[3:NB]                                # .05 ... .99, eight thresholds
BAND_LABELS <- sprintf("[%g,%g%s", EDGES[1:NB], pmin(EDGES[2:(NB + 1L)], 1), c(rep(")", NB - 1L), "]"))
MIN_N <- 50L; MIN_UNITS <- 20L
W <- 7.2; H <- 11

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
relf <- function(x) factor(REL_LAB[x], REL_LAB[REL])

# ---- figure 1: average precision ------------------------------------------------
L2 <- readRDS(file.path(I6, "combined_replicate_metrics.rds"))
stopifnot(identical(sort(unique(L2$annotation_type)), "continuous"), setequal(unique(L2$relationship), REL))
present <- sort(unique(L2$method))
HL  <- intersect(M7, present)
OTH <- sprintf("other methods (%d)", length(setdiff(present, HL)))
ap <- do.call(rbind, lapply(split(L2, list(L2$relationship, L2$method), drop = TRUE), function(r) {
  v <- r$ap[is.finite(r$ap)]
  data.frame(relationship = r$relationship[1], method = r$method[1],
             m = mean(v), se = sd(v) / sqrt(length(v)), n = length(v))
}))
stopifnot(all(ap$n == 100L))
ord <- ap[ap$relationship == "additive", ]; ord <- ord$method[order(ord$m)]
ap$method <- factor(ap$method, levels = ord)
ap$grp <- factor(ifelse(as.character(ap$method) %in% HL, as.character(ap$method), OTH), c(HL, OTH))
ap$rel <- relf(ap$relationship)
ap <- ap[order(ap$grp != OTH), ]                          # highlighted drawn over grey
p1 <- ggplot(ap, aes(m, method, colour = grp, shape = grp)) +
  geom_errorbarh(aes(xmin = m - 1.96 * se, xmax = m + 1.96 * se), height = .42, linewidth = .45) +
  geom_point(size = 1.9, stroke = .6) +
  facet_wrap(~ rel, ncol = 1) +
  scale_colour_manual(values = c(PAL7[HL], setNames(GREY, OTH)), limits = c(HL, OTH), name = NULL) +
  scale_shape_manual(values = c(SHP7[HL], setNames(16, OTH)), limits = c(HL, OTH), name = NULL) +
  scale_x_continuous(breaks = seq(.3, .9, .2)) +
  labs(x = "Mean average precision over iterations (95% CI)", y = NULL) +
  th + theme(panel.grid.major.y = element_line(linewidth = .2, colour = "grey93"))
suppressWarnings(ggsave(file.path(FIG, "fig_iter006_auprc.pdf"), p1, width = W, height = H))
message("figure 1 written")

# ---- per-draw band and threshold counts from the full PIPs --------------------
L1 <- readRDS(file.path(I6, "combined_fit_metrics.rds"))
u <- unique(L1[, c("job_dir", "scenario_id", "S", "phi", "iter")])
if (any(duplicated(u[, c("job_dir", "scenario_id")])) || any(duplicated(u[, c("job_dir", "S", "phi", "iter")])))
  stop("scenario_id does not index (S, phi, iter) one-to-one", call. = FALSE)
L1key <- paste(L1$job_dir, L1$scenario_id, L1$region_id, L1$method, sep = "\r")
REFCOLS <- c(paste0("n_band_", c("lo", "mid", "hi", "top")), paste0("c_band_", c("lo", "mid", "hi", "top")),
             paste0("nsel_at_", c(50, 80, 90, 95, 99)), paste0("tp_at_", c(50, 80, 90, 95, 99)))
revcum <- function(M) { for (k in (ncol(M) - 1L):1L) M[, k] <- M[, k] + M[, k + 1L]; M }
drawL <- list(); n_fit <- 0L; max_dev <- 0
for (f in sort(list.files(file.path(I6, "piptail"), "^piptail_.*[.]rds$", full.names = TRUE))) {
  o <- readRDS(f)
  if (!identical(as.numeric(o$floor), 0)) stop(basename(f), " was not extracted at floor 0", call. = FALSE)
  tl <- Filter(function(z) length(z$pip) > 0L, o$tail); nf <- length(tl)
  if (any(vapply(tl, function(z) z$n_below != 0, TRUE))) stop("PIPs missing below the floor in ", basename(f), call. = FALSE)
  mi <- match(paste(o$job_dir, vapply(tl, function(z) as.integer(z$scenario_id), 0L),
                    vapply(tl, function(z) as.integer(z$region_id), 0L), vapply(tl, `[[`, "", "method"), sep = "\r"), L1key)
  if (anyNA(mi)) stop("fits in ", basename(f), " not found in L1", call. = FALSE)
  pl <- lapply(tl, `[[`, "pip"); yl <- lapply(tl, `[[`, "is_causal")
  fid <- rep.int(seq_len(nf), lengths(pl)); p <- unlist(pl, use.names = FALSE); y <- unlist(yl, use.names = FALSE)
  ok <- is.finite(p); p <- p[ok]; y <- y[ok]; fid <- fid[ok]
  ix <- (fid - 1L) * NB + findInterval(p, EDGES)
  rs <- rowsum(cbind(1, y, p), ix); pos <- as.integer(rownames(rs))
  Nv <- Cv <- Pv <- numeric(nf * NB); Nv[pos] <- rs[, 1]; Cv[pos] <- rs[, 2]; Pv[pos] <- rs[, 3]
  N <- matrix(Nv, nf, NB, byrow = TRUE); C <- matrix(Cv, nf, NB, byrow = TRUE); P <- matrix(Pv, nf, NB, byrow = TRUE)
  RN <- revcum(N); RC <- revcum(C)                    # {pip >= EDGES[k]} = bands k..NB
  mine <- cbind(rowSums(N[, 1:3]), rowSums(N[, 4:5]), rowSums(N[, 6:7]), rowSums(N[, 8:10]),
                rowSums(C[, 1:3]), rowSums(C[, 4:5]), rowSums(C[, 6:7]), rowSums(C[, 8:10]),
                RN[, 6:10], RC[, 6:10])
  max_dev <- max(max_dev, max(abs(mine - as.matrix(L1[mi, REFCOLS]))))
  r <- L1[mi, ]
  drawL[[f]] <- rowsum(cbind(N, C, P, RN[, 3:NB], RC[, 3:NB]),
                       paste(r$relationship, r$method, r$job_dir, r$scenario_id, sep = "\r"))
  n_fit <- n_fit + nf
}
message(sprintf("fits read: %d; rebuilt band and threshold counts against L1: max |diff| %g", n_fit, max_dev))
if (n_fit != nrow(L1) || max_dev != 0) stop("rebuilt counts do not reproduce the pipeline's; not plotting", call. = FALSE)

DRAW <- do.call(rbind, drawL)
kp <- do.call(rbind, strsplit(rownames(DRAW), "\r", fixed = TRUE))
grp <- paste(kp[, 1], kp[, 2], sep = "\r")
cal <- fdr <- list()
for (g in unique(grp)) {
  ii <- which(grp == g); sp <- strsplit(g, "\r", fixed = TRUE)[[1]]
  for (bb in seq_len(NB)) {
    n <- DRAW[ii, bb]; k <- n > 0
    if (!any(k)) next
    yv <- DRAW[ii, NB + bb][k] / n[k]; xv <- DRAW[ii, 2L * NB + bb][k] / n[k]
    cal[[length(cal) + 1L]] <- data.frame(relationship = sp[1], method = sp[2], band = BAND_LABELS[bb],
      x = mean(xv), y = mean(yv), se = if (sum(k) > 1L) sd(yv) / sqrt(sum(k)) else NA_real_,
      n_total = sum(n), units_used = sum(k), units_dropped = sum(!k))
  }
  for (j in seq_along(THRESH)) {
    ns <- DRAW[ii, 3L * NB + j]; k <- ns > 0
    if (!any(k)) next
    fv <- (ns[k] - DRAW[ii, 3L * NB + length(THRESH) + j][k]) / ns[k]
    fdr[[length(fdr) + 1L]] <- data.frame(relationship = sp[1], method = sp[2], t = THRESH[j],
      m = mean(fv), se = if (sum(k) > 1L) sd(fv) / sqrt(sum(k)) else NA_real_,
      units_used = sum(k), units_dropped = sum(!k))
  }
}
cal <- do.call(rbind, cal); fdr <- do.call(rbind, fdr)
cal$lo <- cal$y - 2 * cal$se; cal$hi <- cal$y + 2 * cal$se
cal$reportable <- cal$n_total >= MIN_N & cal$units_used >= MIN_UNITS
fdr$reportable <- fdr$units_used >= MIN_UNITS
saveRDS(list(ap = ap, cal = cal, fdr = fdr, units = DRAW,
             meta = list(edges = EDGES, thresholds = THRESH, min_n = MIN_N, min_units = MIN_UNITS,
                         highlighted = HL, fits = n_fit, validated = TRUE, created = Sys.time())),
        file.path(I6, "analysis", "iter006_figures_data.rds"))

hl_scales <- function() list(scale_colour_manual(values = PAL7[HL], limits = HL, name = NULL),
                             scale_shape_manual(values = SHP7[HL], limits = HL, name = NULL))

# ---- figure 2: false discovery rate ----------------------------------------
fd <- fdr[fdr$reportable & fdr$method %in% HL, ]
fd$method <- factor(fd$method, HL); fd$rel <- relf(fd$relationship)
fd$pos <- match(fd$t, THRESH)
bound <- data.frame(pos = seq_along(THRESH), m = 1 - THRESH)
p2 <- ggplot(fd, aes(pos, m, group = method, colour = method, shape = method)) +
  geom_line(data = bound, aes(pos, m), inherit.aes = FALSE, colour = REF, linetype = "22", linewidth = .45) +
  geom_ribbon(aes(ymin = m - 2 * se, ymax = m + 2 * se, fill = method),
              colour = NA, alpha = .12, show.legend = FALSE, na.rm = TRUE) +
  geom_line(linewidth = .55) +
  geom_point(size = 1.4, stroke = .6) +
  facet_wrap(~ rel, ncol = 1) + hl_scales() +
  scale_fill_manual(values = PAL7[HL], guide = "none") +
  scale_x_continuous(breaks = seq_along(THRESH), labels = sub("^0", "", sprintf("%.2f", THRESH))) +
  labs(x = "PIP threshold t", y = "False discovery rate (mean over iterations, ±2 SE)") + th
ggsave(file.path(FIG, "fig_iter006_fdr.pdf"), p2, width = W, height = H)
message("figure 2 written")

# ---- figure 3: PIP calibration -----------------------------------------------
ca <- cal[cal$reportable & cal$method %in% HL, ]
ca$method <- factor(ca$method, HL); ca$rel <- relf(ca$relationship)
ca <- ca[order(ca$method, ca$x), ]
p3 <- ggplot(ca, aes(x, y, group = method, colour = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = REF, linewidth = .45) +
  geom_line(linewidth = .5) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .35, na.rm = TRUE) +
  geom_point(size = 1.5, stroke = .6) +
  facet_wrap(~ rel, ncol = 1) + hl_scales() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean PIP assigned within band", y = "Proportion of those variants causal") + th
ggsave(file.path(FIG, "fig_iter006_calibration.pdf"), p3, width = W, height = H)
message("figure 3 written")
