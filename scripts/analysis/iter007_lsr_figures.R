#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter007_lsr_figures.R
#
# >>> ITERATION 007 ONLY - TEMPORARY. See docs/autoresearch/iteration-007-REVERT.md
#
# Three Iteration 007 figures in the LSR's format - average precision, false
# discovery rate and PIP calibration - with the LD level down the rows and the
# annotation type across the columns. Every grid row is the sparse model, so each
# panel is one stratum; nothing is pooled across LD levels.
#
#   Rscript scripts/analysis/iter007_lsr_figures.R [iter007_dir] [fig_dir]
#
# Reads <dir>/analysis/iter007_results.rds from iter007_ld_results.R, which rebuilds
# the band and threshold counts from the floor-0 PIPs and stops unless they
# reproduce L1 exactly. A file not marked validated is refused.
#
# Format as iter006_lsr_figures.R: lsr_palette.R colours and shapes for the
# displayed methods that were run (sbayesrc was not), every other method grey;
# AP shows every method, FDR and calibration the displayed ones; no title or
# caption in the PDF.
#   * AP: each draw's mean AP over its valid regions, averaged over the 60 draws
#     of a stratum, with a 95% CI.
#   * FDR and calibration: the per-iteration estimator, +- 2 SE; reportable points
#     only (calibration: at least 50 variants and 20 draws; FDR: 20 draws).
# =============================================================================
suppressMessages({ library(ggplot2); library(grid) })
args <- commandArgs(trailingOnly = TRUE)
I7  <- if (length(args) >= 1L) args[1] else "results/iter007"
FIG <- if (length(args) >= 2L) args[2] else file.path(I7, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
source("scripts/analysis/lsr_palette.R")
R <- readRDS(file.path(I7, "analysis", "iter007_results.rds"))
stopifnot(isTRUE(R$meta$validated))

LD_LAB  <- c("in-sample" = "In-sample LD", "2000" = "Panel of 2,000", "1500" = "Panel of 1,500",
             "750" = "Panel of 750", "500" = "Panel of 500")
ANN_LAB <- c(binary = "Binary annotations", continuous = "Continuous annotations")
strata <- function(d) {
  d$ldf  <- factor(LD_LAB[d$ld], LD_LAB)
  d$annf <- factor(ANN_LAB[d$annotation_type], ANN_LAB)
  stopifnot(!anyNA(d$ldf), !anyNA(d$annf))
  d
}
W <- 7.2; H <- 9.6
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

present <- sort(unique(R$ap$method))
HL  <- intersect(M7, present)
OTH <- sprintf("other methods (%d)", length(setdiff(present, HL)))
hl_scales <- function() list(scale_colour_manual(values = PAL7[HL], limits = HL, name = NULL),
                             scale_shape_manual(values = SHP7[HL], limits = HL, name = NULL))

# ---- figure 1: average precision ----------------------------------------------
ap <- strata(R$ap)
stopifnot(all(ap$draws == 60L))
ord <- ap[ap$ld == "in-sample" & ap$annotation_type == "binary", ]
ord <- ord$method[order(ord$m)]
ap$method <- factor(ap$method, levels = ord)
ap$grp <- factor(ifelse(as.character(ap$method) %in% HL, as.character(ap$method), OTH), c(HL, OTH))
ap <- ap[order(ap$grp != OTH), ]                           # highlighted drawn over grey
p1 <- ggplot(ap, aes(m, method, colour = grp, shape = grp)) +
  geom_errorbarh(aes(xmin = m - 1.96 * se, xmax = m + 1.96 * se), height = .42, linewidth = .45) +
  geom_point(size = 1.9, stroke = .6) +
  facet_grid(ldf ~ annf) +
  scale_colour_manual(values = c(PAL7[HL], setNames(GREY, OTH)), limits = c(HL, OTH), name = NULL) +
  scale_shape_manual(values = c(SHP7[HL], setNames(16, OTH)), limits = c(HL, OTH), name = NULL) +
  labs(x = "Mean average precision over iterations (95% CI)", y = NULL) +
  th + theme(panel.grid.major.y = element_line(linewidth = .2, colour = "grey93"))
suppressWarnings(ggsave(file.path(FIG, "fig_ld_auprc.pdf"), p1, width = W, height = H))
message("figure 1 written")

# ---- figure 2: false discovery rate ----------------------------------------------
THRESH <- R$meta$thresholds
fd <- strata(R$fdr[R$fdr$reportable & R$fdr$method %in% HL, ])
fd$method <- factor(fd$method, HL)
fd$pos <- match(round(fd$t, 6), round(THRESH, 6))
stopifnot(!anyNA(fd$pos))
bound <- data.frame(pos = seq_along(THRESH), m = 1 - THRESH)
p2 <- ggplot(fd, aes(pos, m, group = method, colour = method, shape = method)) +
  geom_line(data = bound, aes(pos, m), inherit.aes = FALSE, colour = REF, linetype = "22", linewidth = .45) +
  geom_ribbon(aes(ymin = m - 2 * se, ymax = m + 2 * se, fill = method),
              colour = NA, alpha = .12, show.legend = FALSE, na.rm = TRUE) +
  geom_line(linewidth = .55) +
  geom_point(size = 1.4, stroke = .6) +
  facet_grid(ldf ~ annf) + hl_scales() +
  scale_fill_manual(values = PAL7[HL], guide = "none") +
  scale_x_continuous(breaks = seq_along(THRESH), labels = sub("^0", "", sprintf("%.2f", THRESH))) +
  labs(x = "PIP threshold t", y = "False discovery rate (mean over iterations, ±2 SE)") + th
suppressWarnings(ggsave(file.path(FIG, "fig_ld_fdr.pdf"), p2, width = W, height = H, device = cairo_pdf))
message("figure 2 written")

# ---- figure 3: PIP calibration ----------------------------------------------------
ca <- strata(R$cal[R$cal$reportable & R$cal$method %in% HL, ])
ca$method <- factor(ca$method, HL)
ca$lo <- ca$y - 2 * ca$se; ca$hi <- ca$y + 2 * ca$se
ca <- ca[order(ca$method, ca$x), ]
p3 <- ggplot(ca, aes(x, y, group = method, colour = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = REF, linewidth = .45) +
  geom_line(linewidth = .5) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .35, na.rm = TRUE) +
  geom_point(size = 1.5, stroke = .6) +
  facet_grid(ldf ~ annf) + hl_scales() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean PIP assigned within band", y = "Proportion of those variants causal") + th
suppressWarnings(ggsave(file.path(FIG, "fig_ld_calibration.pdf"), p3, width = W, height = H))
message("figure 3 written")
message(sprintf("highlighted: %s; others: %s", paste(HL, collapse = ", "),
                paste(setdiff(present, HL), collapse = ", ")))
