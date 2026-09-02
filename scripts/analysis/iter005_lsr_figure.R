#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter005_lsr_figure.R
#
# Two Iteration 005 figures: PIP calibration and false discovery control, each
# split across the six annotation-to-causality relationships.
#
#   Rscript scripts/analysis/iter005_lsr_figure.R <out_dir> [l3_rds]
#
# Layout is relationship down the rows and annotation type across the columns,
# which gives each panel roughly twice the width it had transposed. Annotation
# type is never pooled with the relationships - they are separate strata and
# pooling has been shown to invert conclusions.
#
# The null row comes first and is shaded: there the annotations are supplied but
# carry NO information about causality, so it is the arm against which every
# other row should be read.
#
# Calibration rests on the four PIP bands frozen into iter004_collect.R
# (lo/mid/hi/top). Finer bands would need the raw PIPs, which for Iteration 005
# still live only in results.rds on ephemeral.
# =============================================================================

.a  <- commandArgs(TRUE)
FIG <- if (length(.a) >= 1) .a[1] else "figures"
L3F <- if (length(.a) >= 2) .a[2] else "results/iter005/combined_scenario_metrics.rds"
suppressMessages(library(ggplot2))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

DROP  <- c("carma", "fb_pooled", "polyfun_est")
REL   <- c("null", "additive", "threshold", "mixed", "cooccur", "nonmono")
BANDS <- c("lo", "mid", "hi", "top")
TH    <- c("50", "80", "90", "95", "99"); TV <- c(.50, .80, .90, .95, .99)
MIN_N <- 50L
REFC  <- "#C1272D"
GREY  <- "#9AA0A6"

# sbayesrc IS included here. Its reported quantity is mixture membership rather
# than a fine-mapping PIP, which is exactly why it belongs on a calibration
# figure: the gap it opens is the visible consequence of that difference.
HI <- c(susie               = "#0072B2",
        polyfun_ldsc        = "#009E73",
        funmap              = "#E69F00",
        fb_xregion          = "#D55E00",
        functional_beatrice = "#CC79A7",
        sbayesrc            = "#56B4E9",
        polyfun_oracle      = "#000000")

L3 <- readRDS(L3F)
L3 <- L3[L3$n_failed < L3$n_fits & !L3$method %in% DROP, ]
L3$rel <- factor(L3$relationship, REL)
L3$arm <- factor(L3$annotation_type, c("binary", "continuous"),
                 c("Binary annotations", "Continuous annotations"))

shade <- data.frame(rel = factor("null", REL))

th <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey70", linewidth = .35),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 7.2, face = "bold"),
        axis.title = element_text(size = 7.8), axis.text = element_text(size = 6.6),
        legend.text = element_text(size = 6.9), legend.key.size = unit(8, "pt"),
        legend.position = "bottom",
        legend.box.margin = margin(t = -4),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7.2, colour = "grey35"))

grp <- function(m) factor(ifelse(m %in% names(HI), as.character(m), "other methods"),
                          c(names(HI), "other methods"))
pal <- c(HI, `other methods` = GREY)

# ---------------------------------------------------------------------------
# Figure: PIP calibration. Counts pooled to the stratum before any rate is
# formed - a reliability estimate inside one cell rests on too few variants.
# ---------------------------------------------------------------------------
cal <- do.call(rbind, lapply(split(L3, list(L3$arm, L3$rel, L3$method), drop = TRUE),
  function(r) do.call(rbind, lapply(BANDS, function(b) {
    n <- sum(r[[paste0("n_band_", b)]]); c_ <- sum(r[[paste0("c_band_", b)]])
    sp <- sum(r[[paste0("sum_pip_band_", b)]])
    if (!is.finite(n) || n < MIN_N) return(NULL)
    data.frame(arm = r$arm[1], rel = r$rel[1], method = as.character(r$method[1]),
               x = sp / n, y = c_ / n,
               lo = qbeta(.025, c_ + .5, n - c_ + .5),
               hi = qbeta(.975, c_ + .5, n - c_ + .5))}))))
cal$grp <- grp(cal$method)
cal <- cal[order(cal$grp == "other methods", decreasing = TRUE), ]

p1 <- ggplot(cal, aes(x, y, group = method, colour = grp)) +
  geom_rect(data = shade, inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            fill = "grey50", alpha = .07) +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              colour = REFC, linewidth = .4) +
  geom_line(data = ~subset(.x, grp == "other methods"), linewidth = .3, alpha = .8) +
  geom_line(data = ~subset(.x, grp != "other methods"), linewidth = .55) +
  geom_linerange(data = ~subset(.x, grp != "other methods"),
                 aes(ymin = lo, ymax = hi), linewidth = .35) +
  geom_point(data = ~subset(.x, grp != "other methods"), size = 1.1) +
  facet_grid(rel ~ arm) +
  scale_colour_manual(values = pal, limits = names(pal), name = NULL) +
  scale_x_continuous(breaks = c(0, .5, 1)) + scale_y_continuous(breaks = c(0, .5, 1)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean PIP assigned within band", y = "Proportion of those variants causal",
       title = "PIP calibration under each annotation-to-causality relationship",
       subtitle = paste0("Red dashed line is y = x; below it is overconfident.\n",
                         "Counts pooled within stratum, bars are 95% Jeffreys intervals.\n",
                         "The shaded null row carries no annotation-causality relationship.")) + th
ggsave(file.path(FIG, "fig_rel_calibration.pdf"), p1, width = 6.6, height = 9.0)

# ---------------------------------------------------------------------------
# Figure: false discovery control. Rates formed within a cell then averaged;
# cells selecting nothing are excluded rather than scored zero.
# ---------------------------------------------------------------------------
fdr <- do.call(rbind, lapply(split(L3, list(L3$arm, L3$rel, L3$method), drop = TRUE),
  function(r) do.call(rbind, lapply(seq_along(TH), function(i) {
    k <- r[[paste0("nsel_at_", TH[i])]] > 0
    v <- r[[paste0("fdr_at_", TH[i])]][k]; v <- v[is.finite(v)]
    if (length(v) < 2L) return(NULL)
    data.frame(arm = r$arm[1], rel = r$rel[1], method = as.character(r$method[1]),
               pos = i, m = mean(v), se = sd(v) / sqrt(length(v)))}))))
fdr$grp <- grp(fdr$method)
fdr <- fdr[order(fdr$grp == "other methods", decreasing = TRUE), ]
bound <- data.frame(pos = seq_along(TV), m = 1 - TV)

p2 <- ggplot(fdr, aes(pos, m, group = method, colour = grp)) +
  geom_rect(data = shade, inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            fill = "grey50", alpha = .07) +
  geom_line(data = bound, aes(pos, m), inherit.aes = FALSE,
            colour = REFC, linetype = "22", linewidth = .4) +
  geom_line(data = ~subset(.x, grp == "other methods"), linewidth = .3, alpha = .8) +
  geom_line(data = ~subset(.x, grp != "other methods"), linewidth = .55) +
  geom_linerange(data = ~subset(.x, grp != "other methods"),
                 aes(ymin = pmax(0, m - 2 * se), ymax = m + 2 * se), linewidth = .3) +
  geom_point(data = ~subset(.x, grp != "other methods"), size = 1.1) +
  facet_grid(rel ~ arm) +
  scale_colour_manual(values = pal, limits = names(pal), name = NULL) +
  scale_x_continuous(breaks = seq_along(TV),
                     labels = c(".50", ".80", ".90", ".95", ".99")) +
  labs(x = "PIP threshold t", y = "False discovery rate (mean over cells, ±2 SE)",
       title = "False discovery control under each annotation-to-causality relationship",
       subtitle = paste0("Red dashed line is y = 1 - t, the highest rate a calibrated ",
                         "method can attain at that threshold.\nBars are 2 SE over four cells ",
                         "and wide by design; cells selecting nothing are excluded.\n",
                         "The shaded null row carries no relationship at all.")) + th
ggsave(file.path(FIG, "fig_rel_fdr.pdf"), p2, width = 6.6, height = 9.0)

cat("fig_rel_calibration and fig_rel_fdr written\n")
