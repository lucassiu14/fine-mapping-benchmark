#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter005_lsr_figure.R
#
# The Iteration 005 figure: what each annotation-aware method gains over
# annotation-blind susie, as the SHAPE of the annotation -> causality
# relationship changes.
#
#   Rscript scripts/analysis/iter005_lsr_figure.R <out_dir> [l3_rds]
#
# Everything is a PAIRED difference against susie on the same cells. The level
# of average precision is not comparable between relationships - the six forms
# make the problem different amounts of hard, even after the concentration
# calibration - but the gap over an annotation-blind baseline is.
#
# The null arm is drawn first and deliberately: it is the arm where the
# annotations carry NO information about causality, so every gap there should
# be zero. A method whose line starts above zero has an advantage that is not
# annotation information, and its gains elsewhere have to be read net of it.
# =============================================================================

.a  <- commandArgs(TRUE)
FIG <- if (length(.a) >= 1) .a[1] else "figures"
L3F <- if (length(.a) >= 2) .a[2] else "results/iter005/combined_scenario_metrics.rds"
suppressMessages(library(ggplot2))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

DROP <- c("carma", "fb_pooled", "polyfun_est")
REL  <- c("null", "additive", "threshold", "mixed", "cooccur", "nonmono")
KEY  <- c("S", "phi", "region_size")     # job_dir differs between relationships

L3 <- readRDS(L3F)
L3 <- L3[L3$n_failed < L3$n_fits & !L3$method %in% DROP, ]

gap <- function(m, arm, rel) {
  d <- L3[L3$annotation_type == arm & L3$relationship == rel, ]
  x <- d[d$method == m, ]; y <- d[d$method == "susie", ]
  i <- match(do.call(paste, c(x[KEY], sep = "\r")),
             do.call(paste, c(y[KEY], sep = "\r")))
  v <- x$ap[!is.na(i)] - y$ap[i[!is.na(i)]]; v <- v[is.finite(v)]
  if (length(v) < 2L) return(c(NA_real_, NA_real_))
  c(mean(v), sd(v) / sqrt(length(v)))
}

METH <- setdiff(sort(unique(L3$method)), "susie")
d <- do.call(rbind, lapply(METH, function(m) do.call(rbind, lapply(
  c("binary", "continuous"), function(a) do.call(rbind, lapply(REL, function(r) {
    g <- gap(m, a, r)
    data.frame(method = m, arm = a, rel = r, m = g[1], se = g[2])}))))))
d <- d[is.finite(d$m), ]
d$rel <- factor(d$rel, REL)
d$arm <- factor(d$arm, c("binary", "continuous"),
                c("Binary annotations", "Continuous annotations"))

HI <- c(polyfun_oracle = "#000000", polyfun_ldsc = "#009E73", funmap = "#E69F00",
        fb_xregion = "#D55E00", functional_beatrice = "#CC79A7")
d$grp <- ifelse(d$method %in% names(HI), d$method, "other methods")
d$grp <- factor(d$grp, c(names(HI), "other methods"))
d <- d[order(d$grp == "other methods", decreasing = TRUE), ]   # focus drawn last

p <- ggplot(d, aes(rel, m, colour = grp, group = method)) +
  annotate("rect", xmin = .4, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = "grey50", alpha = .07) +
  geom_hline(yintercept = 0, linetype = "22", colour = "#C1272D", linewidth = .45) +
  geom_line(data = ~subset(.x, grp == "other methods"), linewidth = .32, alpha = .8) +
  geom_line(data = ~subset(.x, grp != "other methods"), linewidth = .6) +
  geom_errorbar(data = ~subset(.x, grp != "other methods"),
                aes(ymin = m - 2 * se, ymax = m + 2 * se), width = .13, linewidth = .4) +
  geom_point(data = ~subset(.x, grp != "other methods"), size = 1.5) +
  facet_wrap(~ arm, nrow = 1) +
  # sbayesrc sits near -0.2 in every arm and would compress everything else into
  # the top fifth of the panel; it is named in the caption instead.
  coord_cartesian(ylim = c(-0.09, 0.17)) +
  scale_colour_manual(values = c(HI, `other methods` = "#9AA0A6"),
                      limits = c(names(HI), "other methods"), name = NULL) +
  labs(x = NULL, y = "Average precision, gain over annotation-blind susie",
       title = "What annotations buy, as the shape of the relationship changes",
       subtitle = paste("Paired differences on the same cells, ±2 SE. The shaded null arm",
                        "carries no annotation-causality link,\nso every method should sit on",
                        "the red line there. Iteration 005: four cells per point.",
                        "sbayesrc is off scale near -0.19.")) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "grey70", linewidth = .35),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 7.4, face = "bold"),
        axis.title = element_text(size = 7.8), axis.text = element_text(size = 6.9),
        legend.text = element_text(size = 6.9), legend.key.size = unit(8, "pt"),
        legend.position = "right",
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7.2, colour = "grey35"))
ggsave(file.path(FIG, "fig_relationships.pdf"), p, width = 9.2, height = 3.9)
cat("fig_relationships written\n")
