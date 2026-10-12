#!/usr/bin/env Rscript
# =============================================================================
# fig1_auprc.R  --  Figure 1 (average precision, six strata), revised.
#
#   Rscript fig1_auprc.R <out_pdf> <iter004_L3> [iter008_L3]
#
# Three changes against iter004_lsr_figures.R (figure 1):
#   1. Overarching strip labels: "Annotation type" above the column strips and
#      "Generative model" beside the row strips, added by gtable surgery
#      (ggh4x is not installed).
#   2. Only the seven displayed methods are plotted; the grey "other methods"
#      group is dropped.
#   3. Methods are labelled with their published names.
#
# If an iter008 L3 table is given, polyfun_ldsc, polyfun_oracle and sbayesrc
# are taken from it and the remaining methods from iter004. The two runs are
# exactly paired: susie's per-cell AP is identical in both.
#
# Everything else follows iter004_lsr_figures.R: the same L3 table, the same
# exclusions, the mean over cells with a 95% interval, the lsr_palette.R
# colours and shapes, and the ordering by the sparse / continuous stratum.
# =============================================================================

suppressMessages({library(ggplot2); library(grid); library(gtable)})
.a   <- commandArgs(TRUE)
OUT  <- if (length(.a) >= 1) .a[1] else "fig_auprc.pdf"
L3F  <- if (length(.a) >= 2) .a[2] else "results/iter004/combined_scenario_metrics_with_power.rds"
L3F8 <- if (length(.a) >= 3) .a[3] else NA_character_
source("scripts/analysis/lsr_palette.R")

DROP <- c("carma", "fb_pooled", "polyfun_est")
NEW  <- c("polyfun_ldsc", "polyfun_oracle", "sbayesrc")

L3 <- readRDS(L3F)
L3 <- L3[!L3$method %in% DROP & L3$n_failed < L3$n_fits, ]

if (!is.na(L3F8)) {                      # swap in the paper-faithful comparators
  L8 <- readRDS(L3F8)
  L8 <- L8[L8$method %in% NEW & L8$n_failed < L8$n_fits, ]
  keep <- intersect(names(L3), names(L8))
  L3 <- rbind(L3[!L3$method %in% NEW, keep], L8[, keep])
  cat("swapped in from iter008:", paste(NEW, collapse = ", "), "\n")
}

L3 <- L3[L3$method %in% M7, ]
stopifnot(all(M7 %in% L3$method))
L3$model <- factor(L3$model, c("sparse", "sparse_inf"),
                   c("Sparse", "Sparse + infinitesimal"))
L3$arm   <- factor(L3$annotation_type, c("none", "binary", "continuous"),
                   c("No annotations", "Binary annotations", "Continuous annotations"))

## ---- published names ------------------------------------------------------
NAMES <- c(susie          = "SuSiE",
           polyfun_ldsc   = "PolyFun + SuSiE",
           finemap_inf    = "FINEMAP-inf",
           polyfun_oracle = "Oracle PolyFun + SuSiE",
           sbayesrc       = "SBayesRC",
           beatrice       = "BEATRICE",
           fb_xregion     = "Functional BEATRICE")
stopifnot(setequal(names(NAMES), M7))

## ---- mean over cells, as in iter004_lsr_figures.R ------------------------
agg <- function(d, f) do.call(rbind, lapply(
  split(d, list(d$model, d$arm, d$method), drop = TRUE), f))
ap <- agg(L3, function(r) {
  v <- r$ap[is.finite(r$ap)]; if (!length(v)) return(NULL)
  data.frame(model = r$model[1], arm = r$arm[1], method = as.character(r$method[1]),
             m = mean(v), se = sd(v) / sqrt(length(v)), n = length(v))})

ord <- ap[ap$model == "Sparse" & ap$arm == "Continuous annotations", ]
ord <- ord$method[order(ord$m)]
ap$method <- factor(NAMES[ap$method], levels = NAMES[ord])
pal <- setNames(PAL7[M7], NAMES[M7])
shp <- setNames(SHP7[M7], NAMES[M7])
lim <- unname(NAMES[M7])          # unnamed, or ggplot labels the legend with the names

th <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linewidth = .25, colour = "grey90"),
        panel.border = element_rect(colour = "grey70", linewidth = .35),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 6.6, face = "bold"),
        axis.title = element_text(size = 6.8),
        axis.text  = element_text(size = 6.2),
        legend.position = "none")

p1 <- ggplot(ap, aes(m, method, colour = method, shape = method)) +
  geom_errorbarh(aes(xmin = m - 1.96 * se, xmax = m + 1.96 * se),
                 height = .42, linewidth = .45) +
  geom_point(size = 1.6, stroke = .55) +
  facet_grid(model ~ arm) +
  scale_colour_manual(values = pal, limits = lim, name = NULL) +
  scale_shape_manual(values = shp, limits = lim, name = NULL) +
  scale_x_continuous(breaks = seq(.3, .9, .2)) +
  labs(x = "Mean average precision over cells (95% CI)", y = NULL) +
  th + theme(panel.grid.major.y = element_line(linewidth = .2, colour = "grey93"))

## ---- overarching strip labels --------------------------------------------
# A row above the column strips and a column outside the row strips, each
# holding one spanning label drawn in the strip's own style.
outer_strips <- function(p, top, right, fill = "grey94", size = 6.6) {
  g  <- ggplotGrob(p)
  lab <- function(txt, rot = 0) grobTree(
    rectGrob(gp = gpar(fill = fill, col = NA)),
    textGrob(txt, rot = rot, gp = gpar(fontsize = size, fontface = "bold")))

  st <- min(g$layout$t[grepl("^strip-t", g$layout$name)])
  pc <- range(g$layout$l[grepl("^panel", g$layout$name)])
  g  <- gtable_add_rows(g, unit(1.15, "lines"), pos = st - 1)
  g  <- gtable_add_grob(g, lab(top), t = st, l = pc[1], r = pc[2], name = "outer-top")

  sr <- max(g$layout$r[grepl("^strip-r", g$layout$name)])
  pr <- range(g$layout$t[grepl("^panel", g$layout$name)])
  g  <- gtable_add_cols(g, unit(1.15, "lines"), pos = sr)
  g  <- gtable_add_grob(g, lab(right, rot = -90),
                        t = pr[1], b = pr[2], l = sr + 1, name = "outer-right")
  g
}
g1 <- outer_strips(p1, "Annotation type", "Generative model")

if (grepl("\\.png$", OUT)) png(OUT, width = 7.28, height = 3.0, units = "in", res = 200, type = "cairo") else PDF_DEVICE(OUT, width = 7.28, height = 3.0)
grid.draw(g1)
invisible(dev.off())
cat("wrote", OUT, "\n")

## ---- the plotted numbers, for the text -----------------------------------
tab <- reshape(ap[, c("model", "arm", "method", "m")],
               idvar = c("model", "method"), timevar = "arm", direction = "wide")
names(tab) <- sub("^m\\.", "", names(tab))
tab[, 3:5] <- round(tab[, 3:5], 3)
print(tab[order(tab$model, tab$method), ], row.names = FALSE)
