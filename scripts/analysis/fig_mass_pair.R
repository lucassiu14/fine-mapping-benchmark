#!/usr/bin/env Rscript
# =============================================================================
# fig_mass_pair.R  --  where the posterior mass goes, and how much there is.
#
#   Rscript fig_mass_pair.R <out_pdf> <mass_iter004> <mass_iter008>
#
# Two quantities per fit, from the PIP tails:
#
#   mean PIP of a causal variant = (sum of PIPs over causal variants) / S
#       bounded in [0,1]; what probability the truth actually receives.
#   total mass / S               = (sum of PIPs over all variants) / S
#       1 when the declared mass matches the truth, >1 when it is inflated.
#
# Together they separate what the single ratio causal_mass/total_mass mixes:
# a method is not penalised here for being diffuse in the first panel, and
# cannot hide inflation in the second.
#
# Rows are the number of causal variants; columns nest metric inside generative
# model. Methods that use annotations are read on the binary and continuous
# strata together, those that cannot on the unannotated stratum. Fits are
# aggregated to cells and each point is the mean over cells with 1.96 SE.
# =============================================================================

suppressMessages({library(ggplot2); library(grid); library(gtable)})
.a  <- commandArgs(TRUE)
OUT <- .a[1]
PAL_F <- Sys.getenv("LSR_PALETTE", "scripts/analysis/lsr_palette.R")
if (!file.exists(PAL_F)) stop("cannot find lsr_palette.R at '", PAL_F, "'")
source(PAL_F)

d <- readRDS(.a[2])
if (length(.a) >= 3 && file.exists(.a[3])) {
  d8 <- readRDS(.a[3]); d <- rbind(d[!d$method %in% unique(d8$method), ], d8)
}
ANNOT <- c("polyfun_oracle", "polyfun_ldsc", "sbayesrc", "fb_xregion")
BLIND <- c("susie", "finemap_inf", "beatrice")
d <- rbind(d[d$method %in% ANNOT & d$annotation_type %in% c("binary", "continuous"), ],
           d[d$method %in% BLIND & d$annotation_type == "none", ])

d$causal_pip <- d$causal_mass / d$S
d$mass_ratio <- d$total_mass  / d$S
d$cell <- paste(d$job_dir, d$S, d$phi, d$region_size, sep = "|")

long <- rbind(
  data.frame(d[, c("model", "S", "method", "cell")], metric = "Mean PIP of a causal variant",
             v = d$causal_pip),
  data.frame(d[, c("model", "S", "method", "cell")], metric = "Posterior mass ratio",
             v = d$mass_ratio))
cells <- aggregate(v ~ model + S + method + metric + cell, long, mean)
out <- do.call(rbind, lapply(
  split(cells, list(cells$model, cells$S, cells$method, cells$metric), drop = TRUE),
  function(r) data.frame(model = r$model[1], S = r$S[1], method = r$method[1],
                         metric = r$metric[1], m = mean(r$v),
                         se = sd(r$v) / sqrt(nrow(r)))))

NAMES <- c(susie = "SuSiE", polyfun_ldsc = "PolyFun + SuSiE", finemap_inf = "FINEMAP-inf",
           polyfun_oracle = "Oracle PolyFun + SuSiE", sbayesrc = "SBayesRC",
           beatrice = "BEATRICE", fb_xregion = "Functional BEATRICE")
have_a <- ANNOT[ANNOT %in% out$method]; have_b <- BLIND[BLIND %in% out$method]
lev <- rev(c(NAMES[have_a], " ", NAMES[have_b]))   # " " is the gap between the groups
out$lab  <- factor(NAMES[out$method], levels = lev)
out$Slab <- factor(out$S, sort(unique(out$S)), sprintf("S = %s", sort(unique(out$S))))
out$mlab <- factor(out$model, c("sparse", "sparse_inf"), c("Sparse", "Sparse + infinitesimal"))
out$metric <- factor(out$metric, c("Mean PIP of a causal variant",
                                   "Posterior mass ratio"))
pal <- setNames(PAL7[names(NAMES)], NAMES)

# The mass ratio spans 0.5 to 21, so that column is drawn on a log2 axis with
# its baseline at 1: a bar to the right over-declares mass, to the left under-
# declares. The PIP column stays linear on [0,1]. Break and label functions
# switch on the column's own range, since facet_grid shares one scale object.
MASS <- levels(out$metric)[2]
tx <- function(v, metric) ifelse(metric == MASS, log2(pmax(v, 1e-6)), v)
out$mp <- tx(out$m, out$metric)
out$lo <- tx(pmax(out$m - 1.96 * out$se, 1e-6), out$metric)
out$hi <- tx(out$m + 1.96 * out$se, out$metric)

RT <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32)
brk <- function(l) if (max(l, na.rm = TRUE) <= 1.2) seq(0, 1, .5) else log2(RT)
lbl <- function(b) if (max(b, na.rm = TRUE) <= 1.2) format(b) else
  format(2^b, drop0trailing = TRUE)

rng <- range(c(out$lo, out$hi)[out$metric == MASS], na.rm = TRUE)
lims <- do.call(rbind, lapply(levels(out$mlab), function(ml) rbind(
  data.frame(mlab = ml, metric = levels(out$metric)[1], x = c(0, 1)),
  data.frame(mlab = ml, metric = MASS, x = rng))))
lims$Slab <- factor(levels(out$Slab)[1], levels(out$Slab))
lims$mlab <- factor(lims$mlab, levels(out$mlab))
lims$metric <- factor(lims$metric, levels(out$metric))
lims$lab <- factor(levels(out$lab)[1], levels(out$lab))
vline <- data.frame(mlab = factor(levels(out$mlab), levels(out$mlab)),
                    metric = factor(MASS, levels(out$metric)))

th <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = .25, colour = "grey90"),
        panel.border = element_rect(colour = "grey70", linewidth = .35),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 6.5, face = "bold"),
        axis.title = element_text(size = 6.8), axis.text = element_text(size = 6.1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 5.8),
        legend.key.size = unit(6.5, "pt"), legend.box.margin = margin(t = -4),
        legend.margin = margin(t = 0, b = 0), legend.spacing.x = unit(1, "pt"),
        legend.key.spacing.x = unit(3, "pt"),
        panel.spacing.x = unit(5, "pt"))

p <- ggplot(out, aes(mp, lab, fill = lab)) +
  geom_vline(data = vline, aes(xintercept = 0), inherit.aes = FALSE,
             colour = REF, linetype = "22", linewidth = .4) +
  geom_col(width = .68) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .28, linewidth = .3,
                 colour = "grey25") +
  geom_blank(data = lims, aes(x = x, y = lab), inherit.aes = FALSE) +
  facet_grid(Slab ~ mlab + metric, scales = "free_x") +
  scale_fill_manual(values = pal, na.value = NA, name = NULL,
                    breaks = unname(c(NAMES[have_a], NAMES[have_b])),
                    guide = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = brk, labels = lbl, expand = expansion(mult = c(.03, .03))) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) + th

g <- ggplotGrob(p)
lab <- function(txt, rot = 0) grobTree(
  rectGrob(gp = gpar(fill = "grey94", col = NA)),
  textGrob(txt, rot = rot, gp = gpar(fontsize = 6.5, fontface = "bold")))

# facet_grid draws the model strip once per column; replace that row with two
# grobs, each spanning the two metric columns of its model.
pl  <- sort(unique(g$layout$l[grepl("^panel", g$layout$name)]))
top <- min(g$layout$t[grepl("^strip-t", g$layout$name)])
keep <- !(grepl("^strip-t", g$layout$name) & g$layout$t == top)
g$grobs <- g$grobs[keep]; g$layout <- g$layout[keep, ]
g <- gtable_add_grob(g, lab(levels(out$mlab)[1]), t = top, l = pl[1], r = pl[2], name = "m1")
g <- gtable_add_grob(g, lab(levels(out$mlab)[2]), t = top, l = pl[3], r = pl[4], name = "m2")

# this ggplot2 stacks both faceting variables into one strip row, so the metric
# labels went with it; add them back underneath, one per column.
g <- gtable_add_rows(g, unit(1.5, "lines"), pos = top)
for (i in seq_along(pl))
  g <- gtable_add_grob(g, lab(levels(out$metric)[(i - 1) %% 2 + 1]),
                       t = top + 1, l = pl[i], name = paste0("k", i))

sr <- max(g$layout$r[grepl("^strip-r", g$layout$name)])
pr <- range(g$layout$t[grepl("^panel", g$layout$name)])
g  <- gtable_add_cols(g, unit(1.15, "lines"), pos = sr)
g  <- gtable_add_grob(g, lab("Causal variants per region", rot = -90),
                      t = pr[1], b = pr[2], l = sr + 1, name = "or")

if (grepl("\\.png$", OUT)) png(OUT, width = 7.28, height = 3.8, units = "in", res = 200, type = "cairo") else PDF_DEVICE(OUT, width = 7.28, height = 3.8)
grid.draw(g); invisible(dev.off())
cat("wrote", OUT, "\n")

w <- out[out$model == "sparse", c("S", "method", "metric", "m")]
w <- reshape(w, idvar = c("method", "metric"), timevar = "S", direction = "wide")
names(w) <- sub("^m\\.", "S=", names(w))
cat("\nSparse model:\n"); print(w[order(w$metric, w$method), ], row.names = FALSE, digits = 2)
