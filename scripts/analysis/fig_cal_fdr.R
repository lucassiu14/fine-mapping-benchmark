#!/usr/bin/env Rscript
# =============================================================================
# fig_cal_fdr.R  --  PIP calibration and false discovery rate, six strata,
# per-iteration (draw) estimator, with the iteration 008 comparators.
#
#   Rscript fig_cal_fdr.R <out_dir> <periter_iter004> <periter_iter008>
#
# polyfun_ldsc, polyfun_oracle and sbayesrc are taken from the iteration 008
# tables; susie, finemap_inf, beatrice and fb_xregion from iteration 004. The
# two runs are exactly paired (susie identical fit for fit).
#
# Three changes against iter004_lsr_figures.R, matching Figure 1:
#   1. Overarching strip labels, by gtable surgery (ggh4x is not installed).
#   2. Only the seven displayed methods (these figures already showed seven).
#   3. Published method names.
# Everything else - the estimator, the reportability filter, the reference
# lines, the palette - is unchanged.
# =============================================================================

suppressMessages({library(ggplot2); library(grid); library(gtable)})
.a  <- commandArgs(TRUE)
FIG <- .a[1]; P4 <- .a[2]; P8 <- .a[3]
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
PAL_F <- Sys.getenv("LSR_PALETTE", "scripts/analysis/lsr_palette.R")
if (!file.exists(PAL_F)) stop("cannot find lsr_palette.R at '", PAL_F, "'")
source(PAL_F)

NEW <- c("polyfun_ldsc", "polyfun_oracle", "sbayesrc")
p4 <- readRDS(P4); p8 <- readRDS(P8)
stopifnot(isTRUE(p4$meta$validated), isTRUE(p8$meta$validated))
merge_tab <- function(a, b) {
  k <- intersect(names(a), names(b))
  rbind(a[!a$method %in% NEW, k], b[b$method %in% NEW, k])
}
cal <- merge_tab(p4$cal, p8$cal)
fdr <- merge_tab(p4$fdr, p8$fdr)
cal <- cal[cal$estimator == "draw" & cal$reportable, ]
fdr <- fdr[fdr$estimator == "draw" & fdr$reportable, ]
cat("calibration rows:", nrow(cal), " fdr rows:", nrow(fdr), "\n")
cat("methods:", paste(sort(unique(cal$method)), collapse = ", "), "\n")

NAMES <- c(susie = "SuSiE", polyfun_ldsc = "PolyFun + SuSiE", finemap_inf = "FINEMAP-inf",
           polyfun_oracle = "Oracle PolyFun + SuSiE", sbayesrc = "SBayesRC",
           beatrice = "BEATRICE", fb_xregion = "Functional BEATRICE")
pal <- setNames(PAL7[M7], NAMES[M7]); shp <- setNames(SHP7[M7], NAMES[M7])
lim <- unname(NAMES[M7])
prep <- function(d) {
  d$model <- factor(as.character(d$model), c("sparse", "sparse_inf"),
                    c("Sparse", "Sparse + infinitesimal"))
  d$arm <- factor(as.character(d$annotation_type), c("none", "binary", "continuous"),
                  c("No annotations", "Binary annotations", "Continuous annotations"))
  d$method <- factor(NAMES[as.character(d$method)], levels = lim)
  d
}
cal <- prep(cal); fdr <- prep(fdr)

th <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linewidth = .25, colour = "grey90"),
        panel.border = element_rect(colour = "grey70", linewidth = .35),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 6.6, face = "bold"),
        axis.title = element_text(size = 6.8), axis.text = element_text(size = 6.2),
        legend.text = element_text(size = 6.4), legend.key.size = unit(8, "pt"),
        legend.position = "bottom", legend.box.margin = margin(t = -4),
        legend.margin = margin(t = 0, b = 0))
msc <- list(scale_colour_manual(values = pal, limits = lim, name = NULL,
                                guide = guide_legend(nrow = 1)),
            scale_shape_manual(values = shp, limits = lim, name = NULL,
                               guide = guide_legend(nrow = 1)))

outer_strips <- function(p, top, right, fill = "grey94", size = 6.6) {
  g <- ggplotGrob(p)
  lab <- function(txt, rot = 0) grobTree(
    rectGrob(gp = gpar(fill = fill, col = NA)),
    textGrob(txt, rot = rot, gp = gpar(fontsize = size, fontface = "bold")))
  st <- min(g$layout$t[grepl("^strip-t", g$layout$name)])
  pc <- range(g$layout$l[grepl("^panel", g$layout$name)])
  g  <- gtable_add_rows(g, unit(1.15, "lines"), pos = st - 1)
  g  <- gtable_add_grob(g, lab(top), t = st, l = pc[1], r = pc[2], name = "ot")
  sr <- max(g$layout$r[grepl("^strip-r", g$layout$name)])
  pr <- range(g$layout$t[grepl("^panel", g$layout$name)])
  g  <- gtable_add_cols(g, unit(1.15, "lines"), pos = sr)
  g  <- gtable_add_grob(g, lab(right, rot = -90), t = pr[1], b = pr[2], l = sr + 1, name = "or")
  g
}
save_g <- function(p, file, w = 7.28, h = 3.6) {
  g <- outer_strips(p, "Annotation type", "Generative model")
  f <- file.path(FIG, file)
  if (nzchar(Sys.getenv("AS_PNG"))) { f <- sub("\\.pdf$", ".png", f)
    png(f, width = w, height = h, units = "in", res = 200, type = "cairo") }
  else PDF_DEVICE(f, width = w, height = h)
  grid.draw(g); invisible(dev.off())
  cat("wrote", file.path(FIG, file), "\n")
}

## --- calibration ----------------------------------------------------------
cal <- cal[order(cal$method, cal$x), ]
p_cal <- ggplot(cal, aes(x, y, group = method, colour = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = REF, linewidth = .45) +
  geom_line(linewidth = .5) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .35, na.rm = TRUE) +
  geom_point(size = 1.2, stroke = .5) +
  facet_grid(model ~ arm) + msc +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean PIP assigned within band",
       y = "Proportion of those variants causal") + th
save_g(p_cal, "fig_calibration.pdf")

## --- false discovery rate -------------------------------------------------
TV <- sort(unique(fdr$t)); fdr$pos <- match(fdr$t, TV)
bound <- data.frame(pos = seq_along(TV), m = 1 - TV)
p_fdr <- ggplot(fdr, aes(pos, m, group = method, colour = method, shape = method)) +
  geom_line(data = bound, aes(pos, m), inherit.aes = FALSE,
            colour = REF, linetype = "22", linewidth = .45) +
  geom_ribbon(aes(ymin = m - 2 * se, ymax = m + 2 * se, fill = method),
              colour = NA, alpha = .12, show.legend = FALSE, na.rm = TRUE) +
  geom_line(linewidth = .45) + geom_point(size = 1.2, stroke = .5) +
  facet_grid(model ~ arm) + msc +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(breaks = seq_along(TV), labels = sub("^0", "", sprintf("%.2f", TV))) +
  labs(x = "PIP threshold t",
       y = "False discovery rate (mean over iterations, ±2 SE)") + th
save_g(p_fdr, "fig_fdr.pdf")

## --- the numbers the text quotes -----------------------------------------
f <- fdr[fdr$t %in% c(0.9, 0.95), c("model", "arm", "method", "t", "m")]
f$m <- round(f$m, 3)
cat("\nFalse discovery rate at t = 0.9 and 0.95:\n")
print(reshape(f, idvar = c("model", "arm", "method"), timevar = "t", direction = "wide"),
      row.names = FALSE)
