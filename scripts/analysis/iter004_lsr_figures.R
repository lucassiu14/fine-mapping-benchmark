#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_lsr_figures.R
#
# The five figures of the LSR results section (§5.3-5.6), all six strata.
#
#   Rscript scripts/analysis/iter004_lsr_figures.R <out_dir> [l3_rds] [cal9_rds]
#
# Every figure is a 2x3 facet grid - generative model down the rows, annotation
# regime across the columns - because model and annotation type are never
# pooled. See docs/evaluation.md.
#
# Excluded throughout: carma, fb_pooled, polyfun_est. polyfun_ldsc is the
# canonical PolyFun (S-LDSC + leave-one-region-out); polyfun_est regresses
# chi^2 on each variant's own annotations and is documented as an error in
# wrapper_polyfun_ldsc.R.
#
# Cells with n_failed == n_fits are dropped: funmap has 625 such cells in the
# unannotated arm, with set_hit_95 == 0 rather than NA, and they would plot as
# coverage 0.000 instead of being absent.
#
# Calibration uses the nine bands from iter004_calibration_bands.R when that
# table is present, and falls back to the four bands frozen into
# iter004_collect.R when it is not.
# =============================================================================

suppressMessages({library(ggplot2); library(grid)})
.a  <- commandArgs(TRUE)
FIG <- if (length(.a) >= 1) .a[1] else "figures"
L3F <- if (length(.a) >= 2) .a[2] else "results/iter004/combined_scenario_metrics_with_power.rds"
CAL9<- if (length(.a) >= 3) .a[3] else "results/iter004/calibration_bands9.rds"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

DROP <- c("carma", "fb_pooled", "polyfun_est")

L3 <- readRDS(L3F)
stopifnot("power" %in% names(L3))   # add it with iter004_add_power.R
L3 <- L3[!L3$method %in% DROP, ]
L3 <- L3[L3$n_failed < L3$n_fits, ]
stopifnot(!any(DROP %in% L3$method))          # funmap/none is 625 all-failed cells
L3$model <- factor(L3$model, c("sparse","sparse_inf"),
                   c("Sparse","Sparse + infinitesimal"))
L3$arm   <- factor(L3$annotation_type, c("none","binary","continuous"),
                   c("No annotations","Binary annotations","Continuous annotations"))
FITS <- 20

## ---- Okabe-Ito highlights, identical across 5.4 / 5.5 / 5.6 ---------------
FOCUS <- c(susie="#0072B2", susie_inf="#56B4E9",
           fb_xregion="#D55E00", sbayesrc="#009E73")
GREY  <- "#9AA0A6"
REF   <- "#C1272D"          # reference / bound lines
focus_col <- function(m) ifelse(m %in% names(FOCUS), FOCUS[as.character(m)], GREY)

th <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linewidth = .25, colour = "grey90"),
        panel.border = element_rect(colour = "grey70", linewidth = .35),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 7.4, face = "bold"),
        axis.title = element_text(size = 7.8),
        axis.text  = element_text(size = 6.9),
        legend.title = element_text(size = 7.2),
        legend.text  = element_text(size = 6.9),
        legend.key.size = unit(8,"pt"),
        plot.title = element_text(size = 9, face = "bold"),
        plot.subtitle = element_text(size = 7.2, colour = "grey35"))

agg <- function(d, f) do.call(rbind, lapply(
  split(d, list(d$model, d$arm, d$method), drop = TRUE), f))

## =========================================================================
## Figure 1  -  average precision, six strata
## =========================================================================
ap <- agg(L3, function(r) {
  v <- r$ap[is.finite(r$ap)]; if (!length(v)) return(NULL)
  data.frame(model=r$model[1], arm=r$arm[1], method=as.character(r$method[1]),
             m=mean(v), se=sd(v)/sqrt(length(v)), n=length(v))})
ord <- ap[ap$model=="Sparse" & ap$arm=="Continuous annotations", ]
ord <- ord$method[order(ord$m)]
ap$method <- factor(ap$method, levels = ord)
HI1 <- c(polyfun_oracle="#000000", fb_xregion="#D55E00",
         functional_beatrice="#CC79A7", beatrice="#E69F00",
         susie="#0072B2", polyfun_ldsc="#009E73")
n_other <- length(setdiff(unique(as.character(ap$method)), names(HI1)))
OTH1 <- sprintf("other methods (%d)", n_other)
ap$grp <- ifelse(as.character(ap$method) %in% names(HI1),
                 as.character(ap$method), OTH1)
ap$grp <- factor(ap$grp, c(names(HI1), OTH1))

p1 <- ggplot(ap, aes(m, method, colour = grp)) +
  geom_errorbarh(aes(xmin = m-1.96*se, xmax = m+1.96*se),
                 height = .42, linewidth = .45) +
  geom_point(size = 1.75) +
  facet_grid(model ~ arm) +
  scale_colour_manual(values = setNames(c(HI1, GREY), c(names(HI1), OTH1)), name = NULL) +
  scale_x_continuous(breaks = seq(.3,.9,.2)) +
  labs(x = "Mean average precision over cells (95% CI)", y = NULL,
       title = "Average precision by method, six strata",
       subtitle = paste("Ordered by the sparse / continuous stratum.",
                        "polyfun_oracle reads the simulated truth and is a ceiling, not a competitor.")) +
  th + theme(legend.position = "right",
             panel.grid.major.y = element_line(linewidth = .2, colour = "grey93"))
ggsave(file.path(FIG,"fig_auprc.pdf"), p1, width = 9.2, height = 5.6)

## =========================================================================
## Figures 2 and 3  -  credible-set coverage, and power, six strata
## =========================================================================
cs <- agg(L3, function(r) data.frame(
  model=r$model[1], arm=r$arm[1], method=as.character(r$method[1]),
  cov = mean(r$set_hit_95, na.rm=TRUE)/FITS,
  size= mean(r$set_size_95, na.rm=TRUE),
  pow = mean(r$power, na.rm=TRUE)))
cs <- cs[is.finite(cs$size) & cs$method != "marginal_z", ]   # 740-variant sets
meths <- sort(unique(cs$method))
cs$id <- match(cs$method, meths)
cs$grp <- ifelse(cs$method %in% names(FOCUS), cs$method, "other methods")
cs$grp <- factor(cs$grp, c(names(FOCUS), "other methods"))
.lab <- sprintf("%d %s", seq_along(meths), meths)
.h   <- ceiling(length(.lab)/2)
KEY  <- paste(paste(.lab[1:.h], collapse = "   ·   "),
              paste(.lab[(.h+1):length(.lab)], collapse = "   ·   "), sep = "\n")

scatter <- function(yv, ylab, hline, title, sub) {
  g <- ggplot(cs, aes(size, .data[[yv]], colour = grp))
  if (!is.na(hline)) g <- g + geom_hline(yintercept = hline, linetype = 2,
                                         colour = "grey55", linewidth = .3)
  g + geom_text(aes(label = id), size = 2.25, fontface = "bold",
                show.legend = FALSE) +
    geom_point(alpha = 0, size = 1.6) +
    facet_grid(model ~ arm) +
    scale_colour_manual(values = c(FOCUS, `other methods` = GREY),
                        limits = c(names(FOCUS), "other methods"), name = NULL) +
    scale_x_log10(breaks = c(1,1.5,2,3,5,10,20,50),
                  labels = c("1","1.5","2","3","5","10","20","50"),
                  expand = expansion(mult = c(.06,.06))) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    labs(x = "Mean 95% credible-set size (variants, log scale)", y = ylab,
         title = title, subtitle = sub, caption = KEY) +
    th + theme(legend.position = "right",
               plot.caption = element_text(size = 6.2, colour = "grey30",
                                           hjust = 0, margin = margin(t = 6)))
}
p2 <- scatter("cov", "Coverage: sets containing a causal variant", 0.95,
  "Credible-set coverage against set size, six strata",
  "Dashed line is the nominal 0.95 level. marginal_z omitted: coverage 1.00 at ~740 variants per set.")
ggsave(file.path(FIG,"fig_coverage.pdf"), p2, width = 9.2, height = 5.8)

p3 <- scatter("pow", "Power: causal variants captured by a set", NA,
  "Credible-set power against set size, six strata",
  "Power is the proportion of causal variants falling in the set. marginal_z omitted: power 0.96 at ~740 variants per set.")
ggsave(file.path(FIG,"fig_power.pdf"), p3, width = 9.2, height = 5.8)

cat("figs 1-3 written\n")

## =========================================================================
## Figure 4  -  PIP calibration (counts pooled to stratum), six strata
## =========================================================================
if (file.exists(CAL9)) {
  message("calibration: nine bands from ", CAL9)
  cal <- readRDS(CAL9)
  cal <- cal[cal$reportable & !cal$method %in% DROP, ]
  cal$model <- factor(cal$model, c("sparse","sparse_inf"),
                      c("Sparse","Sparse + infinitesimal"))
  cal$arm   <- factor(cal$annotation_type, c("none","binary","continuous"),
                      c("No annotations","Binary annotations","Continuous annotations"))
  cal <- cal[order(cal$method, cal$x), ]
} else {
  message("calibration: falling back to the four frozen bands (no ", CAL9, ")")
  BANDS <- c("lo","mid","hi","top")
  cal <- do.call(rbind, lapply(split(L3, list(L3$model, L3$arm, L3$method), drop = TRUE),
    function(r) do.call(rbind, lapply(BANDS, function(b) {
      n <- sum(r[[paste0("n_band_",b)]], na.rm=TRUE)
      c_<- sum(r[[paste0("c_band_",b)]], na.rm=TRUE)
      s <- sum(r[[paste0("sum_pip_band_",b)]], na.rm=TRUE)
      if (!is.finite(n) || n < 50) return(NULL)
      lo <- qbeta(.025, c_+.5, n-c_+.5); hi <- qbeta(.975, c_+.5, n-c_+.5)
      data.frame(model=r$model[1], arm=r$arm[1], method=as.character(r$method[1]),
                 band=b, x=s/n, y=c_/n, lo=lo, hi=hi, n=n)}))))
}
cal$grp <- ifelse(cal$method %in% names(FOCUS), cal$method, "other methods")
cal$grp <- factor(cal$grp, c(names(FOCUS), "other methods"))
cal <- cal[order(cal$grp != "other methods", cal$method, cal$x), ]  # focus drawn last

p4 <- ggplot(cal, aes(x, y, group = method, colour = grp)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              colour = REF, linewidth = .45) +
  geom_line(data = subset(cal, grp == "other methods"), linewidth = .32, alpha = .75) +
  geom_line(data = subset(cal, grp != "other methods"), linewidth = .62) +
  geom_linerange(data = subset(cal, grp != "other methods"),
                 aes(ymin = lo, ymax = hi), linewidth = .42) +
  geom_point(data = subset(cal, grp != "other methods"), size = 1.35) +
  facet_grid(model ~ arm) +
  scale_colour_manual(values = c(FOCUS, `other methods` = GREY),
                      limits = c(names(FOCUS), "other methods"), name = NULL) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "Mean PIP assigned within band", y = "Proportion of those variants causal",
       title = "Reliability of reported posterior inclusion probabilities, six strata",
       subtitle = paste("Red dashed line is y = x, where a perfectly calibrated method sits.",
                        "Counts pooled to the stratum before any rate is formed;",
                        "bars are 95% Jeffreys intervals. Bands under 50 variants are not drawn.")) +
  th + theme(legend.position = "right")
ggsave(file.path(FIG,"fig_calibration.pdf"), p4, width = 9.2, height = 5.6)

## =========================================================================
## Figure 5  -  false discovery rate against threshold, six strata
## =========================================================================
TH <- c("50","80","90","95","99"); TV <- c(.50,.80,.90,.95,.99)
fdr <- do.call(rbind, lapply(split(L3, list(L3$model, L3$arm, L3$method), drop = TRUE),
  function(r) do.call(rbind, lapply(seq_along(TH), function(i) {
    keep <- r[[paste0("nsel_at_", TH[i])]] > 0        # abstaining cells excluded
    v <- r[[paste0("fdr_at_", TH[i])]][keep]; v <- v[is.finite(v)]
    if (length(v) < 20) return(NULL)
    data.frame(model=r$model[1], arm=r$arm[1], method=as.character(r$method[1]),
               t=TV[i], m=mean(v), se=sd(v)/sqrt(length(v)), ncell=length(v))}))))
fdr$grp <- ifelse(fdr$method %in% names(FOCUS), fdr$method, "other methods")
fdr$grp <- factor(fdr$grp, c(names(FOCUS), "other methods"))

bound <- data.frame(t = seq(.50, .99, length.out = 60))
bound$m <- 1 - bound$t
p5 <- ggplot(fdr, aes(t, m, group = method, colour = grp)) +
  geom_line(data = bound, aes(t, m), inherit.aes = FALSE,
            colour = REF, linetype = "22", linewidth = .45) +
  geom_line(data = subset(fdr, grp == "other methods"), linewidth = .32, alpha = .75) +
  geom_ribbon(data = subset(fdr, grp != "other methods"),
              aes(ymin = m-2*se, ymax = m+2*se, fill = grp),
              colour = NA, alpha = .18, show.legend = FALSE) +
  geom_line(data = subset(fdr, grp != "other methods"), linewidth = .62) +
  geom_point(data = subset(fdr, grp != "other methods"), size = 1.25) +
  facet_grid(model ~ arm) +
  scale_colour_manual(values = c(FOCUS, `other methods` = GREY),
                      limits = c(names(FOCUS), "other methods"), name = NULL) +
  scale_fill_manual(values = FOCUS, guide = "none") +
  scale_x_continuous(breaks = TV, labels = c(".50",".80",".90",".95",".99")) +
  labs(x = "PIP threshold t", y = "False discovery rate (mean over cells, ±2 SE)",
       title = "False discovery rate against PIP threshold, six strata",
       subtitle = paste("Red dashed line is y = 1 - t, the highest false discovery rate a",
                        "calibrated method can have at that threshold. Rates formed within a",
                        "cell then averaged over cells;\ncells selecting nothing are excluded,",
                        "not scored zero. marginal_z selects nothing anywhere and is absent.")) +
  th + theme(legend.position = "right")
ggsave(file.path(FIG,"fig_fdr.pdf"), p5, width = 9.2, height = 5.6)

cat("figs 4-5 written\n")
