suppressPackageStartupMessages({library(ggplot2)})
`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
FIG <- "results/iter004_sparse/report/figs"
L1 <- readRDS("/tmp/l1_sparse/L1_sparsearm.rds")
L2 <- readRDS("results/iter004_sparse/combined_replicate_metrics.rds")
L3 <- readRDS("results/iter004_sparse/combined_scenario_metrics.rds")
for (d in list(L2,L3)) NULL
L2 <- L2[!L2$method %in% EXCLUDE_METHODS, ]; L3 <- L3[!L3$method %in% EXCLUDE_METHODS, ]
L1 <- L1[!L1$method %in% EXCLUDE_METHODS, ]

STRAT <- c(none="No annotations", binary="Binary annotations", continuous="Continuous annotations")
L3$stratum <- factor(STRAT[L3$annotation_type], levels = STRAT)
L2$stratum <- factor(STRAT[L2$annotation_type], levels = STRAT)

# --- palette: Okabe-Ito, colourblind-safe. FB family carries hue; the other
# --- 13 methods are context and share one muted ink so they never compete.
FAM <- c(beatrice="#E69F00", functional_beatrice="#D55E00",
         fb_pooled="#56B4E9", fb_xregion="#0072B2",
         polyfun_oracle="#000000", susie="#009E73")
COL <- function(m) ifelse(m %in% names(FAM), FAM[m], "#BBBBBB")
INK <- "#222222"; MUTED <- "#666666"
th <- theme_minimal(base_size = 9) + theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour = "#E8E8E8", linewidth = .3),
  axis.title = element_text(colour = INK), axis.text = element_text(colour = MUTED),
  strip.text = element_text(colour = INK, face = "bold", size = 8),
  plot.title = element_text(colour = INK, face = "bold", size = 10),
  plot.subtitle = element_text(colour = MUTED, size = 8),
  legend.title = element_blank(), legend.text = element_text(colour = INK, size = 7))
sv <- function(p, f, w = 7, h = 4.2) ggsave(file.path(FIG, f), p, width = w, height = h, device = "pdf")

agg <- function(d, by, v, f = function(z) mean(z, na.rm = TRUE)) {
  k <- do.call(paste, c(d[by], sep = "\r"))
  s <- tapply(d[[v]], k, f)
  parts <- do.call(rbind, strsplit(names(s), "\r", fixed = TRUE))
  out <- as.data.frame(parts, stringsAsFactors = FALSE); names(out) <- by
  out[[v]] <- as.numeric(s); out
}

# ---------------------------------------------------------------- FIG 1: AP ranking
d <- agg(L3, c("stratum","method"), "ap")
se <- agg(L3, c("stratum","method"), "ap", function(z) sd(z,na.rm=TRUE)/sqrt(sum(!is.na(z))))
d$se <- se$ap; d$stratum <- factor(d$stratum, levels = STRAT)
ord <- agg(L3[L3$annotation_type=="continuous",], "method", "ap")
d$method <- factor(d$method, levels = ord$method[order(ord$ap)])
d$col <- COL(as.character(d$method))
p <- ggplot(d, aes(ap, method)) +
  geom_errorbarh(aes(xmin=ap-1.96*se, xmax=ap+1.96*se), height=0, colour="#CCCCCC", linewidth=.4) +
  geom_point(aes(colour=col), size=1.9) + scale_colour_identity() +
  facet_wrap(~stratum) + th +
  labs(title="Average precision by method", x="Mean AP (95% CI)", y=NULL,
       subtitle="Sparse arm, 1,125 cells per method per stratum. Coloured: BEATRICE family, PolyFun-oracle (ceiling), SuSiE.")
sv(p, "fig01_ap_ranking.pdf", 8, 4.6)

# ---------------------------------------------------------------- FIG 2-4: AP vs design factors
mk <- function(xv, xlab, fname) {
  dd <- agg(L3, c("stratum","method",xv), "ap"); dd[[xv]] <- as.numeric(dd[[xv]])
  dd$col <- COL(dd$method); dd$hl <- dd$method %in% names(FAM)
  p <- ggplot(dd, aes(.data[[xv]], ap, group=method, colour=col)) +
    geom_line(data=dd[!dd$hl,], linewidth=.3, alpha=.55) +
    geom_line(data=dd[dd$hl,], linewidth=.8) +
    geom_point(data=dd[dd$hl,], size=1) +
    scale_colour_identity() + facet_wrap(~stratum) + th +
    labs(title=paste("Average precision against", xlab), x=xlab, y="Mean AP",
         subtitle="Grey: the other 13 methods. Every line is a mean over all other design factors.")
  sv(p, fname, 8, 3.8)
}
mk("S","number of causal variants S","fig02_ap_vs_S.pdf")
mk("phi","region heritability phi","fig03_ap_vs_phi.pdf")
mk("region_size","region size p (variants)","fig04_ap_vs_p.pdf")

# ---------------------------------------------------------------- FIG 5: enrichment
de <- L3[L3$annotation_type!="none",]; de$enrichment_fold <- as.numeric(de$enrichment_fold)
dd <- agg(de, c("stratum","method","enrichment_fold"), "ap")
dd$enrichment_fold <- as.numeric(dd$enrichment_fold); dd$col <- COL(dd$method)
dd$hl <- dd$method %in% names(FAM)
p <- ggplot(dd, aes(enrichment_fold, ap, group=method, colour=col)) +
  geom_line(data=dd[!dd$hl,], linewidth=.3, alpha=.55) +
  geom_line(data=dd[dd$hl,], linewidth=.8) + geom_point(data=dd[dd$hl,], size=1) +
  scale_colour_identity() + facet_wrap(~stratum) + th +
  labs(title="Average precision against annotation enrichment", x="enrichment fold", y="Mean AP",
       subtitle="Near-flat lines are the finding: the PRESENCE of annotations matters, their STRENGTH barely does.")
sv(p, "fig05_ap_vs_enrich.pdf", 6.5, 3.6)
cat("figs 1-5 done\n")
