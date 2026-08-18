suppressPackageStartupMessages({library(ggplot2)})
`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
FIG <- "results/iter004_sparse/report/figs"
L1 <- readRDS("/tmp/l1_sparse/L1_sparsearm.rds"); L1 <- L1[!L1$method %in% EXCLUDE_METHODS,]
L3 <- readRDS("results/iter004_sparse/combined_scenario_metrics.rds"); L3 <- L3[!L3$method %in% EXCLUDE_METHODS,]
STRAT <- c(none="No annotations", binary="Binary annotations", continuous="Continuous annotations")
L3$stratum <- factor(STRAT[L3$annotation_type], levels=STRAT)

# Legend series in FIXED order, never cycled. The 12 unhighlighted methods
# collapse to one muted entry rather than each taking a hue - 18 distinct hues
# would be unreadable and is exactly what a categorical palette must not do.
FAM <- c("polyfun_oracle"="#000000", "fb_xregion"="#0072B2", "fb_pooled"="#56B4E9",
         "functional_beatrice"="#D55E00", "beatrice"="#E69F00", "susie"="#009E73")
OTHER <- "other methods (12)"
PAL <- c(FAM, setNames("#BBBBBB", OTHER))
LEV <- c(names(FAM), OTHER)
grp <- function(m) factor(ifelse(m %in% names(FAM), m, OTHER), levels=LEV)
INK <- "#222222"; MUTED <- "#666666"
th <- theme_minimal(base_size=9) + theme(
  panel.grid.minor=element_blank(), panel.grid.major=element_line(colour="#E8E8E8", linewidth=.3),
  axis.title=element_text(colour=INK), axis.text=element_text(colour=MUTED),
  strip.text=element_text(colour=INK, face="bold", size=8),
  plot.title=element_text(colour=INK, face="bold", size=10),
  plot.subtitle=element_text(colour=MUTED, size=8),
  legend.title=element_blank(), legend.text=element_text(colour=INK, size=7),
  legend.position="right", legend.key.height=unit(9,"pt"))
sv <- function(p,f,w=7,h=4.2) ggsave(file.path(FIG,f), p, width=w, height=h, device="pdf")
agg <- function(d, by, v, f=function(z) mean(z,na.rm=TRUE)) {
  k <- do.call(paste, c(d[by], sep="\r")); s <- tapply(d[[v]], k, f)
  o <- as.data.frame(do.call(rbind, strsplit(names(s),"\r",fixed=TRUE)), stringsAsFactors=FALSE)
  names(o) <- by; o[[v]] <- as.numeric(s); o }
SC <- scale_colour_manual(values=PAL, breaks=LEV, drop=FALSE)

# --------------------------------------------------- FIG 1
d  <- agg(L3, c("stratum","method"), "ap")
se <- agg(L3, c("stratum","method"), "ap", function(z) sd(z,na.rm=TRUE)/sqrt(sum(!is.na(z))))
d$se <- se$ap; d$stratum <- factor(d$stratum, levels=STRAT); d$g <- grp(d$method)
ord <- agg(L3[L3$annotation_type=="continuous",], "method", "ap")
d$method <- factor(d$method, levels=ord$method[order(ord$ap)])
p <- ggplot(d, aes(ap, method, colour=g)) +
  geom_errorbarh(aes(xmin=ap-1.96*se, xmax=ap+1.96*se), height=0, colour="#CCCCCC", linewidth=.4) +
  geom_point(size=1.9) + SC + facet_wrap(~stratum) + th +
  labs(title="Average precision by method", x="Mean AP (95% CI)", y=NULL,
       subtitle="Sparse arm, 1,125 cells per method per stratum.")
sv(p, "fig01_ap_ranking.pdf", 9, 4.6)

# --------------------------------------------------- FIG 2-4
mk <- function(xv, xlab, fname, sub) {
  dd <- agg(L3, c("stratum","method",xv), "ap"); dd[[xv]] <- as.numeric(dd[[xv]])
  dd$g <- grp(dd$method); dd$hl <- dd$method %in% names(FAM)
  p <- ggplot(dd, aes(.data[[xv]], ap, group=method, colour=g)) +
    geom_line(data=dd[!dd$hl,], linewidth=.3, alpha=.55) +
    geom_line(data=dd[dd$hl,], linewidth=.8) +
    geom_point(data=dd[dd$hl,], size=1) +
    SC + facet_wrap(~stratum) + th +
    labs(title=paste("Average precision against", xlab), x=xlab, y="Mean AP", subtitle=sub)
  sv(p, fname, 9, 3.8)
}
mk("S","number of causal variants S","fig02_ap_vs_S.pdf",
   "Each line is one method, averaged over all other design factors.")
mk("phi","region heritability phi","fig03_ap_vs_phi.pdf",
   "Each line is one method, averaged over all other design factors.")
mk("region_size","region size p (variants)","fig04_ap_vs_p.pdf",
   "Nearly flat: over 500-2000 variants, region size is not what makes this hard.")

# --------------------------------------------------- FIG 5
de <- L3[L3$annotation_type!="none",]; de$enrichment_fold <- as.numeric(de$enrichment_fold)
dd <- agg(de, c("stratum","method","enrichment_fold"), "ap")
dd$enrichment_fold <- as.numeric(dd$enrichment_fold); dd$g <- grp(dd$method)
dd$hl <- dd$method %in% names(FAM)
p <- ggplot(dd, aes(enrichment_fold, ap, group=method, colour=g)) +
  geom_line(data=dd[!dd$hl,], linewidth=.3, alpha=.55) +
  geom_line(data=dd[dd$hl,], linewidth=.8) + geom_point(data=dd[dd$hl,], size=1) +
  SC + facet_wrap(~stratum) + th +
  labs(title="Average precision against annotation enrichment", x="enrichment fold", y="Mean AP",
       subtitle="Near-flat lines are the finding: annotation PRESENCE matters, its STRENGTH barely does.")
sv(p, "fig05_ap_vs_enrich.pdf", 7.6, 3.6)

# --------------------------------------------------- FIG 6 (labels removed)
d <- merge(agg(L3, c("stratum","method"), "total_mass_ratio"),
           agg(L3, c("stratum","method"), "fdr_at_90"), by=c("stratum","method"))
d$g <- grp(d$method); d$stratum <- factor(d$stratum, levels=STRAT)
p <- ggplot(d, aes(total_mass_ratio, fdr_at_90, colour=g)) +
  geom_vline(xintercept=1, colour="#AAAAAA", linewidth=.4, linetype="22") +
  geom_point(size=2.1) + SC + scale_x_log10() + facet_wrap(~stratum) + th +
  labs(title="Calibration: how much mass a method claims, versus how often it is wrong",
       x="total mass ratio (sum of PIPs / number of causals), log scale",
       y="FDR at PIP >= 0.9",
       subtitle="Dashed line = perfect mass calibration. Bottom-left is good: honest mass AND low false discovery.")
sv(p, "fig06_calibration_tradeoff.pdf", 9, 4)

# --------------------------------------------------- FIG 7
bands <- c(lo="[0,0.1)", mid="[0.1,0.5)", hi="[0.5,0.9)", top="[0.9,1]")
rows <- lapply(names(bands), function(b) {
  n  <- tapply(L1[[paste0("n_band_",b)]],       L1$method, sum, na.rm=TRUE)
  cc <- tapply(L1[[paste0("c_band_",b)]],       L1$method, sum, na.rm=TRUE)
  sp <- tapply(L1[[paste0("sum_pip_band_",b)]], L1$method, sum, na.rm=TRUE)
  data.frame(method=names(n), observed=as.numeric(cc)/pmax(as.numeric(n),1),
             predicted=as.numeric(sp)/pmax(as.numeric(n),1), n=as.numeric(n)) })
rb <- do.call(rbind, rows); rb <- rb[rb$n>0,]; rb$g <- grp(rb$method)
rb$hl <- rb$method %in% names(FAM)
p <- ggplot(rb, aes(predicted, observed, group=method, colour=g)) +
  geom_abline(slope=1, intercept=0, colour="#AAAAAA", linewidth=.4, linetype="22") +
  geom_line(data=rb[!rb$hl,], linewidth=.3, alpha=.5) +
  geom_line(data=rb[rb$hl,], linewidth=.8) + geom_point(data=rb[rb$hl,], size=1.2) +
  SC + coord_equal(xlim=c(0,1), ylim=c(0,1)) + th +
  labs(title="Reliability diagram: claimed PIP versus realised causal rate",
       x="mean PIP within band (what the method claims)",
       y="fraction truly causal (what happened)",
       subtitle="Four PIP bands pooled over the sparse arm. On the diagonal = calibrated; below = overconfident.")
sv(p, "fig07_reliability.pdf", 6.6, 4.8)

# --------------------------------------------------- FIG 8
ths <- c(50,80,90,95,99)
fd <- do.call(rbind, lapply(ths, function(t) {
  a <- agg(L3, c("stratum","method"), paste0("fdr_at_",t)); names(a)[3] <- "fdr"; a$thr <- t/100; a }))
fd$g <- grp(fd$method); fd$hl <- fd$method %in% names(FAM)
fd$stratum <- factor(fd$stratum, levels=STRAT)
p <- ggplot(fd, aes(thr, fdr, group=method, colour=g)) +
  geom_line(data=fd[!fd$hl,], linewidth=.3, alpha=.5) +
  geom_line(data=fd[fd$hl,], linewidth=.8) + geom_point(data=fd[fd$hl,], size=1) +
  SC + facet_wrap(~stratum) + th +
  labs(title="False discovery rate across PIP thresholds", x="PIP threshold t",
       y="FDR among variants with PIP >= t",
       subtitle="A calibrated method's FDR at threshold t should be at most 1-t (e.g. <=0.10 at t=0.9).")
sv(p, "fig08_fdr_thresholds.pdf", 9, 3.8)

# --------------------------------------------------- FIG 9 (labels removed)
d <- merge(agg(L3, c("stratum","method"), "set_size_95"),
           agg(L3, c("stratum","method"), "set_hit_95"), by=c("stratum","method"))
d$g <- grp(d$method); d$stratum <- factor(d$stratum, levels=STRAT)
p <- ggplot(d, aes(set_size_95, set_hit_95, colour=g)) +
  geom_hline(yintercept=.95, colour="#AAAAAA", linewidth=.4, linetype="22") +
  geom_point(size=2.1) + SC + scale_x_log10() + facet_wrap(~stratum) + th +
  labs(title="95% credible sets: follow-up burden versus how often the causal is inside",
       x="mean set size (variants), log scale", y="coverage (fraction of fits containing a causal)",
       subtitle="Top-left is good: small sets that still contain the truth. Dashed line = nominal 95%.")
sv(p, "fig09_credible_sets.pdf", 9, 4)
cat("figs 1-9 regenerated with legends\n")
