suppressPackageStartupMessages({library(ggplot2)})
`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
FIG <- "results/iter004_sparse/report/figs"
L1 <- readRDS("/tmp/l1_sparse/L1_sparsearm.rds"); L1 <- L1[!L1$method %in% EXCLUDE_METHODS,]
L2 <- readRDS("results/iter004_sparse/combined_replicate_metrics.rds"); L2 <- L2[!L2$method %in% EXCLUDE_METHODS,]
L3 <- readRDS("results/iter004_sparse/combined_scenario_metrics.rds"); L3 <- L3[!L3$method %in% EXCLUDE_METHODS,]
STRAT <- c(none="No annotations", binary="Binary annotations", continuous="Continuous annotations")
L3$stratum <- factor(STRAT[L3$annotation_type], levels=STRAT)
FAM <- c(beatrice="#E69F00", functional_beatrice="#D55E00", fb_pooled="#56B4E9",
         fb_xregion="#0072B2", polyfun_oracle="#000000", susie="#009E73")
COL <- function(m) ifelse(m %in% names(FAM), FAM[m], "#BBBBBB")
INK <- "#222222"; MUTED <- "#666666"
th <- theme_minimal(base_size=9) + theme(
  panel.grid.minor=element_blank(), panel.grid.major=element_line(colour="#E8E8E8", linewidth=.3),
  axis.title=element_text(colour=INK), axis.text=element_text(colour=MUTED),
  strip.text=element_text(colour=INK, face="bold", size=8),
  plot.title=element_text(colour=INK, face="bold", size=10),
  plot.subtitle=element_text(colour=MUTED, size=8),
  legend.title=element_blank(), legend.text=element_text(colour=INK, size=7))
sv <- function(p,f,w=7,h=4.2) ggsave(file.path(FIG,f), p, width=w, height=h, device="pdf")
agg <- function(d, by, v, f=function(z) mean(z,na.rm=TRUE)) {
  k <- do.call(paste, c(d[by], sep="\r")); s <- tapply(d[[v]], k, f)
  o <- as.data.frame(do.call(rbind, strsplit(names(s),"\r",fixed=TRUE)), stringsAsFactors=FALSE)
  names(o) <- by; o[[v]] <- as.numeric(s); o }

# ------------------------------------------- FIG 6: mass ratio vs FDR@90
d <- merge(agg(L3, c("stratum","method"), "total_mass_ratio"),
           agg(L3, c("stratum","method"), "fdr_at_90"), by=c("stratum","method"))
d$col <- COL(d$method); d$stratum <- factor(d$stratum, levels=STRAT)
d$lab <- ifelse(d$method %in% names(FAM), d$method, "")
p <- ggplot(d, aes(total_mass_ratio, fdr_at_90)) +
  geom_vline(xintercept=1, colour="#AAAAAA", linewidth=.4, linetype="22") +
  geom_point(aes(colour=col), size=2) +
  geom_text(aes(label=lab), hjust=-0.12, size=2.1, colour=INK) +
  scale_colour_identity() + scale_x_log10() + facet_wrap(~stratum) + th +
  labs(title="Calibration: how much mass a method claims, versus how often it is wrong",
       x="total mass ratio (sum of PIPs / number of causals), log scale",
       y="FDR at PIP >= 0.9",
       subtitle="Dashed line = perfect mass calibration. Bottom-left is good: honest mass AND low false discovery.")
sv(p, "fig06_calibration_tradeoff.pdf", 8.4, 4)

# ------------------------------------------- FIG 7: reliability diagram from band counts
bands <- data.frame(band=c("lo","mid","hi","top"),
                    lab=factor(c("[0,0.1)","[0.1,0.5)","[0.5,0.9)","[0.9,1]"),
                               levels=c("[0,0.1)","[0.1,0.5)","[0.5,0.9)","[0.9,1]")))
rows <- list()
for (i in seq_len(nrow(bands))) {
  b <- bands$band[i]
  n <- tapply(L1[[paste0("n_band_",b)]], L1$method, sum, na.rm=TRUE)
  cc <- tapply(L1[[paste0("c_band_",b)]], L1$method, sum, na.rm=TRUE)
  sp <- tapply(L1[[paste0("sum_pip_band_",b)]], L1$method, sum, na.rm=TRUE)
  rows[[i]] <- data.frame(method=names(n), band=bands$lab[i],
                          observed=as.numeric(cc)/pmax(as.numeric(n),1),
                          predicted=as.numeric(sp)/pmax(as.numeric(n),1),
                          n=as.numeric(n))
}
rb <- do.call(rbind, rows); rb <- rb[rb$n > 0, ]; rb$col <- COL(rb$method)
rb$hl <- rb$method %in% names(FAM)
p <- ggplot(rb, aes(predicted, observed, group=method, colour=col)) +
  geom_abline(slope=1, intercept=0, colour="#AAAAAA", linewidth=.4, linetype="22") +
  geom_line(data=rb[!rb$hl,], linewidth=.3, alpha=.5) +
  geom_line(data=rb[rb$hl,], linewidth=.8) + geom_point(data=rb[rb$hl,], size=1.2) +
  scale_colour_identity() + coord_equal(xlim=c(0,1), ylim=c(0,1)) + th +
  labs(title="Reliability diagram: claimed PIP versus realised causal rate",
       x="mean PIP within band (what the method claims)",
       y="fraction truly causal (what happened)",
       subtitle="Pooled over the whole sparse arm, four PIP bands. On the diagonal = calibrated; below = overconfident.")
sv(p, "fig07_reliability.pdf", 5.2, 5)

# ------------------------------------------- FIG 8: FDR across thresholds
ths <- c(50,80,90,95,99)
rows <- lapply(ths, function(t) {
  a <- agg(L3, c("stratum","method"), paste0("fdr_at_",t)); names(a)[3] <- "fdr"; a$thr <- t/100; a })
fd <- do.call(rbind, rows); fd$col <- COL(fd$method); fd$hl <- fd$method %in% names(FAM)
fd$stratum <- factor(fd$stratum, levels=STRAT)
p <- ggplot(fd, aes(thr, fdr, group=method, colour=col)) +
  geom_line(data=fd[!fd$hl,], linewidth=.3, alpha=.5) +
  geom_line(data=fd[fd$hl,], linewidth=.8) + geom_point(data=fd[fd$hl,], size=1) +
  scale_colour_identity() + facet_wrap(~stratum) + th +
  labs(title="False discovery rate across PIP thresholds",
       x="PIP threshold t", y="FDR among variants with PIP >= t",
       subtitle="A calibrated method's FDR at threshold t should be at most 1-t (e.g. <=0.10 at t=0.9).")
sv(p, "fig08_fdr_thresholds.pdf", 8.4, 3.8)

# ------------------------------------------- FIG 9: credible set size vs coverage
d <- merge(agg(L3, c("stratum","method"), "set_size_95"),
           agg(L3, c("stratum","method"), "set_hit_95"), by=c("stratum","method"))
d$col <- COL(d$method); d$lab <- ifelse(d$method %in% names(FAM), d$method, "")
d$stratum <- factor(d$stratum, levels=STRAT)
p <- ggplot(d, aes(set_size_95, set_hit_95)) +
  geom_hline(yintercept=.95, colour="#AAAAAA", linewidth=.4, linetype="22") +
  geom_point(aes(colour=col), size=2) +
  geom_text(aes(label=lab), hjust=-.12, size=2.1, colour=INK) +
  scale_colour_identity() + scale_x_log10() + facet_wrap(~stratum) + th +
  labs(title="95% credible sets: how many variants must be followed up, and how often the causal is inside",
       x="mean set size (variants), log scale", y="coverage (fraction of fits containing a causal)",
       subtitle="Top-left is good: small sets that still contain the truth. Dashed line = nominal 95%.")
sv(p, "fig09_credible_sets.pdf", 8.4, 4)

# ------------------------------------------- FIG 10: Murphy decomposition
mu <- rbind(
  transform(agg(L3, c("stratum","method"), "rel"), part="REL (miscalibration, lower better)", v=rel)[,c("stratum","method","part","v")],
  transform(agg(L3, c("stratum","method"), "res"), part="RES (resolution, higher better)",  v=res)[,c("stratum","method","part","v")])
mu$col <- COL(mu$method); mu$stratum <- factor(mu$stratum, levels=STRAT)
ordm <- agg(L3[L3$annotation_type=="continuous",], "method", "res")
mu$method <- factor(mu$method, levels=ordm$method[order(ordm$res)])
p <- ggplot(mu, aes(v, method, fill=col)) + geom_col(width=.65) + scale_fill_identity() +
  facet_grid(stratum ~ part, scales="free_x") + th +
  labs(title="Murphy decomposition of the Brier score", x="component", y=NULL,
       subtitle="RES exposes the predictor that games calibration: a method emitting S/p everywhere has REL = 0 but RES = 0 too.")
sv(p, "fig10_murphy.pdf", 8.4, 7)
cat("figs 6-10 done\n")
