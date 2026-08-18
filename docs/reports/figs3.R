suppressPackageStartupMessages({library(ggplot2)})
`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
FIG <- "results/iter004_sparse/report/figs"
L2 <- readRDS("results/iter004_sparse/combined_replicate_metrics.rds"); L2 <- L2[!L2$method %in% EXCLUDE_METHODS,]
L3 <- readRDS("results/iter004_sparse/combined_scenario_metrics.rds"); L3 <- L3[!L3$method %in% EXCLUDE_METHODS,]
STRAT <- c(none="No annotations", binary="Binary annotations", continuous="Continuous annotations")
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

# ------------------------------------------- FIG 11: Sense V variance shares
ln <- readLines("results/iter004_sparse/iter004_sense_v/sense_v.txt")
i0 <- grep("^TOP TERMS PER", ln)[1]; i1 <- grep("^INTERACTION INVOLVEMENT", ln)[1]
ln <- ln[i0:i1]
anyh <- grepl("^[a-z_]+[.][a-z_0-9]+[.][a-z_0-9]+$", ln); aph <- grepl("^[a-z_]+[.][a-z_0-9]+[.]ap$", ln)
cur <- NA_character_; rows <- list(); k <- 0L
pat <- "^  ([A-Za-z_:]+) +([0-9]+) +(yes|-) +(-?[0-9.]+)%"
for (i in seq_along(ln)) {
  if (anyh[i]) { cur <- if (aph[i]) ln[i] else NA_character_; next }
  if (is.na(cur)) next
  m <- regmatches(ln[i], regexec(pat, ln[i]))[[1]]
  if (length(m)==5) { k<-k+1; p <- strsplit(cur,".",fixed=TRUE)[[1]]
    rows[[k]] <- data.frame(stratum=p[1], method=p[2], term=m[2], omega=as.numeric(m[5])) }
}
v <- do.call(rbind, rows); v <- v[!grepl(":", v$term), ]
v <- v[v$term %in% c("S","phi","region_size","enrich"), ]
v$stratum <- factor(STRAT[sub("^sparse_?","",v$stratum)], levels=STRAT)
v$stratum[is.na(v$stratum)] <- STRAT[1]
v$col <- COL(v$method)
v$term <- factor(v$term, levels=c("S","phi","region_size","enrich"))
p <- ggplot(v, aes(omega, reorder(method, omega), fill=col)) +
  geom_col(width=.7) + scale_fill_identity() + facet_grid(. ~ term) + th +
  labs(title="Sense V: share of cell-mean variance in AP explained by each design factor",
       x=expression(omega^2~" (%), noise-corrected"), y=NULL,
       subtitle="Averaged over strata. S dominates everywhere; enrichment fold barely registers.")
sv(p, "fig11_sense_v.pdf", 8.4, 5)

# ------------------------------------------- FIG 12: Sense D decidability
dl <- readLines("results/iter004_sparse/iter004_sense_d/sense_d_ap.txt")
rx <- "^(sparse_[a-z]+) +([a-z_]+) +([0-9.]+)% +([0-9.NaN]+)% +([0-9.]+)% +([0-9.]+)%"
mm <- regmatches(dl, regexec(rx, dl)); mm <- mm[lengths(mm)==7]
sd_ <- do.call(rbind, lapply(mm, function(m) data.frame(
  stratum=m[2], factor=m[3], flip_raw=as.numeric(m[4]),
  flip_gated=suppressWarnings(as.numeric(m[5])), decided=as.numeric(m[7]))))
sd2 <- rbind(transform(sd_, kind="apparent (raw flip rate)", val=flip_raw),
             transform(sd_, kind="survives a paired test", val=ifelse(is.na(flip_gated),0,flip_gated)))
sd2$stratum <- factor(STRAT[sub("^sparse_?","",sd2$stratum)], levels=STRAT)
sd2$stratum[is.na(sd2$stratum)] <- STRAT[1]
p <- ggplot(sd2, aes(val, factor, fill=kind)) +
  geom_col(position="dodge", width=.7) +
  scale_fill_manual(values=c("apparent (raw flip rate)"="#BBBBBB","survives a paired test"="#0072B2")) +
  facet_wrap(~stratum) + th + theme(legend.position="bottom") +
  labs(title="Sense D: how much of the apparent 'which method wins' structure is noise",
       x="% of adjacent level-pairs where the winner changes", y=NULL,
       subtitle="The gap between the bars IS the noise. Grey far exceeding blue means the winner is effectively arbitrary.")
sv(p, "fig12_sense_d.pdf", 8.4, 3.6)

# ------------------------------------------- FIG 13: ablation forest plot
key <- paste(L2$job_dir, L2$S, L2$phi, L2$region_size, L2$iter, sep="|")
res <- list(); k <- 0L
for (st in c("none","binary","continuous")) {
  sel <- L2$annotation_type == st
  sub <- function(m) { v <- L2$ap[sel & L2$method==m]; names(v) <- key[sel & L2$method==m]; v }
  for (p3 in list(c("functional_beatrice","beatrice","annotations\n(FB - BEATRICE)"),
                  c("fb_xregion","functional_beatrice","region sharing\n(FB-xregion - FB)"),
                  c("fb_xregion","beatrice","combined\n(FB-xregion - BEATRICE)"))) {
    a <- sub(p3[1]); b <- sub(p3[2]); ix <- intersect(names(a), names(b))
    d <- a[ix]-b[ix]; d <- d[is.finite(d)]; if (!length(d)) next
    k <- k+1
    res[[k]] <- data.frame(stratum=STRAT[st], contrast=p3[3], mean=mean(d),
                           lo=mean(d)-1.96*sd(d)/sqrt(length(d)),
                           hi=mean(d)+1.96*sd(d)/sqrt(length(d)), n=length(d))
  }
}
fr <- do.call(rbind, res); fr$stratum <- factor(fr$stratum, levels=STRAT)
fr$contrast <- factor(fr$contrast, levels=rev(unique(fr$contrast)))
p <- ggplot(fr, aes(mean, contrast)) +
  geom_vline(xintercept=0, colour="#AAAAAA", linewidth=.4) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0, colour="#0072B2", linewidth=.7) +
  geom_point(colour="#0072B2", size=2.2) +
  facet_wrap(~stratum) + th +
  labs(title="The ablation, paired within replicate",
       x="difference in AP (95% CI)", y=NULL,
       subtitle="Each contrast compares two methods on the SAME simulated data. The 'No annotations' panel is the negative control: exactly zero.")
sv(p, "fig13_ablation.pdf", 8.4, 3.2)
cat("figs 11-13 done\n")
