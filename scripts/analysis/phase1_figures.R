#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/phase1_figures.R
#
# Figures for the Iteration 001 (Phase 1) findings report
# (results/iteration-001-findings.pdf).
#
# Every metric is computed WITHIN one of the four non-poolable cells
# (sparse/none, sparse/binary, sparse_inf/none, sparse_inf/binary) and counts
# are pooled before rates are computed. Pooling across cells inverts
# conclusions rather than blurring them -- see finding 8 in the report.
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(grid)})
# Usage: Rscript scripts/analysis/phase1_figures.R [results_dir] [out_dir]
args <- commandArgs(trailingOnly = TRUE)
RES  <- if (length(args) >= 1) args[1] else "results"
OUT  <- if (length(args) >= 2) args[2] else file.path(RES, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

BLUE<-"#2a78d6"; RED<-"#e34948"; GREEN<-"#008300"; VIOLET<-"#4a3aa7"; ORANGE<-"#eb6834"
INK<-"#0b0b0b"; INK2<-"#52514e"; GRID<-"#e6e6e3"; SURF<-"#fcfcfb"

th <- theme_minimal(base_size = 9) + theme(
  plot.background=element_rect(fill=SURF,colour=NA),
  panel.background=element_rect(fill=SURF,colour=NA),
  panel.grid.minor=element_blank(),
  panel.grid.major=element_line(colour=GRID,linewidth=.3),
  axis.title=element_text(colour=INK2,size=8),
  axis.text=element_text(colour=INK2,size=7.5),
  plot.title=element_text(colour=INK,face="bold",size=10),
  plot.subtitle=element_text(colour=INK2,size=8),
  strip.text=element_text(colour=INK,face="bold",size=8),
  legend.position="bottom", legend.title=element_blank(),
  legend.text=element_text(colour=INK2,size=7.5),
  legend.key.height=unit(9,"pt"), legend.margin=margin(t=-4),
  panel.spacing=unit(11,"pt"))

fdr <- readRDS(file.path(RES, "combined_fdr_curves.rds"))
cal <- readRDS(file.path(RES, "combined_pip_calibration.rds"))
tagit <- function(d){
  d$arm   <- ifelse(grepl("anNone", d$job_dir), "none", "binary")
  d$model <- ifelse(grepl("sparse_inf", d$job_dir), "sparse_inf", "sparse")
  d$cell  <- paste(d$model, d$arm, sep="/"); d }
fdr <- tagit(fdr); cal <- tagit(cal)

mets <- function(d){
  a <- aggregate(cbind(tp,fp,fn) ~ threshold, data=d, FUN=sum)
  a <- a[order(a$threshold),]; ns <- a$tp+a$fp
  prec <- ifelse(ns>0,a$tp/ns,1); rec <- ifelse((a$tp+a$fn)>0,a$tp/(a$tp+a$fn),0)
  o <- order(rec,-prec); f <- ifelse(ns>0,a$fp/ns,0)
  c(ap=sum(diff(c(0,rec[o]))*prec[o]), viol=max(pmax(0,f-(1-a$threshold)))) }

grid_cmp <- do.call(rbind, lapply(split(fdr, list(fdr$cell,fdr$method,fdr$phi), drop=TRUE),
  function(x){ m<-mets(x); data.frame(cell=x$cell[1],method=x$method[1],
    phi=as.numeric(x$phi[1]), ap=m["ap"], viol=m["viol"], row.names=NULL)}))

# ---- FIG 1: violation vs phi, by cell -------------------------------------
sel <- c("susie","susie_inf","beatrice","sbayesrc","paintor")
d1 <- grid_cmp[grid_cmp$method %in% sel,]
d1$cell <- factor(d1$cell, levels=c("sparse/none","sparse/binary",
                                    "sparse_inf/none","sparse_inf/binary"))
lab1 <- d1[d1$phi==0.4,]
p1 <- ggplot(d1, aes(phi, viol, colour=method)) +
  geom_hline(yintercept=.05, linetype="22", colour=INK2, linewidth=.3) +
  geom_line(linewidth=.7) + geom_point(size=1.5) +
  scale_x_log10(breaks=c(.0075,.05,.1,.2,.4), labels=c("0.0075","0.05","0.1","0.2","0.4"),
                expand=expansion(mult=c(.05,.05))) +
  scale_colour_manual(values=c(susie=BLUE,susie_inf=GREEN,beatrice=VIOLET,
                               sbayesrc=ORANGE,paintor=RED)) +
  facet_wrap(~cell, nrow=1) +
  labs(title="Calibration failure is concentrated at weak signal",
       subtitle="Max FDR violation vs the FDR ≤ 1−t bound. Dashed line = 0.05 tolerance.",
       x="per-variant heritability φ (log scale)", y="max FDR violation") + th
ggsave(file.path(OUT,"fig1_viol_phi.png"), p1, width=9.4, height=3.3, dpi=200)

# ---- FIG 2: the U-shape, sparse_inf/none ----------------------------------
sel2 <- c("susie","susie_inf","finemap","sparsepro")
d2 <- grid_cmp[grid_cmp$cell=="sparse_inf/none" & grid_cmp$method %in% sel2,]
lab2 <- d2[d2$phi==0.4,]
p2 <- ggplot(d2, aes(phi, viol, colour=method)) +
  geom_hline(yintercept=.05, linetype="22", colour=INK2, linewidth=.3) +
  geom_line(linewidth=.8) + geom_point(size=1.8) +
  scale_x_log10(breaks=c(.0075,.05,.1,.2,.4), labels=c("0.0075","0.05","0.1","0.2","0.4"),
                expand=expansion(mult=c(.05,.05))) +
  scale_colour_manual(values=c(susie=BLUE,susie_inf=GREEN,finemap=RED,sparsepro=VIOLET)) +
  labs(title="Three methods break down again at high heritability — susie_inf does not",
       subtitle="sparse_inf / none. Violation is U-shaped: worst at φ=0.0075, recovers, then rises at φ=0.4.",
       x="per-variant heritability φ (log scale)", y="max FDR violation") + th
ggsave(file.path(OUT,"fig2_ushape.png"), p2, width=5.6, height=3.5, dpi=200)

# ---- FIG 3: annotation ceiling vs learned gain ----------------------------
d3 <- grid_cmp[grid_cmp$cell=="sparse/binary" &
               grid_cmp$method %in% c("susie","polyfun_est","polyfun_oracle"),]
d3$method <- factor(d3$method, levels=c("polyfun_oracle","polyfun_est","susie"),
                    labels=c("oracle annotations","learned (S-LDSC, pooled)","no annotations"))
lab3 <- d3[d3$phi==0.4,]
p3 <- ggplot(d3, aes(phi, ap, colour=method)) +
  geom_line(linewidth=.8) + geom_point(size=1.8) +
  scale_x_log10(breaks=c(.0075,.05,.1,.2,.4), labels=c("0.0075","0.05","0.1","0.2","0.4"),
                expand=expansion(mult=c(.05,.05))) +
  scale_colour_manual(values=c(BLUE,GREEN,VIOLET)) +
  labs(title="Learned annotations capture almost none of the available signal",
       subtitle="sparse/binary. The gap between oracle and learned is the unrealised opportunity.",
       x="per-variant heritability φ (log scale)", y="average precision") + th
ggsave(file.path(OUT,"fig3_ceiling.png"), p3, width=6.4, height=3.5, dpi=200)

# ---- FIG 4: % of ceiling captured -----------------------------------------
# Index columns BY NAME. reshape() orders wide columns by first appearance, not
# by factor level, so positional naming silently swaps oracle and learned.
d3b <- grid_cmp[grid_cmp$cell=="sparse/binary" &
                grid_cmp$method %in% c("susie","polyfun_est","polyfun_oracle"),]
w <- reshape(d3b[,c("method","phi","ap")], idvar="phi", timevar="method", direction="wide")
w <- data.frame(phi     = w$phi,
                oracle  = w[["ap.polyfun_oracle"]],
                learned = w[["ap.polyfun_est"]],
                susie   = w[["ap.susie"]])
stopifnot(all(w$oracle >= w$learned), all(w$learned >= w$susie))
w$pct <- 100*(w$learned-w$susie)/(w$oracle-w$susie)
w <- w[order(w$phi),]
cat("fig4 pct captured:", paste(round(w$pct), collapse=", "), "\n")
p4 <- ggplot(w, aes(factor(phi), pct)) +
  geom_col(fill=BLUE, width=.62) +
  geom_text(aes(label=sprintf("%.0f%%", pct)), vjust=-.5, size=2.8, colour=INK) +
  scale_y_continuous(limits=c(0,32), expand=expansion(mult=c(0,.06))) +
  labs(title="Share of the annotation ceiling that pooled learning actually captures",
       subtitle="sparse/binary. Worst exactly where annotations matter most.",
       x="per-variant heritability φ", y="% of ceiling captured") + th +
  theme(legend.position="none")
ggsave(file.path(OUT,"fig4_pct.png"), p4, width=5.2, height=2.7, dpi=200)

# ---- FIG 5: FB - beatrice, diverging --------------------------------------
d5 <- grid_cmp[grid_cmp$cell=="sparse/binary" &
               grid_cmp$method %in% c("beatrice","functional_beatrice"),]
w5 <- reshape(d5[,c("method","phi","ap")], idvar="phi", timevar="method", direction="wide")
# By name, never positionally -- see the fig 4 note above.
w5 <- data.frame(phi      = w5$phi,
                 beatrice = w5[["ap.beatrice"]],
                 fb       = w5[["ap.functional_beatrice"]])
stopifnot(!anyNA(w5$beatrice), !anyNA(w5$fb))
w5$gain <- w5$fb - w5$beatrice
w5 <- w5[order(w5$phi),]
cat("fig5 gain:", paste(sprintf("%+.3f", w5$gain), collapse=", "), "\n")
p5 <- ggplot(w5, aes(factor(phi), gain, fill=gain>0)) +
  geom_hline(yintercept=0, colour=INK2, linewidth=.4) +
  geom_col(width=.62) +
  geom_text(aes(label=sprintf("%+.3f",gain), vjust=ifelse(gain>0,-.5,1.4)),
            size=2.7, colour=INK) +
  scale_fill_manual(values=c(`FALSE`=RED,`TRUE`=BLUE), guide="none") +
  scale_y_continuous(expand=expansion(mult=c(.18,.18))) +
  labs(title="Per-locus annotation learning: harmful at weak signal, helpful at strong",
       subtitle="functional BEATRICE minus BEATRICE — identical inference, annotations the only difference (sparse/binary).",
       x="per-variant heritability φ", y="Δ average precision") + th +
  theme(legend.position="none")
ggsave(file.path(OUT,"fig5_fb.png"), p5, width=6.0, height=2.9, dpi=200)

# ---- FIG 6: implicit causal-count assumption ------------------------------
# Split by cell as well as method x S -- pooling mass across sparse/sparse_inf or
# none/binary would average over structures that behave differently.
mr <- do.call(rbind, lapply(split(cal, list(cal$cell,cal$method,cal$S), drop=TRUE), function(d){
  b <- aggregate(cbind(n,n_causal,sum_pip) ~ bin, data=d, FUN=sum)
  data.frame(cell=d$cell[1], method=d$method[1], S=as.numeric(d$S[1]),
             sumpip=sum(b$sum_pip)/ (sum(b$n_causal)/as.numeric(d$S[1])), row.names=NULL)}))
sel6 <- c("abf","susie","susie_inf","sbayesrc","beatrice")
d6 <- mr[mr$method %in% sel6,]
d6$cell <- factor(d6$cell, levels=c("sparse/none","sparse/binary",
                                    "sparse_inf/none","sparse_inf/binary"))
p6 <- ggplot(d6, aes(S, sumpip, colour=method)) +
  facet_wrap(~cell, nrow=1) +
  geom_abline(slope=1, intercept=0, linetype="22", colour=INK2, linewidth=.35) +
  geom_line(linewidth=.8) + geom_point(size=1.8) +
  scale_x_continuous(breaks=c(1,2,3,5,10), expand=expansion(mult=c(.05,.05))) +
  scale_colour_manual(values=c(abf=BLUE,susie=VIOLET,susie_inf=GREEN,
                               sbayesrc=ORANGE,beatrice=RED)) +
  labs(title="Every method encodes a fixed prior belief about how many causals exist",
       subtitle="Total PIP mass per region vs the true causal count. Dashed diagonal = truth (mass = S); flat lines ignore it.",
       x="true number of causal variants S", y="total PIP mass emitted") + th
ggsave(file.path(OUT,"fig6_mass.png"), p6, width=9.4, height=3.4, dpi=200)

cat("figures written to", OUT, "\n"); print(list.files(OUT))

# ---- FIG 7: reliability at the sharp end (PIP >= 0.8), within cells ---------
# A variant given PIP >= 0.8 should be causal >= 80% of the time. Within cell,
# pooling bin COUNTS then computing the rate.
hi <- do.call(rbind, lapply(split(cal, list(cal$cell, cal$method), drop=TRUE), function(d){
  b <- aggregate(cbind(n,n_causal) ~ bin, data=d, FUN=sum); top <- b[b$bin>=9,]
  if (!sum(top$n)) return(NULL)
  data.frame(cell=d$cell[1], method=d$method[1], n_hi=sum(top$n),
             prec=sum(top$n_causal)/sum(top$n), row.names=NULL)}))
hi$cell <- factor(hi$cell, levels=c("sparse/none","sparse/binary",
                                    "sparse_inf/none","sparse_inf/binary"))
ord <- aggregate(prec ~ method, data=hi, FUN=mean)
hi$method <- factor(hi$method, levels=ord$method[order(ord$prec)])
hi$ok <- hi$prec >= 0.8
p7 <- ggplot(hi, aes(prec, method, colour=ok)) +
  geom_vline(xintercept=.8, linetype="22", colour=INK2, linewidth=.35) +
  geom_point(size=1.9) +
  scale_colour_manual(values=c(`FALSE`=RED,`TRUE`=BLUE), guide="none") +
  scale_x_continuous(limits=c(.4,1), breaks=c(.4,.6,.8,1)) +
  facet_wrap(~cell, nrow=1) +
  labs(title="Reliability at the sharp end: not every method means 0.8 when it says 0.8",
       subtitle="Share of PIP ≥ 0.8 variants that are truly causal. Dashed line = the guarantee; red = broken.",
       x="precision among PIP ≥ 0.8 calls", y=NULL) + th +
  theme(panel.grid.major.y=element_line(colour=GRID, linewidth=.3))
ggsave(file.path(OUT,"fig7_sharp.png"), p7, width=9.4, height=3.0, dpi=200)
cat("fig7 rows:", nrow(hi), "\n")
