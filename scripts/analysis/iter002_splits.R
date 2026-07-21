#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter002_splits.R
#
# Iteration 002 report by data split -> results/iteration-002-splits.pdf figures.
#
# Splits the 16,875 scenarios into LD (in-sample / ref500 / ref200) x model
# (sparse / sparse_inf) x enrichment (2.7 / 5.4 / 8.1 / 10.8 / none) = 30 splits.
# NB: splitting by enrichment (not annotation type) means binary + continuous
# are pooled WITHIN a split; annotation type is instead shown as one of the
# "other variables" in the per-split graphs, so it stays visible.
#
# Per split:
#   * one large all-method table: method (row, ordered worst->best overall) x S
#     (col), faceted by phi. Colour = median AP as a fraction of the best
#     PRACTICAL method in that scenario; overlaid ! / x flag the GUARDED FDR
#     violation (>=20 selected) at 0.25 / 0.50.
#   * three metric graphs (AUPRC, guarded FDR violation, high-PIP reliability),
#     each faceted by EVERY other varying variable in the split
#     (S, phi, region length, and p_causal / annotation type where they vary).
#
# Uses the GUARDED metrics from the fixed collect (max_fdr_violation_n20), and
# leads calibration with high-PIP reliability, not ECE (ECE is dominated by the
# ~94% of observations with PIP < 0.1 - see the metric audit).
#
# Usage: Rscript scripts/analysis/iter002_splits.R [metrics_rds] [fig_dir]
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(grid); library(scales)})
args <- commandArgs(trailingOnly=TRUE)
MET  <- if (length(args)>=1) args[1] else "results/iter002_fixed/combined_scenario_metrics.rds"
OUT  <- if (length(args)>=2) args[2] else "results/iter002_fixed/figures_splits"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

s <- readRDS(MET)
num <- function(x) as.numeric(as.character(x))
s$S <- num(s$S); s$phi <- num(s$phi); s$rs <- num(s$region_size)
s$ld  <- factor(ifelse(is.na(s$n_ref),"in-sample",paste0("ref",s$n_ref)),
                levels=c("in-sample","ref500","ref200"))
s$enr <- ifelse(is.na(s$enrichment_fold),"none",as.character(s$enrichment_fold))
s$ann <- s$annotation_type
s$pc  <- ifelse(is.na(s$p_causal),"sparse",as.character(s$p_causal))
s$scen<- paste(s$job_dir,s$S,s$phi,s$region_size,sep="|")
PRACT <- setdiff(unique(s$method),"polyfun_oracle")

# ap relative to best PRACTICAL method in each scenario (for the table colour)
bp <- tapply(s$ap[s$method %in% PRACT], s$scen[s$method %in% PRACT], max, na.rm=TRUE)
s$best_ap <- as.numeric(bp[s$scen]); s$ap_rel <- ifelse(s$best_ap>0, s$ap/s$best_ap, NA_real_)
# global method order (worst -> best by overall relative AP; best plotted on top)
MORD <- names(sort(tapply(s$ap_rel, s$method, median, na.rm=TRUE)))

INK<-"#0b0b0b"; INK2<-"#52514e"; SURF<-"#fcfcfb"; GRIDC<-"#e6e6e3"
RAMP <- c("#f7f7f5","#d9ecd9","#a8d5a8","#6bbb6b","#2f9e2f","#006300")
# 13-method palette (distinct hues)
PAL <- c(abf="#2a78d6", beatrice="#e34948", finemap="#006300",
         functional_beatrice="#eb6834", funmap="#9467bd", marginal_z="#8c8c8c",
         paintor="#c9a227", polyfun_est="#1baf7a", polyfun_ldsc="#17becf",
         sbayesrc="#d6277a", sparsepro="#4a3aa7", susie="#111111", susie_inf="#7fb800")
flag <- function(v) ifelse(is.na(v),"", ifelse(v>0.5,"×", ifelse(v>0.25,"!","")))
th <- theme_minimal(base_size=8)+theme(
  plot.background=element_rect(fill=SURF,colour=NA), panel.background=element_rect(fill=SURF,colour=NA),
  panel.grid.minor=element_blank(), panel.grid.major=element_line(colour=GRIDC,linewidth=.3),
  axis.title=element_text(colour=INK2,size=7.5), axis.text=element_text(colour=INK2,size=6.8),
  plot.title=element_text(colour=INK,face="bold",size=9.5), plot.subtitle=element_text(colour=INK2,size=7.2),
  strip.text=element_text(colour=INK,face="bold",size=7.3), legend.position="bottom",
  legend.title=element_blank(), legend.text=element_text(colour=INK2,size=6.6),
  legend.key.height=unit(7,"pt"), legend.margin=margin(t=-3), panel.spacing=unit(6,"pt"))
th_tab <- th + theme(panel.grid=element_blank())

sid <- function(ld,md,e) sprintf("%s__%s__enr-%s", ld, md, gsub("\\.","_",e))
VARLAB <- c(S="number of causal variants S", phi="per-variant heritability φ",
            rs="region length (variants)", pc="p_causal (sparse_inf mix)",
            ann="annotation type")

splits <- expand.grid(ld=levels(s$ld), model=c("sparse","sparse_inf"),
                      enr=c("none","2.7","5.4","8.1","10.8"), stringsAsFactors=FALSE)
meta <- list()

for (k in seq_len(nrow(splits))) {
  ld <- splits$ld[k]; md <- splits$model[k]; e <- splits$enr[k]
  d <- s[s$ld==ld & s$model==md & s$enr==e,]
  if (!nrow(d)) next
  id <- sid(ld,md,e)
  dp <- d[d$method %in% PRACT,]
  # which "other variables" vary here
  varset <- c("S","phi","rs")
  if (length(unique(d$pc))>1)  varset <- c(varset,"pc")
  if (length(unique(d$ann))>1) varset <- c(varset,"ann")
  nscen <- length(unique(d$scen))
  aggnote <- paste(setdiff(c("region length","annotation type","p_causal"),
                           c())[c(TRUE, "ann"%in%varset, "pc"%in%varset)], collapse=", ")
  meta[[id]] <- list(ld=ld, model=md, enr=e, nscen=nscen, varset=varset)

  # ---- TABLE: method x S faceted by phi (relative AP + guarded FDR flag) ----
  a <- aggregate(cbind(ap_rel, max_fdr_violation_n20) ~ method+phi+S, d, median)
  a$flag <- flag(a$max_fdr_violation_n20); a$method <- factor(a$method, levels=MORD)
  tab <- ggplot(a, aes(factor(S), method, fill=ap_rel))+
    geom_tile(colour=SURF, linewidth=.8)+
    geom_text(aes(label=flag), size=2.5, colour="#b00000", fontface="bold")+
    facet_wrap(~phi, nrow=1, labeller=labeller(phi=function(x) paste0("φ = ",x)))+
    scale_fill_gradientn(colours=RAMP, limits=c(0.3,1.0), oob=squish,
      name="median AP as a fraction of the best practical method")+
    labs(title=sprintf("%s / %s / enrichment %s : accuracy across heritability and causal count",
                       ld, md, e),
         subtitle="A run of dark cells across a row = the method holds up across every S at that heritability. Symbols = GUARDED FDR violation (>=20 selected): ! > 0.25, × > 0.5.\nEach square aggregates over region length and any other varying axis in this split (annotation type / p_causal where present).",
         x="number of causal variants S", y=NULL)+th_tab
  ggsave(file.path(OUT,sprintf("T_%s.png",id)), tab, width=9.4, height=4.4, dpi=190)

  # ---- long df for the metric graphs: median metric per (variable, level, method)
  buildlong <- function(metric){
    do.call(rbind, lapply(varset, function(v){
      ag <- aggregate(dp[[metric]], by=list(method=dp$method, lvl=dp[[v]]),
                      FUN=function(x) median(x,na.rm=TRUE))
      names(ag)[3] <- "y"
      ag$panel <- factor(VARLAB[[v]], levels=unname(VARLAB[varset]))
      # globally-unique, per-panel-ordered x factor: "<v>:<value>"
      lv <- ag$lvl
      ord <- if (v=="ann") sort(unique(as.character(lv))) else as.character(sort(unique(num(lv))))
      ag$xf <- factor(paste0(v,":",as.character(lv)), levels=paste0(v,":",ord))
      ag[,c("panel","xf","method","y","lvl")]
    }))
  }
  stripx <- function(x) sub("^[^:]+:","",x)
  mkgraph <- function(metric, title, ylab, ylim=NULL, hline=NULL){
    L <- buildlong(metric); L$method <- factor(L$method, levels=names(PAL))
    g <- ggplot(L, aes(xf, y, colour=method, group=method))
    if (!is.null(hline)) g <- g + geom_hline(yintercept=hline, linetype="22", colour=INK2, linewidth=.3)
    g <- g + geom_line(linewidth=.5, na.rm=TRUE)+geom_point(size=1.1, na.rm=TRUE)+
      facet_wrap(~panel, nrow=1, scales="free_x")+
      scale_x_discrete(labels=stripx)+
      scale_colour_manual(values=PAL, drop=FALSE)+
      guides(colour=guide_legend(nrow=2, override.aes=list(linewidth=1)))+
      labs(title=sprintf("%s  —  %s / %s / enrichment %s", title, ld, md, e),
           subtitle="One line per practical method (oracle excluded); each panel varies one axis, others aggregated over.",
           x=NULL, y=ylab)+th
    if (!is.null(ylim)) g <- g + coord_cartesian(ylim=ylim)
    g
  }
  wds <- max(7.6, 1.8 + 1.7*length(varset))
  ggsave(file.path(OUT,sprintf("gAUPRC_%s.png",id)),
         mkgraph("ap","AUPRC (average precision)","median AP", ylim=c(0,1)),
         width=wds, height=3.5, dpi=190)
  ggsave(file.path(OUT,sprintf("gFDR_%s.png",id)),
         mkgraph("max_fdr_violation_n20","FDR control (guarded max violation, lower is better)",
                 "median max FDR violation", ylim=c(0,1), hline=0.05),
         width=wds, height=3.5, dpi=190)
  ggsave(file.path(OUT,sprintf("gCAL_%s.png",id)),
         mkgraph("hi_pip_reliab","PIP calibration (high-PIP reliability, higher is better)",
                 "median reliability of PIP ≥ 0.9 calls", ylim=c(0,1), hline=0.9),
         width=wds, height=3.5, dpi=190)
  cat(sprintf("[%2d/30] %-34s n=%d  vars=%s\n", k, id, nscen, paste(varset,collapse=",")))
}
saveRDS(meta, file.path(OUT,"_meta.rds"))
cat("done:", length(meta), "splits\n")

# ---------------------------------------------------------------------------
# Per-split data-backed observations (consumed by the PDF assembler).
# ---------------------------------------------------------------------------
rows <- list()
for (id in names(meta)) {
  m <- meta[[id]]
  d <- s[s$ld==m$ld & s$model==m$model & s$enr==m$enr, ]
  dp <- d[d$method %in% PRACT, ]
  ap  <- sort(tapply(dp$ap, dp$method, median, na.rm=TRUE), decreasing=TRUE)
  vio <- tapply(dp$max_fdr_violation_n20, dp$method, median, na.rm=TRUE)
  rel <- tapply(dp$hi_pip_reliab, dp$method, median, na.rm=TRUE)
  safe <- names(vio)[!is.na(vio) & vio <= 0.05]
  safe_ap <- ap[names(ap) %in% safe]
  best_safe <- if (length(safe_ap)) names(safe_ap)[1] else NA_character_
  holds <- c()
  for (mm in names(ap)[1:min(4,length(ap))]) for (ph in sort(unique(dp$phi))) {
    sub <- dp[dp$method==mm & dp$phi==ph, ]
    byS <- tapply(sub$ap_rel, sub$S, median, na.rm=TRUE)
    if (length(byS)==5 && all(byS >= 0.9, na.rm=TRUE)) holds <- c(holds, sprintf("%s@phi=%s", mm, ph))
  }
  worst <- names(sort(vio, decreasing=TRUE))[1]
  rows[[length(rows)+1]] <- data.frame(id=id,
    best_ap=names(ap)[1], best_ap_v=round(ap[1],3),
    second=names(ap)[2], second_v=round(ap[2],3),
    best_safe=best_safe, best_safe_v=if (!is.na(best_safe)) round(ap[best_safe],3) else NA_real_,
    worst_fdr=worst, worst_fdr_v=round(vio[worst],3),
    best_rel=names(sort(rel, decreasing=TRUE))[1], best_rel_v=round(max(rel, na.rm=TRUE),2),
    holds=paste(head(holds,4), collapse="; "), stringsAsFactors=FALSE)
}
write.csv(do.call(rbind, rows), file.path(OUT, "split_observations.csv"), row.names=FALSE)
cat("wrote split_observations.csv\n")
