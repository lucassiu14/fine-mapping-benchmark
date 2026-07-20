#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter002_per_method.R
#
# Iteration 002 per-method analysis -> results/iteration-002-per-method.pdf
#
# Grades every method on the THREE PRIMARY METRICS (not rankings, which are
# unreliable here: polyfun_est/polyfun_ldsc/susie are numerically IDENTICAL on
# the unannotated arms, so "wins" there are decided by tie-breaking alone):
#   AP (average precision - the correct AUPRC estimator), PIP calibration
#   (ECE / mass ratio / hi-PIP reliability), and FDR control (max violation).
#
# Produces, per method: a phi x S grid and a region-length x LD grid, faceted by
# the six model x annotation settings; plus large all-method tables per setting.
# Colour = median AP as a fraction of the best PRACTICAL method in that scenario;
# overlaid symbols flag untrustworthy probabilities (! > 0.25, x > 0.5 violation).
#
# NEVER aggregates across model (sparse/sparse_inf) or annotation type
# (none/binary/continuous). NA-as-a-level: n_ref=NA is IN-SAMPLE LD,
# enrichment_fold=NA is the none arm, p_causal=NA is the sparse model - all
# recoded explicitly because aggregate()/split() silently drop NA.
#
# Usage: Rscript scripts/analysis/iter002_per_method.R [metrics_rds] [fig_dir]
# =============================================================================
args <- commandArgs(trailingOnly=TRUE)
MET <- if (length(args)>=1) args[1] else "results/iter002/combined_scenario_metrics.rds"
FIGD<- if (length(args)>=2) args[2] else "results/iter002/figures_per_method"
dir.create(FIGD, showWarnings=FALSE, recursive=TRUE)
s <- readRDS(MET)
num <- function(x) as.numeric(as.character(x))
s$S <- num(s$S); s$phi <- num(s$phi); s$rs <- num(s$region_size)
s$ld  <- factor(ifelse(is.na(s$n_ref),"in-sample",paste0("ref",s$n_ref)),
                levels=c("in-sample","ref500","ref200"))
s$enr <- ifelse(is.na(s$enrichment_fold),"none-arm",as.character(s$enrichment_fold))
s$cell<- paste(s$model,s$annotation_type,sep="/")
s$scen<- paste(s$job_dir,s$S,s$phi,s$region_size,sep="|")
PRACT <- setdiff(unique(s$method),"polyfun_oracle")
# best AP among PRACTICAL methods, per scenario
bp <- tapply(s$ap[s$method %in% PRACT], s$scen[s$method %in% PRACT], max, na.rm=TRUE)
s$best_ap <- as.numeric(bp[s$scen])
s$ap_rel  <- ifelse(s$best_ap>0, s$ap/s$best_ap, NA_real_)


cat("=== PER-METHOD x CELL SUMMARY (medians; win% among 13 practical) ===\n")
r <- s[s$method %in% PRACT,]
o <- order(r$scen,-r$ap,na.last=TRUE); r <- r[o,]
g <- rle(r$scen)$lengths; r$rank <- unlist(lapply(g,seq_len)); r$is_win <- r$rank==1L
s$is_win <- FALSE; s$is_win[match(paste(r$scen,r$method)[r$is_win], paste(s$scen,s$method))] <- TRUE

CELLS <- c("sparse/none","sparse/binary","sparse/continuous",
           "sparse_inf/none","sparse_inf/binary","sparse_inf/continuous")
for (m in sort(unique(s$method))) {
  d <- s[s$method==m,]
  cat(sprintf("\n---- %s ----\n", m))
  out <- do.call(rbind, lapply(CELLS, function(cl){
    x <- d[d$cell==cl,]; if(!nrow(x)) return(NULL)
    data.frame(cell=cl, n=nrow(x),
      AP=round(median(x$ap,na.rm=TRUE),3), APrel=round(median(x$ap_rel,na.rm=TRUE),2),
      FDRviol=round(median(x$max_fdr_violation,na.rm=TRUE),3),
      ECE=round(median(x$ece,na.rm=TRUE),3),
      mass=round(median(x$total_mass_ratio,na.rm=TRUE),2),
      hiPIP=round(median(x$hi_pip_reliab,na.rm=TRUE),2),
      win=round(100*mean(x$is_win),1),
      broken=round(100*mean(x$max_fdr_violation>0.5,na.rm=TRUE)))
  }))
  print(out, row.names=FALSE)
}

# ---- per-method grids + large tables ----
suppressPackageStartupMessages({library(ggplot2); library(grid)})

OUT <- FIGD
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
CELLS <- c("sparse/none","sparse/binary","sparse/continuous",
           "sparse_inf/none","sparse_inf/binary","sparse_inf/continuous")
s$cell <- factor(s$cell, levels=CELLS)
INK<-"#0b0b0b"; INK2<-"#52514e"; SURF<-"#fcfcfb"
th <- theme_minimal(base_size=8)+theme(
  plot.background=element_rect(fill=SURF,colour=NA), panel.background=element_rect(fill=SURF,colour=NA),
  panel.grid=element_blank(), axis.title=element_text(colour=INK2,size=7.5),
  axis.text=element_text(colour=INK2,size=7), plot.title=element_text(colour=INK,face="bold",size=9.5),
  plot.subtitle=element_text(colour=INK2,size=7.3), strip.text=element_text(colour=INK,face="bold",size=7.3),
  legend.position="bottom", legend.title=element_text(colour=INK2,size=7),
  legend.text=element_text(colour=INK2,size=6.8), legend.key.height=unit(7,"pt"),
  legend.key.width=unit(20,"pt"), legend.margin=margin(t=-3), panel.spacing=unit(6,"pt"))
RAMP <- c("#f7f7f5","#d9ecd9","#a8d5a8","#6bbb6b","#2f9e2f","#006300")
flagsym <- function(v) ifelse(is.na(v),"", ifelse(v>0.5,"×", ifelse(v>0.25,"!","")))

# ---------- A) per-method grid: phi (y) x S (x), faceted by cell ----------
agg <- aggregate(cbind(ap_rel,max_fdr_violation) ~ method+cell+phi+S, s, median)
agg$flag <- flagsym(agg$max_fdr_violation)
for (m in sort(unique(agg$method))) {
  d <- agg[agg$method==m,]
  p <- ggplot(d, aes(factor(S), factor(phi), fill=ap_rel))+
    geom_tile(colour=SURF, linewidth=.9)+
    geom_text(aes(label=flag), size=2.9, colour="#b00000", fontface="bold")+
    facet_wrap(~cell, nrow=2)+
    scale_fill_gradientn(colours=RAMP, limits=c(0.3,1.0), oob=scales::squish,
      name="accuracy: median AP as a fraction of the best practical method in that scenario")+
    labs(title=sprintf("%s - accuracy across heritability x causal count, per simulation setting", m),
         subtitle="Darker green = closer to the best method. Symbols mark unreliable probabilities: ! = median max FDR violation > 0.25, × = > 0.5.\nEach square aggregates over region length (5), LD regime (3), enrichment (4 where applicable) and p_causal (4 for sparse_inf).",
         x="number of causal variants S", y="per-variant heritability phi")+th
  ggsave(file.path(OUT,sprintf("m_%s.png",m)), p, width=8.6, height=4.5, dpi=200)
}
# ---------- B) per-method grid: region length (x) x LD (y) ----------
agg2 <- aggregate(cbind(ap_rel,max_fdr_violation) ~ method+cell+ld+rs, s, median)
agg2$flag <- flagsym(agg2$max_fdr_violation)
for (m in sort(unique(agg2$method))) {
  d <- agg2[agg2$method==m,]
  p <- ggplot(d, aes(factor(rs), ld, fill=ap_rel))+
    geom_tile(colour=SURF, linewidth=.9)+
    geom_text(aes(label=flag), size=2.9, colour="#b00000", fontface="bold")+
    facet_wrap(~cell, nrow=2)+
    scale_fill_gradientn(colours=RAMP, limits=c(0.3,1.0), oob=scales::squish,
      name="accuracy: median AP as a fraction of the best practical method in that scenario")+
    labs(title=sprintf("%s - accuracy across region length x LD quality, per simulation setting", m),
         subtitle="Darker green = closer to the best method. ! = median max FDR violation > 0.25, × = > 0.5.\nEach square aggregates over S (5), heritability phi (5), enrichment (4 where applicable) and p_causal (4 for sparse_inf).",
         x="region length (variants)", y="LD regime")+th
  ggsave(file.path(OUT,sprintf("r_%s.png",m)), p, width=8.6, height=4.0, dpi=200)
}
cat("per-method grids:", length(list.files(OUT)), "files\n")

OUT <- FIGD
CELLS <- c("sparse/none","sparse/binary","sparse/continuous",
           "sparse_inf/none","sparse_inf/binary","sparse_inf/continuous")
INK<-"#0b0b0b"; INK2<-"#52514e"; SURF<-"#fcfcfb"
RAMP <- c("#f7f7f5","#d9ecd9","#a8d5a8","#6bbb6b","#2f9e2f","#006300")
flagsym <- function(v) ifelse(is.na(v),"", ifelse(v>0.5,"×", ifelse(v>0.25,"!","")))
th <- theme_minimal(base_size=8)+theme(
  plot.background=element_rect(fill=SURF,colour=NA), panel.background=element_rect(fill=SURF,colour=NA),
  panel.grid=element_blank(), axis.title=element_text(colour=INK2,size=7.5),
  axis.text=element_text(colour=INK2,size=7), axis.text.y=element_text(colour=INK,size=7.2),
  plot.title=element_text(colour=INK,face="bold",size=9.5), plot.subtitle=element_text(colour=INK2,size=7.2),
  strip.text=element_text(colour=INK,face="bold",size=7.6), legend.position="bottom",
  legend.title=element_text(colour=INK2,size=7), legend.text=element_text(colour=INK2,size=6.8),
  legend.key.height=unit(7,"pt"), legend.key.width=unit(22,"pt"), legend.margin=margin(t=-3),
  panel.spacing=unit(5,"pt"))
# order methods by overall accuracy so the table reads top=best
ord <- aggregate(ap_rel ~ method, s, median); ord <- ord$method[order(ord$ap_rel)]
slug <- function(x) gsub("[^a-z_]","", gsub("/","_",x))

# ---- TABLE 1 per cell: method (y) x S (x), faceted by phi ----
a1 <- aggregate(cbind(ap_rel,max_fdr_violation) ~ method+cell+phi+S, s, median)
a1$flag <- flagsym(a1$max_fdr_violation); a1$method <- factor(a1$method, levels=ord)
for (cl in CELLS) {
  d <- a1[a1$cell==cl,]
  p <- ggplot(d, aes(factor(S), method, fill=ap_rel))+
    geom_tile(colour=SURF, linewidth=.8)+
    geom_text(aes(label=flag), size=2.6, colour="#b00000", fontface="bold")+
    facet_wrap(~phi, nrow=1, labeller=labeller(phi=function(x) paste0("phi = ",x)))+
    scale_fill_gradientn(colours=RAMP, limits=c(0.3,1.0), oob=scales::squish,
      name="median AP as a fraction of the best practical method")+
    labs(title=sprintf("ALL METHODS - %s : accuracy across heritability and causal count", cl),
         subtitle="A run of dark cells across a row = the method holds up across every S at that heritability. ! = median max FDR violation > 0.25, × = > 0.5 (probabilities not trustworthy).\nEach square aggregates over region length (5), LD regime (3), enrichment (4 where applicable) and p_causal (4 for sparse_inf).\npolyfun_oracle is a ceiling (it is given the true annotation weights), not an achievable method.",
         x="number of causal variants S", y=NULL)+th
  ggsave(file.path(OUT,sprintf("T1_%s.png", slug(cl))), p, width=9.4, height=4.5, dpi=200)
}
# ---- TABLE 2 per cell: method (y) x region length (x), faceted by LD ----
a2 <- aggregate(cbind(ap_rel,max_fdr_violation) ~ method+cell+ld+rs, s, median)
a2$flag <- flagsym(a2$max_fdr_violation); a2$method <- factor(a2$method, levels=ord)
for (cl in CELLS) {
  d <- a2[a2$cell==cl,]
  p <- ggplot(d, aes(factor(rs), method, fill=ap_rel))+
    geom_tile(colour=SURF, linewidth=.8)+
    geom_text(aes(label=flag), size=2.6, colour="#b00000", fontface="bold")+
    facet_wrap(~ld, nrow=1)+
    scale_fill_gradientn(colours=RAMP, limits=c(0.3,1.0), oob=scales::squish,
      name="median AP as a fraction of the best practical method")+
    labs(title=sprintf("ALL METHODS - %s : accuracy across region length and LD quality", cl),
         subtitle="A run of dark cells across a row = the method holds up across every region length in that LD regime. ! = median max FDR violation > 0.25, × = > 0.5.\nEach square aggregates over S (5), heritability phi (5), enrichment (4 where applicable) and p_causal (4 for sparse_inf).",
         x="region length (variants)", y=NULL)+th
  ggsave(file.path(OUT,sprintf("T2_%s.png", slug(cl))), p, width=9.4, height=4.5, dpi=200)
}
cat("tables written\n")

# ---- export summary tables consumed by the PDF builder ----

CELLS <- c("sparse/none","sparse/binary","sparse/continuous",
           "sparse_inf/none","sparse_inf/binary","sparse_inf/continuous")
# per-method x cell primary metrics
out <- do.call(rbind, lapply(sort(unique(s$method)), function(m){
  d <- s[s$method==m,]
  do.call(rbind, lapply(CELLS, function(cl){ x<-d[d$cell==cl,]; if(!nrow(x)) return(NULL)
    data.frame(method=m, cell=cl,
      AP=round(median(x$ap,na.rm=TRUE),3), APrel=round(median(x$ap_rel,na.rm=TRUE),2),
      FDR=round(median(x$max_fdr_violation,na.rm=TRUE),3),
      ECE=round(median(x$ece,na.rm=TRUE),3),
      mass=round(median(x$total_mass_ratio,na.rm=TRUE),2),
      hiPIP=round(median(x$hi_pip_reliab,na.rm=TRUE),2),
      broken=round(100*mean(x$max_fdr_violation>0.5,na.rm=TRUE)))})) }))
write.csv(out,"/tmp/method_summary.csv",row.names=FALSE)
# per-method x LD (the dominant axis)
out2 <- do.call(rbind, lapply(sort(unique(s$method)), function(m){
  d <- s[s$method==m,]
  do.call(rbind, lapply(levels(d$ld), function(l){ x<-d[d$ld==l,]
    data.frame(method=m, ld=l, AP=round(median(x$ap,na.rm=TRUE),3),
      APrel=round(median(x$ap_rel,na.rm=TRUE),2),
      FDR=round(median(x$max_fdr_violation,na.rm=TRUE),3),
      mass=round(median(x$total_mass_ratio,na.rm=TRUE),1),
      hiPIP=round(median(x$hi_pip_reliab,na.rm=TRUE),2))})) }))
write.csv(out2,"/tmp/method_ld.csv",row.names=FALSE)
# combined: top-tier AND controlled
tc <- aggregate(cbind(good=(ap_rel>=0.95 & max_fdr_violation<=0.25)) ~ method+cell, s, function(x) round(100*mean(x)))
q <- reshape(tc, idvar="method", timevar="cell", direction="wide"); names(q)<-sub("good.","",names(q))
q <- q[order(-rowMeans(q[,-1],na.rm=TRUE)),]
write.csv(q,"/tmp/combined_good.csv",row.names=FALSE)
cat("exported\n"); cat(nrow(out),"method-cell rows\n")
