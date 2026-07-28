#!/usr/bin/env Rscript
suppressWarnings(suppressMessages(library(ggplot2)))
options(stringsAsFactors=FALSE)
OUT <- "/private/tmp/claude-502/-Users-ls1020/cdc231a7-d353-43cf-ac37-fd155cce1deb/scratchpad/figs3"
TAB <- "/private/tmp/claude-502/-Users-ls1020/cdc231a7-d353-43cf-ac37-fd155cce1deb/scratchpad/tabs3"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE); dir.create(TAB, showWarnings=FALSE, recursive=TRUE)

sc <- readRDS("results/iter003/combined_scenario_metrics.rds")
ldlev <- c("in-sample","ref500","ref200")
sc$ld <- factor(ifelse(is.na(sc$n_ref),"in-sample",paste0("ref",sc$n_ref)), levels=ldlev)
CORE <- c("fb_xregion","fb_pooled","functional_beatrice","beatrice","finemap","susie","polyfun_est","polyfun_oracle")
LAB  <- c(fb_xregion="fb_xregion (joint LassoNet)", fb_pooled="fb_pooled (joint linear)",
          functional_beatrice="functional_beatrice (per-locus)", beatrice="beatrice (no annot)",
          finemap="finemap", susie="susie", polyfun_est="polyfun_est", polyfun_oracle="polyfun_oracle (ceiling)")
# Okabe-Ito; joint models get the two most salient colours
COL <- c(fb_xregion="#D55E00", fb_pooled="#E69F00", functional_beatrice="#CC79A7",
         beatrice="#0072B2", finemap="#000000", susie="#56B4E9",
         polyfun_est="#009E73", polyfun_oracle="#999999")
th <- theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),
      strip.background=element_rect(fill="grey92"), legend.position="bottom")
save2 <- function(g,f,w=10,h=5.6) ggsave(file.path(OUT,f),g,width=w,height=h,dpi=130)
mq <- function(df,val,by){
  a<-aggregate(df[[val]],df[by],function(x) c(m=median(x,na.rm=TRUE),
     lo=unname(quantile(x,.25,na.rm=TRUE)),hi=unname(quantile(x,.75,na.rm=TRUE))))
  a<-do.call(data.frame,a); names(a)[(ncol(a)-2):ncol(a)]<-c("m","lo","hi"); a}

d <- sc[sc$method %in% CORE & sc$annotation_type %in% c("binary","continuous"),]
d$method <- factor(d$method, levels=CORE)

## FIG 5 - REAL FDR curves (dumped for iteration 003)
FM <- c("fb_xregion","fb_pooled","functional_beatrice","beatrice","finemap","susie")
jm <- unique(sc[,c("job_dir","annotation_type","n_ref","model")])
# keep only the job_dirs and methods we plot BEFORE merging - the raw file is ~9M rows
keepjob <- jm$job_dir[jm$annotation_type=="binary" & (is.na(jm$n_ref) | jm$n_ref==500)]
fdr <- readRDS("results/iter003/combined_fdr_curves.rds")
fdr <- fdr[fdr$method %in% FM & fdr$job_dir %in% keepjob,]
cat("fdr rows after filter:", nrow(fdr), "\n")
fdr <- merge(fdr, jm, by="job_dir")
fdr$ld <- factor(ifelse(is.na(fdr$n_ref),"in-sample",paste0("ref",fdr$n_ref)), levels=ldlev)
f <- fdr[fdr$method %in% FM & fdr$annotation_type=="binary" & fdr$ld %in% c("in-sample","ref500"),]
agg <- aggregate(cbind(tp,fp,fn) ~ threshold + method + ld, f, sum)
agg$fdr <- agg$fp/pmax(agg$tp+agg$fp,1); agg$nsel <- agg$tp+agg$fp
agg <- agg[agg$nsel >= 20 & agg$threshold >= 0.05 & agg$threshold <= 0.95,]
agg$method <- factor(agg$method, levels=FM)
g <- ggplot(agg, aes(threshold, fdr, color=method)) +
  geom_line(aes(y=1-threshold), color="grey55", linetype=2, linewidth=.7)+
  geom_line(linewidth=1)+facet_wrap(~ld)+th+
  scale_color_manual(values=COL, labels=LAB, name=NULL)+guides(color=guide_legend(nrow=2))+
  labs(title="FDR against PIP threshold (true 201-threshold curves, binary annotations)",
       subtitle="dashed grey = the calibration bound 1-t; a curve above it violates FDR control. >=20 variants selected.",
       x="PIP selection threshold t", y="empirical FDR among PIP >= t")
save2(g,"f5_fdr_curves.png", w=10, h=5.2)

## FIG 6 - reliability curves from the full calibration counts
cal <- readRDS("results/iter003/combined_pip_calibration.rds")
cal <- cal[cal$method %in% FM & cal$job_dir %in% keepjob,]
cat("cal rows after filter:", nrow(cal), "\n")
cal <- merge(cal, jm, by="job_dir")
cal$ld <- factor(ifelse(is.na(cal$n_ref),"in-sample",paste0("ref",cal$n_ref)), levels=ldlev)
cn <- cal[cal$method %in% FM & cal$annotation_type=="binary" & cal$ld %in% c("in-sample","ref500"),]
rel <- aggregate(cbind(n,n_causal,sum_pip)~ld+method+bin, cn, sum)
rel$mean_pip <- rel$sum_pip/rel$n; rel$reliab <- rel$n_causal/rel$n
rel$method <- factor(rel$method, levels=FM)
g <- ggplot(rel, aes(mean_pip, reliab, color=method)) +
  geom_abline(slope=1,intercept=0,linetype=2,color="grey50")+
  geom_line(linewidth=.9)+geom_point(size=1.7)+facet_wrap(~ld)+coord_equal(xlim=c(0,1),ylim=c(0,1))+th+
  scale_color_manual(values=COL,labels=LAB,name=NULL)+guides(color=guide_legend(nrow=2))+
  labs(title="PIP reliability curves (binary annotations)",
       subtitle="P(causal | PIP bin) vs mean PIP; on the diagonal = calibrated, below = over-confident",
       x="mean PIP in bin", y="fraction truly causal")
save2(g,"f6_reliability.png", w=10, h=5.4)

cat("curve figs done\n")
