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

## FIG 1 - mass ratio (the headline calibration metric)
a <- mq(d,"total_mass_ratio",c("annotation_type","ld","method"))
g <- ggplot(a, aes(ld,m,color=method,group=method)) +
  geom_hline(yintercept=1, linetype=2, color="grey50") +
  geom_line(linewidth=1)+geom_point(size=2.2)+facet_wrap(~annotation_type)+
  scale_y_log10()+th+scale_color_manual(values=COL, labels=LAB, name=NULL)+
  guides(color=guide_legend(nrow=3))+
  labs(title="Posterior mass calibration: the joint prior cuts BEATRICE's over-confidence",
       subtitle="median total mass ratio (sum PIP / #causal); 1.0 = calibrated. Log scale.",
       x="LD condition", y="mass ratio (log10)")
save2(g,"f1_mass.png")

## FIG 2 - AUPRC
a <- mq(d,"ap",c("annotation_type","ld","method"))
g <- ggplot(a, aes(ld,m,color=method,group=method)) +
  geom_ribbon(aes(ymin=lo,ymax=hi,fill=method),alpha=.06,color=NA,show.legend=FALSE)+
  geom_line(linewidth=1)+geom_point(size=2.2)+facet_wrap(~annotation_type)+th+
  scale_color_manual(values=COL,labels=LAB,name=NULL)+scale_fill_manual(values=COL)+
  guides(color=guide_legend(nrow=3))+
  labs(title="Ranking accuracy: the shared LassoNet head matches or beats per-locus FB",
       subtitle="median macro-AUPRC with inter-quartile band",
       x="LD condition", y="macro-AUPRC")
save2(g,"f2_auprc.png")

## FIG 3 - high-PIP reliability
a <- mq(d,"hi_pip_reliab",c("annotation_type","ld","method"))
g <- ggplot(a, aes(ld,m,color=method,group=method)) +
  geom_hline(yintercept=1, linetype=2, color="grey50")+
  geom_line(linewidth=1)+geom_point(size=2.2)+facet_wrap(~annotation_type)+th+ylim(0,1.03)+
  scale_color_manual(values=COL,labels=LAB,name=NULL)+guides(color=guide_legend(nrow=3))+
  labs(title="Trustworthiness of confident calls",
       subtitle="median fraction of PIP >= 0.9 calls that are truly causal (1.0 = ideal)",
       x="LD condition", y="high-PIP reliability")
save2(g,"f3_reliab.png")

## FIG 4 - paired win rates vs functional_beatrice
key <- c("job_dir","S","phi","region_size")
pw <- reshape(d[,c(key,"method","ap","total_mass_ratio","hi_pip_reliab","max_fdr_violation_n20")],
              idvar=key, timevar="method", direction="wide")
meta <- unique(sc[,c("job_dir","annotation_type","n_ref","enrichment_fold")])
pw <- merge(pw, meta, by="job_dir")
pw$ld <- factor(ifelse(is.na(pw$n_ref),"in-sample",paste0("ref",pw$n_ref)), levels=ldlev)
mk <- function(j){
  z <- aggregate(cbind(
      `AUPRC`            = pw[[paste0("ap.",j)]] > pw[["ap.functional_beatrice"]],
      `mass ratio`       = abs(pw[[paste0("total_mass_ratio.",j)]]-1) < abs(pw[["total_mass_ratio.functional_beatrice"]]-1),
      `high-PIP reliab`  = pw[[paste0("hi_pip_reliab.",j)]] > pw[["hi_pip_reliab.functional_beatrice"]],
      `FDR violation`    = pw[[paste0("max_fdr_violation_n20.",j)]] < pw[["max_fdr_violation_n20.functional_beatrice"]]
    ) ~ annotation_type + ld, pw, function(x) 100*mean(x,na.rm=TRUE))
  z <- reshape(z, varying=list(3:6), v.names="pct", timevar="metric",
               times=c("AUPRC","mass ratio","high-PIP reliab","FDR violation"), direction="long")
  z$model <- LAB[[j]]; z }
W <- rbind(mk("fb_xregion"), mk("fb_pooled"))
W$metric <- factor(W$metric, levels=c("AUPRC","mass ratio","high-PIP reliab","FDR violation"))
g <- ggplot(W, aes(ld, pct, fill=annotation_type)) +
  geom_hline(yintercept=50, linetype=2, color="grey40")+
  geom_col(position="dodge")+facet_grid(model~metric)+th+
  scale_fill_manual(values=c(binary="#0072B2", continuous="#D55E00"), name=NULL)+
  labs(title="Paired within-cell comparison against per-locus functional BEATRICE",
       subtitle="% of cells where the joint model is better, same cell. Above the dashed 50% line = improvement.",
       x="LD condition", y="% of cells improved")
save2(g,"f4_winrate.png", w=11, h=6)
write.csv(W[,c("model","metric","annotation_type","ld","pct")], file.path(TAB,"winrate_vs_fb.csv"), row.names=FALSE)

## FIG 7 - which head? fb_xregion vs fb_pooled, paired
z <- aggregate(cbind(
    `AUPRC`      = pw[["ap.fb_xregion"]] > pw[["ap.fb_pooled"]],
    `mass ratio` = abs(pw[["total_mass_ratio.fb_xregion"]]-1) < abs(pw[["total_mass_ratio.fb_pooled"]]-1),
    `FDR violation` = pw[["max_fdr_violation_n20.fb_xregion"]] < pw[["max_fdr_violation_n20.fb_pooled"]]
  ) ~ annotation_type + ld, pw, function(x) 100*mean(x,na.rm=TRUE))
z <- reshape(z, varying=list(3:5), v.names="pct", timevar="metric",
             times=c("AUPRC","mass ratio","FDR violation"), direction="long")
g <- ggplot(z, aes(ld,pct,fill=annotation_type))+
  geom_hline(yintercept=50,linetype=2,color="grey40")+
  geom_col(position="dodge")+facet_wrap(~metric)+th+
  scale_fill_manual(values=c(binary="#0072B2",continuous="#D55E00"),name=NULL)+
  labs(title="Which shared prior head? LassoNet (fb_xregion) vs linear (fb_pooled)",
       subtitle="% of cells where fb_xregion is better. Above 50% = the neural head wins.",
       x="LD condition", y="% of cells fb_xregion better")
save2(g,"f7_head.png", w=10, h=4.8)

## FIG 8 - by number of causals
a <- mq(d[d$ld=="ref500",],"total_mass_ratio",c("S","method")); a$S <- as.numeric(a$S)
g <- ggplot(a, aes(S,m,color=method,group=method))+
  geom_hline(yintercept=1,linetype=2,color="grey50")+
  geom_line(linewidth=1)+geom_point(size=2)+th+scale_y_log10()+
  scale_x_continuous(breaks=c(1,2,3,5,10))+
  scale_color_manual(values=COL,labels=LAB,name=NULL)+guides(color=guide_legend(nrow=3))+
  labs(title="Mass inflation is worst when the true number of causals is small",
       subtitle="median mass ratio by S (annotated arms, n_ref=500). Every method assumes more signal than exists at S=1.",
       x="number of causal variants S", y="mass ratio (log10)")
save2(g,"f8_byS.png", w=9.5, h=5.2)

cat("wrote:\n"); print(list.files(OUT))
