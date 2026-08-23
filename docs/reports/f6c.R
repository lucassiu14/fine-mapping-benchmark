suppressPackageStartupMessages(library(ggplot2))
`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
FIG <- "results/iter004_full/report/figs"
L2 <- readRDS("results/iter004/combined_replicate_metrics.rds"); L2<-L2[!L2$method %in% EXCLUDE_METHODS,]
L3 <- readRDS("results/iter004/combined_scenario_metrics.rds"); L3<-L3[!L3$method %in% EXCLUDE_METHODS,]
V <- readRDS("results/iter004/iter004_sense_v/sense_v.rds")
D <- readRDS("results/iter004/iter004_sense_d/sense_d_ap.rds")
G <- readRDS("results/iter004/iter004_sense_g/sense_g_ap.rds")
SN <- c(sparse_none="Sparse / none", sparse_binary="Sparse / binary", sparse_cont="Sparse / continuous",
        sparseinf_none="Sp+inf / none", sparseinf_binary="Sp+inf / binary", sparseinf_cont="Sp+inf / continuous")
LV <- unname(SN)
FAM <- c("polyfun_oracle"="#000000","fb_xregion"="#0072B2","fb_pooled"="#56B4E9",
         "functional_beatrice"="#D55E00","beatrice"="#E69F00","susie"="#009E73","polyfun_est"="#CC79A7")
OTHER<-"other methods (11)"; PAL<-c(FAM,setNames("#BBBBBB",OTHER)); LEV<-c(names(FAM),OTHER)
grp<-function(m) factor(ifelse(m %in% names(FAM),m,OTHER),levels=LEV)
INK<-"#222222";MUTED<-"#666666"
th<-theme_minimal(base_size=8.5)+theme(panel.grid.minor=element_blank(),
  panel.grid.major=element_line(colour="#E8E8E8",linewidth=.3),
  axis.title=element_text(colour=INK),axis.text=element_text(colour=MUTED),
  strip.text=element_text(colour=INK,face="bold",size=7.5),
  plot.title=element_text(colour=INK,face="bold",size=10),
  plot.subtitle=element_text(colour=MUTED,size=7.5),legend.title=element_blank(),
  legend.text=element_text(colour=INK,size=7),legend.position="right",legend.key.height=unit(9,"pt"))
sv<-function(p,f,w=9.5,h=4.6) ggsave(file.path(FIG,f),p,width=w,height=h,device="pdf")
SC<-scale_colour_manual(values=PAL,breaks=LEV,drop=FALSE)

# F10 variance budget
b<-V$budgets[V$budgets$response=="ap" & !is.nan(V$budgets$noise_floor),]
bl<-do.call(rbind,lapply(c("noise_floor","structure_mains","structure_2way","structure_3way_plus"),
  function(k) data.frame(stratum=b$stratum,method=b$method,part=k,val=100*b[[k]])))
bl$part<-factor(bl$part,levels=c("noise_floor","structure_mains","structure_2way","structure_3way_plus"),
  labels=c("noise floor","main effects","2-way","3-way+"))
bl$stratum<-factor(SN[bl$stratum],levels=LV)
om<-b[b$stratum=="sparseinf_cont",]; bl$method<-factor(bl$method,levels=om$method[order(om$structure_mains)])
p<-ggplot(bl,aes(val,method,fill=part))+geom_col(width=.7)+
  scale_fill_manual(values=c("noise floor"="#BBBBBB","main effects"="#0072B2","2-way"="#56B4E9","3-way+"="#E69F00"))+
  facet_wrap(~stratum,ncol=3)+th+
  labs(title="Sense V: the variance budget bounds every claim below it",
       x="% of cell-mean sum of squares",y=NULL,
       subtitle="Strata with a single design row have no identifiable sigma^2_u and are omitted.")
sv(p,"f10_variance_budget.pdf",9.5,7)

# F11 drivers by response
tm<-V$terms[V$terms$order==1 & V$terms$term %in% c("S","phi","region_size","enrich","pc"),]
tm<-tm[!grepl("_none$",tm$stratum),]
ag<-aggregate(omega2~term+response+stratum,tm,function(z) 100*mean(z,na.rm=TRUE))
ag$stratum<-factor(SN[ag$stratum],levels=LV)
ag$response<-factor(ag$response,levels=c("ap","fdr_at_90","total_mass_ratio"),
  labels=c("AP (ranking)","FDR at 0.9","total mass ratio"))
ag$term<-factor(ag$term,levels=c("S","phi","pc","region_size","enrich"))
p<-ggplot(ag,aes(omega2,term,fill=response))+geom_col(position=position_dodge(width=.75),width=.68)+
  scale_fill_manual(values=c("AP (ranking)"="#0072B2","FDR at 0.9"="#D55E00","total mass ratio"="#E69F00"))+
  facet_wrap(~stratum,ncol=2)+th+
  labs(title="Sense V: accuracy and calibration are driven by different factors",
       x=expression("mean noise-corrected "*omega^2*" (%)"),y=NULL,
       subtitle="Negative values are reported unclamped: they diagnose a response too noisy to decompose at this replication.")
sv(p,"f11_sense_v_responses.pdf",9.5,5)

# F12 Sense D winners
w<-D$winners; w$stratum<-factor(SN[w$stratum],levels=LV)
wc<-as.data.frame(table(w$winner,w$decided,w$stratum)); names(wc)<-c("method","decided","stratum","n")
wc<-wc[wc$n>0,]
wc$decided<-factor(ifelse(wc$decided=="TRUE","statistically decided","not decided"),
  levels=c("not decided","statistically decided"))
tot<-tapply(wc$n,wc$method,sum); wc$method<-factor(wc$method,levels=names(sort(tot)))
p<-ggplot(wc,aes(n,method,fill=decided))+geom_col(width=.7)+
  scale_fill_manual(values=c("not decided"="#DDDDDD","statistically decided"="#0072B2"))+
  facet_wrap(~stratum,ncol=3,scales="free_x")+th+
  labs(title="Sense D: which methods are ever best in a cell, and does it survive a test?",
       x="cells where the method scored highest",y=NULL,
       subtitle="Grey dominates everywhere: the nominal winner usually cannot be separated from the runner-up.")
sv(p,"f12_sense_d_winners.pdf",9.5,6.5)

# F13 ablation forest, six strata
key<-paste(L2$job_dir,L2$S,L2$phi,L2$region_size,L2$iter,sep="|")
res<-list(); k<-0
for (mo in c("sparse","sparse_inf")) for (a in c("none","binary","continuous")) {
  s<-L2$model==mo & L2$annotation_type==a; if(!any(s)) next
  sub<-function(m){v<-L2$ap[s&L2$method==m];names(v)<-key[s&L2$method==m];v}
  for (pp in list(c("functional_beatrice","beatrice","annotations"),
                  c("fb_xregion","functional_beatrice","region sharing"),
                  c("fb_xregion","beatrice","combined"))) {
    x<-sub(pp[1]);y<-sub(pp[2]);ix<-intersect(names(x),names(y))
    d<-x[ix]-y[ix]; d<-d[is.finite(d)]; if(!length(d)) next
    k<-k+1
    res[[k]]<-data.frame(arm=ifelse(mo=="sparse","Sparse","Sparse + infinitesimal"),
      ann=c(none="No annotations",binary="Binary",continuous="Continuous")[a],
      contrast=pp[3],mean=mean(d),lo=mean(d)-1.96*sd(d)/sqrt(length(d)),
      hi=mean(d)+1.96*sd(d)/sqrt(length(d))) } }
fr<-do.call(rbind,res)
fr$arm<-factor(fr$arm,levels=c("Sparse","Sparse + infinitesimal"))
fr$ann<-factor(fr$ann,levels=c("No annotations","Binary","Continuous"))
fr$contrast<-factor(fr$contrast,levels=c("combined","region sharing","annotations"))
p<-ggplot(fr,aes(mean,contrast))+geom_vline(xintercept=0,colour="#AAAAAA",linewidth=.4)+
  geom_errorbarh(aes(xmin=lo,xmax=hi),height=0,colour="#0072B2",linewidth=.7)+
  geom_point(colour="#0072B2",size=2)+facet_grid(arm~ann)+th+
  labs(title="The Functional BEATRICE ablation, paired within replicate, all six strata",
       x="difference in AP (95% CI)",y=NULL,
       subtitle="Both unannotated panels are negative controls and are exactly zero over 6,250 paired replicates.")
sv(p,"f13_ablation.pdf",9.5,4)
cat("figures 10-13 done\n")
