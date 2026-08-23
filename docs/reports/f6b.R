suppressPackageStartupMessages(library(ggplot2))
`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
FIG <- "results/iter004_full/report/figs"
L1 <- readRDS("results/iter004/combined_fit_metrics.rds"); L1 <- L1[!L1$method %in% EXCLUDE_METHODS,]
L3 <- readRDS("results/iter004/combined_scenario_metrics.rds"); L3 <- L3[!L3$method %in% EXCLUDE_METHODS,]
L2 <- readRDS("results/iter004/combined_replicate_metrics.rds"); L2 <- L2[!L2$method %in% EXCLUDE_METHODS,]
lab <- function(d){ d$arm <- factor(ifelse(d$model=="sparse","Sparse","Sparse + infinitesimal"),
        levels=c("Sparse","Sparse + infinitesimal"))
  d$ann <- factor(c(none="No annotations",binary="Binary",continuous="Continuous")[d$annotation_type],
        levels=c("No annotations","Binary","Continuous")); d }
L1<-lab(L1); L3<-lab(L3); L2<-lab(L2)
FAM <- c("polyfun_oracle"="#000000","fb_xregion"="#0072B2","fb_pooled"="#56B4E9",
         "functional_beatrice"="#D55E00","beatrice"="#E69F00","susie"="#009E73","polyfun_est"="#CC79A7")
OTHER<-"other methods (11)"; PAL<-c(FAM,setNames("#BBBBBB",OTHER)); LEV<-c(names(FAM),OTHER)
grp <- function(m) factor(ifelse(m %in% names(FAM),m,OTHER),levels=LEV)
INK<-"#222222"; MUTED<-"#666666"
th <- theme_minimal(base_size=8.5)+theme(panel.grid.minor=element_blank(),
  panel.grid.major=element_line(colour="#E8E8E8",linewidth=.3),
  axis.title=element_text(colour=INK),axis.text=element_text(colour=MUTED),
  strip.text=element_text(colour=INK,face="bold",size=8),
  plot.title=element_text(colour=INK,face="bold",size=10),
  plot.subtitle=element_text(colour=MUTED,size=7.5),legend.title=element_blank(),
  legend.text=element_text(colour=INK,size=7),legend.position="right",legend.key.height=unit(9,"pt"))
sv <- function(p,f,w=9.5,h=5.4) ggsave(file.path(FIG,f),p,width=w,height=h,device="pdf")
SC <- scale_colour_manual(values=PAL,breaks=LEV,drop=FALSE); FG <- facet_grid(arm~ann)
agg <- function(d,by,v,f=function(z) mean(z,na.rm=TRUE)){k<-do.call(paste,c(d[by],sep="\r"))
  s<-tapply(d[[v]],k,f); o<-as.data.frame(do.call(rbind,strsplit(names(s),"\r",fixed=TRUE)),stringsAsFactors=FALSE)
  names(o)<-by; o[[v]]<-as.numeric(s); o}
relev <- function(d){ d$arm<-factor(d$arm,levels=c("Sparse","Sparse + infinitesimal"))
  d$ann<-factor(d$ann,levels=c("No annotations","Binary","Continuous")); d}

# F6 calibration tradeoff
d <- merge(agg(L3,c("arm","ann","method"),"total_mass_ratio"),
           agg(L3,c("arm","ann","method"),"fdr_at_90"),by=c("arm","ann","method"))
d$g<-grp(d$method); d<-relev(d)
p <- ggplot(d,aes(total_mass_ratio,fdr_at_90,colour=g))+
  geom_vline(xintercept=1,colour="#AAAAAA",linewidth=.4,linetype="22")+
  geom_point(size=1.9)+SC+scale_x_log10()+FG+th+
  labs(title="Calibration: mass claimed versus errors made",
       x="total mass ratio (sum of PIPs / number of causals), log scale",y="FDR at PIP >= 0.9",
       subtitle="Dashed line = perfect mass calibration. Note the whole field shifts right and up under the infinitesimal model.")
sv(p,"f06_calibration.pdf")

# F7 reliability
bands<-c("lo","mid","hi","top")
rows<-list(); i<-0
for (mo in unique(L1$arm)) for (b in bands) { i<-i+1
  s<-L1$arm==mo
  n<-tapply(L1[[paste0("n_band_",b)]][s],L1$method[s],sum,na.rm=TRUE)
  cc<-tapply(L1[[paste0("c_band_",b)]][s],L1$method[s],sum,na.rm=TRUE)
  sp<-tapply(L1[[paste0("sum_pip_band_",b)]][s],L1$method[s],sum,na.rm=TRUE)
  rows[[i]]<-data.frame(arm=mo,method=names(n),observed=as.numeric(cc)/pmax(as.numeric(n),1),
                        predicted=as.numeric(sp)/pmax(as.numeric(n),1),n=as.numeric(n)) }
rb<-do.call(rbind,rows); rb<-rb[rb$n>0,]; rb$g<-grp(rb$method); rb$hl<-rb$method %in% names(FAM)
rb$arm<-factor(rb$arm,levels=c("Sparse","Sparse + infinitesimal"))
p<-ggplot(rb,aes(predicted,observed,group=method,colour=g))+
  geom_abline(slope=1,intercept=0,colour="#AAAAAA",linewidth=.4,linetype="22")+
  geom_line(data=rb[!rb$hl,],linewidth=.28,alpha=.5)+geom_line(data=rb[rb$hl,],linewidth=.75)+
  geom_point(data=rb[rb$hl,],size=1.1)+SC+coord_equal(xlim=c(0,1),ylim=c(0,1))+facet_wrap(~arm)+th+
  labs(title="Reliability: claimed PIP versus realised causal rate",
       x="mean PIP within band (what the method claims)",y="fraction truly causal (what happened)",
       subtitle="Under the infinitesimal model the truth counts only SPARSE causals, so every method appears overconfident.")
sv(p,"f07_reliability.pdf",8,4.4)

# F8 FDR thresholds
fd<-do.call(rbind,lapply(c(50,80,90,95,99),function(t){
  a<-agg(L3,c("arm","ann","method"),paste0("fdr_at_",t)); names(a)[4]<-"fdr"; a$thr<-t/100; a}))
fd$g<-grp(fd$method); fd$hl<-fd$method %in% names(FAM); fd<-relev(fd)
p<-ggplot(fd,aes(thr,fdr,group=method,colour=g))+
  geom_line(data=fd[!fd$hl,],linewidth=.28,alpha=.5)+geom_line(data=fd[fd$hl,],linewidth=.75)+
  geom_point(data=fd[fd$hl,],size=.9)+SC+FG+th+
  labs(title="False discovery rate across PIP thresholds",x="PIP threshold t",
       y="FDR among variants with PIP >= t",
       subtitle="A calibrated method satisfies FDR <= 1-t. The infinitesimal arm breaks this for the SuSiE/PolyFun family.")
sv(p,"f08_fdr.pdf")

# F9 credible sets
d<-merge(agg(L3,c("arm","ann","method"),"set_size_95"),
         agg(L3,c("arm","ann","method"),"set_hit_95"),by=c("arm","ann","method"))
d$g<-grp(d$method); d<-relev(d)
p<-ggplot(d,aes(set_size_95,set_hit_95,colour=g))+
  geom_hline(yintercept=.95,colour="#AAAAAA",linewidth=.4,linetype="22")+
  geom_point(size=1.9)+SC+scale_x_log10()+FG+th+
  labs(title="95% credible sets: follow-up burden versus coverage",
       x="mean set size (variants), log scale",y="coverage",
       subtitle="Top-left is good. Dashed line = nominal 95%.")
sv(p,"f09_credible_sets.pdf")
cat("figures 6-9 done\n")
