suppressPackageStartupMessages(library(ggplot2))
`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
FIG <- "results/iter004_full/report/figs"
L3 <- readRDS("results/iter004/combined_scenario_metrics.rds"); L3<-L3[!L3$method %in% EXCLUDE_METHODS,]
G  <- readRDS("results/iter004/iter004_sense_g/sense_g_ap.rds")
INK<-"#222222";MUTED<-"#666666"
th<-theme_minimal(base_size=8.5)+theme(panel.grid.minor=element_blank(),
  panel.grid.major=element_line(colour="#E8E8E8",linewidth=.3),
  axis.title=element_text(colour=INK),axis.text=element_text(colour=MUTED),
  strip.text=element_text(colour=INK,face="bold",size=8),
  plot.title=element_text(colour=INK,face="bold",size=10),
  plot.subtitle=element_text(colour=MUTED,size=7.5),legend.title=element_blank(),
  legend.text=element_text(colour=INK,size=7),legend.position="right",legend.key.height=unit(9,"pt"))
sv<-function(p,f,w=9.5,h=5) ggsave(file.path(FIG,f),p,width=w,height=h,device="pdf")

# ---- F14: robustness to misspecification (THE headline)
d <- L3[L3$annotation_type!="none",]
a <- aggregate(ap ~ method + model + annotation_type, d, function(z) mean(z,na.rm=TRUE))
w <- reshape(a, idvar=c("method","annotation_type"), timevar="model", direction="wide")
names(w) <- sub("^ap\\.","",names(w))
w$drop <- w$sparse - w$sparse_inf
w$ann <- factor(c(binary="Binary annotations",continuous="Continuous annotations")[w$annotation_type],
                levels=c("Binary annotations","Continuous annotations"))
KEY <- c("fb_xregion","functional_beatrice","fb_pooled","beatrice",
         "polyfun_est","polyfun_ldsc","polyfun_oracle","susie","susie_inf","finemap_inf","funmap")
w$hl <- w$method %in% KEY
ordm <- aggregate(drop~method,w,mean); w$method <- factor(w$method, levels=ordm$method[order(ordm$drop)])
p <- ggplot(w, aes(y=method)) +
  geom_segment(aes(x=sparse_inf, xend=sparse, yend=method), colour="#CCCCCC", linewidth=.9) +
  geom_point(aes(x=sparse), colour="#0072B2", size=1.9) +
  geom_point(aes(x=sparse_inf), colour="#D55E00", size=1.9) +
  facet_wrap(~ann) + th +
  labs(title="Robustness to model misspecification: AP under the sparse model (blue) versus sparse + infinitesimal (orange)",
       x="Mean AP", y=NULL,
       subtitle="Bar length is the cost of misspecification. Methods are ordered by that cost; short bars degrade gracefully.")
sv(p,"f14_misspecification.pdf",9.5,5)

# ---- F15: ablation conditional on S, both arms
dl <- G$delta
dl$arm <- factor(ifelse(dl$model=="sparse","Sparse","Sparse + infinitesimal"),
                 levels=c("Sparse","Sparse + infinitesimal"))
dl$ann <- factor(c(none="No annotations",binary="Binary",continuous="Continuous")[dl$annotation_type],
                 levels=c("No annotations","Binary","Continuous"))
dl$pair <- factor(dl$pair, levels=c("annot","joint"),
                  labels=c("annotations (FB - BEATRICE)","region sharing (xregion - FB)"))
mk <- function(xv,xl,fn) {
  a <- aggregate(delta ~ arm + ann + pair + get(xv), dl, function(z) c(m=mean(z),se=sd(z)/sqrt(length(z))))
  o <- data.frame(arm=a$arm,ann=a$ann,pair=a$pair,x=a[[4]],m=a$delta[,"m"],se=a$delta[,"se"])
  p <- ggplot(o,aes(x,m,colour=pair,group=pair)) +
    geom_hline(yintercept=0,colour="#AAAAAA",linewidth=.4) +
    geom_errorbar(aes(ymin=m-1.96*se,ymax=m+1.96*se),width=0,linewidth=.45) +
    geom_line(linewidth=.65)+geom_point(size=1.3)+
    scale_colour_manual(values=c("#D55E00","#0072B2"))+facet_grid(arm~ann)+th+
    labs(title=paste("Sense G: where the ablation gains come from -",xl),
         x=xl,y="paired difference in AP (95% CI)",
         subtitle="Unannotated panels are the negative control and sit exactly on zero.")
  sv(p,fn,9.5,4.6) }
mk("S","number of causal variants S","f15_ablation_by_S.pdf")
mk("phi","region heritability phi","f16_ablation_by_phi.pdf")
cat("figures 14-16 done\n")
