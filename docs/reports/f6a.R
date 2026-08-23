suppressPackageStartupMessages(library(ggplot2))
`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
FIG <- "results/iter004_full/report/figs"
L1 <- readRDS("results/iter004/combined_fit_metrics.rds"); L1 <- L1[!L1$method %in% EXCLUDE_METHODS,]
L3 <- readRDS("results/iter004/combined_scenario_metrics.rds"); L3 <- L3[!L3$method %in% EXCLUDE_METHODS,]
lab <- function(d) {
  d$arm <- factor(ifelse(d$model=="sparse","Sparse","Sparse + infinitesimal"),
                  levels=c("Sparse","Sparse + infinitesimal"))
  d$ann <- factor(c(none="No annotations",binary="Binary",continuous="Continuous")[d$annotation_type],
                  levels=c("No annotations","Binary","Continuous")); d }
L1 <- lab(L1); L3 <- lab(L3)
FAM <- c("polyfun_oracle"="#000000","fb_xregion"="#0072B2","fb_pooled"="#56B4E9",
         "functional_beatrice"="#D55E00","beatrice"="#E69F00","susie"="#009E73",
         "polyfun_est"="#CC79A7")
OTHER <- "other methods (11)"; PAL <- c(FAM, setNames("#BBBBBB",OTHER)); LEV <- c(names(FAM),OTHER)
grp <- function(m) factor(ifelse(m %in% names(FAM), m, OTHER), levels=LEV)
INK<-"#222222"; MUTED<-"#666666"
th <- theme_minimal(base_size=8.5) + theme(
  panel.grid.minor=element_blank(), panel.grid.major=element_line(colour="#E8E8E8",linewidth=.3),
  axis.title=element_text(colour=INK), axis.text=element_text(colour=MUTED),
  strip.text=element_text(colour=INK,face="bold",size=8),
  plot.title=element_text(colour=INK,face="bold",size=10),
  plot.subtitle=element_text(colour=MUTED,size=7.5),
  legend.title=element_blank(), legend.text=element_text(colour=INK,size=7),
  legend.position="right", legend.key.height=unit(9,"pt"))
sv <- function(p,f,w=9,h=5.4) ggsave(file.path(FIG,f),p,width=w,height=h,device="pdf")
SC <- scale_colour_manual(values=PAL, breaks=LEV, drop=FALSE)
agg <- function(d,by,v,f=function(z) mean(z,na.rm=TRUE)) {
  k <- do.call(paste,c(d[by],sep="\r")); s <- tapply(d[[v]],k,f)
  o <- as.data.frame(do.call(rbind,strsplit(names(s),"\r",fixed=TRUE)),stringsAsFactors=FALSE)
  names(o) <- by; o[[v]] <- as.numeric(s); o }
FG <- facet_grid(arm ~ ann)

# F1 AP ranking
d <- agg(L3,c("arm","ann","method"),"ap")
se <- agg(L3,c("arm","ann","method"),"ap",function(z) sd(z,na.rm=TRUE)/sqrt(sum(!is.na(z))))
d$se <- se$ap; d$g <- grp(d$method)
o <- agg(L3[L3$model=="sparse_inf"&L3$annotation_type=="continuous",],"method","ap")
d$method <- factor(d$method, levels=o$method[order(o$ap)])
d$arm <- factor(d$arm, levels=levels(L3$arm)); d$ann <- factor(d$ann, levels=levels(L3$ann))
p <- ggplot(d,aes(ap,method,colour=g)) +
  geom_errorbarh(aes(xmin=ap-1.96*se,xmax=ap+1.96*se),height=0,colour="#CCCCCC",linewidth=.35) +
  geom_point(size=1.7) + SC + FG + th +
  labs(title="Average precision by method, all six strata",x="Mean AP (95% CI)",y=NULL,
       subtitle="Ordered by sparse_inf/continuous. PolyFun-oracle reads the simulated truth and is a ceiling, not a competitor.")
sv(p,"f01_ap_ranking.pdf",9.5,7)

# F2-F4 design factors
mk <- function(xv,xl,fn,sub) {
  dd <- agg(L3,c("arm","ann","method",xv),"ap"); dd[[xv]] <- as.numeric(dd[[xv]])
  dd$g <- grp(dd$method); dd$hl <- dd$method %in% names(FAM)
  dd$arm <- factor(dd$arm, levels=levels(L3$arm)); dd$ann <- factor(dd$ann, levels=levels(L3$ann))
  p <- ggplot(dd,aes(.data[[xv]],ap,group=method,colour=g)) +
    geom_line(data=dd[!dd$hl,],linewidth=.28,alpha=.5) +
    geom_line(data=dd[dd$hl,],linewidth=.75) + geom_point(data=dd[dd$hl,],size=.9) +
    SC + FG + th + labs(title=paste("Average precision against",xl),x=xl,y="Mean AP",subtitle=sub)
  sv(p,fn,9.5,5) }
mk("S","number of causal variants S","f02_ap_vs_S.pdf","Each line is one method, averaged over all other design factors.")
mk("phi","region heritability phi","f03_ap_vs_phi.pdf","Note the dip at phi=0.6 for some methods - it is real and method-specific, not a design artefact.")
mk("region_size","nominal region size p","f04_ap_vs_p.pdf","Nominal sizes; the 2000 class is under-filled (median 1590 actual variants).")

# F5 p_causal (sparse_inf only)
di <- L3[L3$model=="sparse_inf",]
dd <- agg(di,c("ann","method","p_causal"),"ap"); dd$p_causal <- as.numeric(dd$p_causal)
dd$g <- grp(dd$method); dd$hl <- dd$method %in% names(FAM)
dd$ann <- factor(dd$ann, levels=levels(L3$ann))
p <- ggplot(dd,aes(p_causal,ap,group=method,colour=g)) +
  geom_line(data=dd[!dd$hl,],linewidth=.28,alpha=.5) +
  geom_line(data=dd[dd$hl,],linewidth=.75) + geom_point(data=dd[dd$hl,],size=.9) +
  SC + facet_wrap(~ann) + th +
  labs(title="Average precision against p_causal (sparse + infinitesimal arm only)",
       x="p_causal = fraction of heritability in the SPARSE component", y="Mean AP",
       subtitle="Higher p_causal means less infinitesimal background, so the sparse signal the metric scores is easier to find.")
sv(p,"f05_ap_vs_pcausal.pdf",9,3.6)
cat("figures 1-5 done\n")
