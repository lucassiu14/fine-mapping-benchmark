#!/usr/bin/env Rscript
suppressWarnings(suppressMessages(library(ggplot2)))
options(stringsAsFactors=FALSE)
OUT <- "/private/tmp/claude-502/-Users-ls1020/cdc231a7-d353-43cf-ac37-fd155cce1deb/scratchpad/figs_ultra"
TAB <- "/private/tmp/claude-502/-Users-ls1020/cdc231a7-d353-43cf-ac37-fd155cce1deb/scratchpad/tabs_ultra"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
th <- theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),
      strip.background=element_rect(fill="grey92"), legend.position="bottom")
sv <- function(g,f,w=10,h=6) ggsave(file.path(OUT,f),g,width=w,height=h,dpi=130)
FOCUS <- c("polyfun_est","polyfun_ldsc","finemap","susie","funmap","sparsepro",
           "fb_xregion","beatrice","functional_beatrice","susie_inf","sbayesrc","abf")

## ---- F1: variance decomposition heatmaps -----------------------------------
vd <- read.csv(file.path(TAB,"variance_decomposition_long.csv"))
vd <- vd[vd$method %in% FOCUS & vd$metric %in% c("ap","hi_pip_reliab"),]
vd$metric <- factor(vd$metric, levels=c("ap","hi_pip_reliab"),
                    labels=c("macro-AUPRC","high-PIP reliability"))
vd$variable <- factor(vd$variable, levels=c("LD condition","number of causals S","heritability phi",
                      "region length","p_causal","effect model","annotation type","annotation enrichment"))
vd$method <- factor(vd$method, levels=rev(FOCUS))
g <- ggplot(vd, aes(variable, method, fill=eta2))+
  geom_tile(color="white",linewidth=.4)+
  geom_text(aes(label=ifelse(eta2>=1, sprintf("%.0f",eta2), ifelse(eta2>=0.1,"<1","~0"))),
            size=2.7, color=ifelse(vd$eta2>35,"white","grey15"))+
  facet_wrap(~metric)+th+
  scale_fill_gradient(low="#f7fbff",high="#08306b",name=expression(eta^2*" (%)"))+
  theme(axis.text.x=element_text(angle=35,hjust=1), legend.position="right")+
  labs(title="Which experimental variables actually move each method?",
       subtitle="eta-squared from a per-method main-effects ANOVA: % of that method's variance explained by each variable",
       x=NULL,y=NULL)
sv(g,"u1_variance_decomp.png", w=11, h=5.6)

## ---- F2: in-sample standing per arm (ranked lollipop; no label collisions)
ins <- do.call(rbind, lapply(c("none","binary","continuous"), function(a){
  z <- read.csv(file.path(TAB,paste0("insample_standing_",a,".csv"))); z$arm <- a; z }))
ins <- ins[ins$method %in% FOCUS,]
ins$arm <- factor(ins$arm, levels=c("none","binary","continuous"),
                  labels=c("no annotations","binary annotations","continuous annotations"))
# order methods by mean in-sample AP across arms
ord <- names(sort(tapply(ins$median_AP, ins$method, mean), decreasing=FALSE))
ins$method <- factor(ins$method, levels=ord)
g <- ggplot(ins, aes(median_AP, method))+
  geom_segment(aes(x=0, xend=median_AP, y=method, yend=method), color="grey80", linewidth=.5)+
  geom_point(aes(fill=pct_FDRok, size=pct_best), shape=21, color="grey25")+
  facet_wrap(~arm)+th+
  scale_fill_gradient2(low="#D55E00",mid="#F0E442",high="#009E73",midpoint=50,
                       name="% cells FDR-controlled", limits=c(0,100))+
  scale_size_continuous(range=c(2,8), name="% cells outright best")+
  theme(legend.position="right", legend.box="vertical")+
  labs(title="IN-SAMPLE LD: each method judged on its own terms, per annotation arm",
       subtitle="bar = median macro-AUPRC; point size = how often it is the best deployable method; colour = how often its FDR is controlled",
       x="median macro-AUPRC", y=NULL)
sv(g,"u2_insample_standing.png", w=12, h=5.6)

## ---- F3: robustness ---------------------------------------------------------
rob <- read.csv(file.path(TAB,"robustness.csv"), check.names=FALSE)
rob <- rob[rob$method %in% FOCUS,]
rl <- data.frame(method=rep(rob$method,3),
                 ld=factor(rep(c("in-sample","ref500","ref200"),each=nrow(rob)),
                           levels=c("in-sample","ref500","ref200")),
                 ap=c(rob$`AP_in-sample`,rob$AP_ref500,rob$AP_ref200))
ord <- rob$method[order(-rob$AP_retained)]
rl$method <- factor(rl$method, levels=ord)
g <- ggplot(rl, aes(ld, ap, group=method, color=method))+
  geom_line(linewidth=.9)+geom_point(size=2)+th+
  facet_wrap(~method, nrow=2)+theme(legend.position="none")+
  labs(title="Robustness to LD mis-specification (panels ordered by fraction of AUPRC retained)",
       subtitle="median macro-AUPRC across the three LD conditions, common y-axis",
       x=NULL,y="median macro-AUPRC")
sv(g,"u3_robustness.png", w=12, h=5.6)

## ---- F4: in-sample sensitivity profiles (S, phi, region) -------------------
mk <- function(f, v, lab){
  w <- read.csv(file.path(TAB,f), check.names=FALSE)
  w <- w[w$method %in% FOCUS,]
  long <- reshape(w, idvar="method", varying=setdiff(names(w),"method"),
                  v.names="ap", timevar="lev", times=setdiff(names(w),"method"), direction="long")
  long$lev <- as.numeric(long$lev); long$variable <- lab; long }
S <- rbind(mk("insample_by_S.csv","S","number of causals S"),
           mk("insample_by_phi.csv","phi","heritability phi"),
           mk("insample_by_region_size.csv","rs","region length"))
S$variable <- factor(S$variable, levels=c("number of causals S","heritability phi","region length"))
hl <- c("sbayesrc","abf","fb_xregion","polyfun_est","finemap")
S$grp <- ifelse(S$method %in% hl, S$method, "other")
g <- ggplot(S, aes(lev, ap, group=method, color=grp))+
  geom_line(data=subset(S,grp=="other"), linewidth=.5, alpha=.45)+
  geom_line(data=subset(S,grp!="other"), linewidth=1.1)+
  geom_point(data=subset(S,grp!="other"), size=1.7)+
  facet_wrap(~variable, scales="free_x")+th+scale_x_log10()+
  scale_color_manual(values=c(sbayesrc="#CC79A7", abf="#999999", fb_xregion="#D55E00",
                              polyfun_est="#009E73", finemap="#000000", other="grey65"), name=NULL)+
  labs(title="IN-SAMPLE LD: how each method responds to the three dominant difficulty axes",
       subtitle="median macro-AUPRC; grey = the remaining methods. Note abf's flat response to phi and sbayesrc's shallow decline in S.",
       x="level (log scale)", y="median macro-AUPRC")
sv(g,"u4_insample_sensitivity.png", w=12, h=5)

## ---- F5: calibration vs volume (the abstention trap) -----------------------
fd <- read.csv(file.path(TAB,"fdr_at_thresholds.csv"))
fd <- fd[fd$threshold==0.95 & fd$method %in% FOCUS & fd$ld %in% c("in-sample","ref500"),]
fd$ld <- factor(fd$ld, levels=c("in-sample","ref500"))
g <- ggplot(fd, aes(n_sel, emp_FDR, color=ld))+
  geom_hline(yintercept=0.05, linetype=2, color="grey40")+
  geom_point(size=3)+geom_text(aes(label=method), size=2.4, vjust=-0.8, show.legend=FALSE)+
  facet_wrap(~ld)+th+theme(legend.position="none")+scale_x_log10()+
  scale_color_manual(values=c(`in-sample`="#0072B2", ref500="#D55E00"))+
  labs(title="Calibration must be read together with volume: FDR at PIP >= 0.95 vs number of calls made",
       subtitle="dashed = the 0.05 promised. Bottom-right = trustworthy AND productive; bottom-left = trustworthy by abstaining. Note sbayesrc calls 8.7M variants at 98% FDR.",
       x="number of variants selected at PIP >= 0.95 (log scale)", y="empirical FDR")
sv(g,"u5_calibration_volume.png", w=11, h=5.2)

## ---- F6: interchangeability ------------------------------------------------
pr <- read.csv(file.path(TAB,"interchangeability.csv"))
pr <- pr[pr$m1 %in% FOCUS & pr$m2 %in% FOCUS,]
sym <- rbind(pr, data.frame(m1=pr$m2,m2=pr$m1,med_abs_diff=pr$med_abs_diff,pct_within_01=pr$pct_within_01))
g <- ggplot(sym, aes(m1,m2,fill=pct_within_01))+
  geom_tile(color="white",linewidth=.4)+
  geom_text(aes(label=sprintf("%.0f",pct_within_01)), size=2.6,
            color=ifelse(sym$pct_within_01>40,"white","grey20"))+
  th+theme(axis.text.x=element_text(angle=40,hjust=1), legend.position="right")+
  scale_fill_gradient(low="#fff7ec",high="#7f0000",name="% cells\nwithin 0.01 AP")+
  labs(title="Which methods are effectively interchangeable?",
       subtitle="% of scenario cells where two methods' macro-AUPRC differs by less than 0.01",
       x=NULL,y=NULL)
sv(g,"u6_interchangeability.png", w=9.5, h=7.5)

cat("wrote:\n"); print(list.files(OUT))
