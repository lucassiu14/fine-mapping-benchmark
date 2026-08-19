suppressPackageStartupMessages({library(ggplot2)})
FIG <- "results/iter004_sparse/report/figs"
V <- readRDS("results/iter004_sparse/iter004_sense_v/sense_v.rds")
D <- readRDS("results/iter004_sparse/iter004_sense_d/sense_d_ap.rds")
G <- readRDS("results/iter004_sparse/iter004_sense_g/sense_g_ap.rds")
SN <- c(sparse_none="No annotations", sparse_binary="Binary annotations", sparse_cont="Continuous annotations")
FAM <- c("polyfun_oracle"="#000000","fb_xregion"="#0072B2","fb_pooled"="#56B4E9",
         "functional_beatrice"="#D55E00","beatrice"="#E69F00","susie"="#009E73")
OTHER <- "other methods (12)"; PAL <- c(FAM, setNames("#BBBBBB",OTHER)); LEV <- c(names(FAM),OTHER)
grp <- function(m) factor(ifelse(m %in% names(FAM), m, OTHER), levels=LEV)
INK<-"#222222"; MUTED<-"#666666"
th <- theme_minimal(base_size=9) + theme(
  panel.grid.minor=element_blank(), panel.grid.major=element_line(colour="#E8E8E8",linewidth=.3),
  axis.title=element_text(colour=INK), axis.text=element_text(colour=MUTED),
  strip.text=element_text(colour=INK,face="bold",size=8),
  plot.title=element_text(colour=INK,face="bold",size=10),
  plot.subtitle=element_text(colour=MUTED,size=8),
  legend.title=element_blank(), legend.text=element_text(colour=INK,size=7),
  legend.position="right", legend.key.height=unit(9,"pt"))
sv <- function(p,f,w=8.6,h=4.4) ggsave(file.path(FIG,f),p,width=w,height=h,device="pdf")
SC <- scale_colour_manual(values=PAL, breaks=LEV, drop=FALSE)

# ---- FIG 14: variance budget (AP). none is EXCLUDED - no noise floor there.
b <- V$budgets[V$budgets$response=="ap" & !is.nan(V$budgets$noise_floor), ]
bl <- do.call(rbind, lapply(c("noise_floor","structure_mains","structure_2way","structure_3way_plus"),
  function(k) data.frame(stratum=b$stratum, method=b$method, part=k, val=100*b[[k]])))
bl$part <- factor(bl$part, levels=c("noise_floor","structure_mains","structure_2way","structure_3way_plus"),
  labels=c("noise floor","main effects","2-way","3-way+"))
bl$stratum <- factor(SN[bl$stratum], levels=SN)
om <- b[b$stratum=="sparse_cont",]; ordm <- om$method[order(om$structure_mains)]
bl$method <- factor(bl$method, levels=ordm)
p <- ggplot(bl, aes(val, method, fill=part)) + geom_col(width=.72) +
  scale_fill_manual(values=c("noise floor"="#BBBBBB","main effects"="#0072B2",
                             "2-way"="#56B4E9","3-way+"="#E69F00")) +
  facet_wrap(~stratum) + th +
  labs(title="Sense V: the variance budget bounds every claim below it",
       x="% of cell-mean sum of squares", y=NULL,
       subtitle="Noise floor = irreducible. A large main-effect block means marginal plots are honest summaries.")
sv(p,"fig14_variance_budget.pdf",9,5)

# ---- FIG 15: do different RESPONSES have different drivers?
tm <- V$terms[V$terms$order==1 & V$terms$term %in% c("S","phi","region_size","enrich"), ]
tm <- tm[tm$stratum != "sparse_none", ]
ag <- aggregate(omega2 ~ term + response + stratum, tm, function(z) 100*mean(z,na.rm=TRUE))
ag$stratum <- factor(SN[ag$stratum], levels=SN)
ag$response <- factor(ag$response, levels=c("ap","fdr_at_90","total_mass_ratio"),
                      labels=c("AP (ranking)","FDR at 0.9","total mass ratio"))
ag$term <- factor(ag$term, levels=c("S","phi","region_size","enrich"))
p <- ggplot(ag, aes(omega2, term, fill=response)) +
  geom_col(position=position_dodge(width=.75), width=.7) +
  scale_fill_manual(values=c("AP (ranking)"="#0072B2","FDR at 0.9"="#D55E00","total mass ratio"="#E69F00")) +
  facet_wrap(~stratum) + th +
  labs(title="Sense V: accuracy and calibration are driven by DIFFERENT factors",
       x=expression("mean noise-corrected "*omega^2*" (%), averaged over methods"), y=NULL,
       subtitle="If the bars differed only in height the metrics would agree; they do not.")
sv(p,"fig15_sense_v_responses.pdf",9,3.8)

# ---- FIG 16: interaction involvement
so <- V$sobol[V$sobol$response=="ap" & V$sobol$stratum!="sparse_none", ]
so <- so[is.finite(so$interaction_involvement), ]
so$stratum <- factor(SN[so$stratum], levels=SN); so$g <- grp(so$method)
p <- ggplot(so, aes(interaction_involvement, factor, colour=g)) +
  geom_vline(xintercept=0, colour="#AAAAAA", linewidth=.4) +
  geom_point(size=1.7, alpha=.85) + SC + facet_wrap(~stratum) + th +
  labs(title="Sense V: how conditional is each factor's effect?",
       x=expression(S[Ti]-S[i]*"   (0 = acts additively; large = effect depends on other factors)"),
       y=NULL, subtitle="Near-zero means the marginal plot for that factor tells the whole story.")
sv(p,"fig16_interaction_involvement.pdf",9,3.4)

# ---- FIG 17: Sense D - who ever wins, and is it decided?
w <- D$winners; w$stratum <- factor(SN[w$stratum], levels=SN)
wc <- as.data.frame(table(w$winner, w$decided, w$stratum))
names(wc) <- c("method","decided","stratum","n"); wc <- wc[wc$n>0,]
wc$decided <- factor(ifelse(wc$decided=="TRUE","statistically decided","not decided"),
                     levels=c("not decided","statistically decided"))
tot <- tapply(wc$n, wc$method, sum); wc$method <- factor(wc$method, levels=names(sort(tot)))
p <- ggplot(wc, aes(n, method, fill=decided)) + geom_col(width=.72) +
  scale_fill_manual(values=c("not decided"="#DDDDDD","statistically decided"="#0072B2")) +
  facet_wrap(~stratum) + th +
  labs(title="Sense D: which methods are ever the best in a cell - and does it survive a test?",
       x="number of cells where the method scored highest", y=NULL,
       subtitle="Grey is the majority everywhere: the nominal winner usually cannot be distinguished from the runner-up.")
sv(p,"fig17_sense_d_winners.pdf",9,4.4)

# ---- FIG 18: margin distribution
w2 <- w[is.finite(w$dbar), ]
p <- ggplot(w2, aes(dbar, fill=decided)) +
  geom_histogram(bins=45, colour=NA) +
  scale_fill_manual(values=c(`FALSE`="#DDDDDD",`TRUE`="#0072B2"),
                    labels=c("not decided","decided"), name=NULL) +
  facet_wrap(~stratum, scales="free_y") + th +
  labs(title="Sense D: the winning margin is tiny almost everywhere",
       x="AP margin between best and runner-up in a cell", y="cells",
       subtitle="Cells only clear the paired test out in the right tail; the mass sits near zero.")
sv(p,"fig18_sense_d_margin.pdf",9,3.4)

# ---- FIG 19/20: ablation conditional on S and on phi
dl <- G$delta; dl$stratum <- factor(SN[dl$stratum], levels=SN)
dl$pair <- factor(dl$pair, levels=c("annot","joint"),
                  labels=c("annotations (FB - BEATRICE)","region sharing (FB-xregion - FB)"))
mk <- function(xv, xlab, fn) {
  a <- aggregate(delta ~ stratum + pair + get(xv), dl,
                 function(z) c(m=mean(z), se=sd(z)/sqrt(length(z))))
  out <- data.frame(stratum=a$stratum, pair=a$pair, x=a[[3]],
                    m=a$delta[,"m"], se=a$delta[,"se"])
  p <- ggplot(out, aes(x, m, colour=pair, group=pair)) +
    geom_hline(yintercept=0, colour="#AAAAAA", linewidth=.4) +
    geom_errorbar(aes(ymin=m-1.96*se, ymax=m+1.96*se), width=0, linewidth=.5) +
    geom_line(linewidth=.7) + geom_point(size=1.6) +
    scale_colour_manual(values=c("#D55E00","#0072B2")) +
    facet_wrap(~stratum) + th +
    labs(title=paste("Sense G: where the ablation gains come from -", xlab),
         x=xlab, y="paired difference in AP (95% CI)",
         subtitle="Zero line is no effect. The unannotated panel is the negative control.")
  sv(p, fn, 9, 3.5)
}
mk("S","number of causal variants S","fig19_ablation_by_S.pdf")
mk("phi","region heritability phi","fig20_ablation_by_phi.pdf")
cat("figs 14-20 written\n")
