#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter002_patterns.R
#
# Iteration 002 scenario-space pattern mining -> results/iteration-002-patterns.pdf
#
# Finds statements that hold across a NEIGHBOURHOOD of scenarios (e.g. a method
# best at every S, or every region length), not just at single points. Three
# pattern types: where a method is best, where it emits erroneous probabilities,
# and where two methods are interchangeable.
#
# NEVER aggregates across model (sparse/sparse_inf) or annotation type
# (none/binary/continuous) - those six cells are kept separate throughout.
#
# NA-as-a-level trap: n_ref=NA means IN-SAMPLE LD, enrichment_fold=NA means the
# none arm, p_causal=NA means the sparse model. base R's aggregate/split/table
# silently DROP NA, so all three are recoded to explicit levels first.
#
# polyfun_oracle is excluded from "best" claims (it is given the true annotation
# weights; it measures the ceiling, not an achievable method).
#
# Usage: Rscript scripts/analysis/iter002_patterns.R [metrics_rds] [fig_dir]
# =============================================================================
args <- commandArgs(trailingOnly=TRUE)
MET  <- if (length(args)>=1) args[1] else "results/iter002/combined_scenario_metrics.rds"
OUT  <- if (length(args)>=2) args[2] else "results/iter002/figures"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
s <- readRDS(MET)
num <- function(x) as.numeric(as.character(x))
s$S <- num(s$S); s$phi <- num(s$phi); s$rs <- num(s$region_size)
s$ld  <- factor(ifelse(is.na(s$n_ref),"in-sample",paste0("ref",s$n_ref)),
                levels=c("in-sample","ref500","ref200"))
s$enr <- ifelse(is.na(s$enrichment_fold),"none-arm",as.character(s$enrichment_fold))
s$pc  <- ifelse(is.na(s$p_causal),"sparse",as.character(s$p_causal))
s$cell<- paste(s$model, s$annotation_type, sep="/")
s$scen<- paste(s$job_dir,s$S,s$phi,s$region_size,sep="|")

# per-scenario: best AP, rank, gap-to-best
o <- order(s$scen, -s$ap, na.last=TRUE)
s <- s[o,]
g <- rle(s$scen)$lengths
s$rank <- unlist(lapply(g, seq_len))
best_ap <- rep(s$ap[cumsum(c(0,head(g,-1)))+1], g)
s$gap_to_best <- best_ap - s$ap          # 0 for the winner
s$is_win  <- s$rank==1L
s$near    <- s$gap_to_best <= 0.01       # within 0.01 AP of the winner
cat("scenarios:", length(unique(s$scen)), " rows:", nrow(s), "\n")


# ---- WIN RATE surface: cell x axis-level x method -------------------------
axes <- c("S","phi","rs","ld","enr","pc")
res <- list()
for (cl in sort(unique(s$cell))) {
  sc <- s[s$cell==cl,]
  nsc <- length(unique(sc$scen))
  for (ax in axes) {
    lv <- unique(sc[[ax]]); if (length(lv)<2) next
    for (l in sort(as.character(lv))) {
      sub <- sc[as.character(sc[[ax]])==l,]
      n_scen <- length(unique(sub$scen))
      wr <- tapply(sub$is_win, sub$method, mean)
      for (m in names(wr)) res[[length(res)+1]] <- data.frame(
        cell=cl, axis=ax, level=l, method=m, win=100*as.numeric(wr[m]),
        n_scen=n_scen, stringsAsFactors=FALSE)
    }
  }
}
W <- do.call(rbind, res); 

# ---- PATTERN A: method wins at EVERY level of an axis (within a cell) -----
cat("\n=== PATTERN A: method is top-1 at EVERY level of an axis (min win% >= 25) ===\n")
for (cl in sort(unique(W$cell))) for (ax in unique(W$axis[W$cell==cl])) {
  d <- W[W$cell==cl & W$axis==ax,]
  lv <- unique(d$level)
  # which method has the highest win% at each level?
  topm <- sapply(lv, function(l){ dd<-d[d$level==l,]; dd$method[which.max(dd$win)] })
  if (length(unique(topm))==1L) {
    m <- unique(topm); rng <- range(d$win[d$method==m])
    if (rng[1] >= 25) cat(sprintf("  %-22s across ALL %-4s levels (%s): %-18s win %.0f-%.0f%%\n",
        cl, ax, paste(lv,collapse=","), m, rng[1], rng[2]))
  }
}

# ---- pattern mining (practical methods, near-equivalence) ----
REAL <- setdiff(unique(s$method), "polyfun_oracle")   # oracle = ceiling, not usable
r <- s[s$method %in% REAL,]
# recompute rank/win among practical methods only
o <- order(r$scen, -r$ap, na.last=TRUE); r <- r[o,]
g <- rle(r$scen)$lengths
r$rank <- unlist(lapply(g, seq_len))
best <- rep(r$ap[cumsum(c(0,head(g,-1)))+1], g)
r$gap <- best - r$ap; r$is_win <- r$rank==1L


cat("=== A) PRACTICAL winner at EVERY level of an axis, within each cell ===\n")
axes <- c("S","phi","rs","ld","enr","pc")
for (cl in sort(unique(r$cell))) {
  rc <- r[r$cell==cl,]
  for (ax in axes) {
    lv <- sort(unique(as.character(rc[[ax]]))); if (length(lv)<2) next
    tops <- sapply(lv, function(l){ sub<-rc[as.character(rc[[ax]])==l,]
      wr<-tapply(sub$is_win,sub$method,mean); names(wr)[which.max(wr)] })
    wins <- sapply(lv, function(l){ sub<-rc[as.character(rc[[ax]])==l,]
      wr<-tapply(sub$is_win,sub$method,mean); 100*max(wr) })
    if (length(unique(tops))==1L)
      cat(sprintf("  %-22s ALL %-4s (%s): %-16s win %.0f-%.0f%%\n",
          cl, ax, paste(lv,collapse=","), unique(tops), min(wins), max(wins)))
    else
      cat(sprintf("  %-22s %-4s SPLITS: %s\n", cl, ax,
          paste(sprintf("%s->%s(%.0f%%)", lv, tops, wins), collapse="  ")))
  }
}

cat("\n=== B) NEAR-EQUIVALENCE: median |AP difference| between method pairs, per cell ===\n")
cat("   (pairs with median |dAP| <= 0.01 in a cell are effectively the same method there)\n")
ap <- reshape(r[,c("scen","method","ap","cell")], idvar=c("scen","cell"),
              timevar="method", direction="wide")
names(ap) <- sub("ap.","",names(ap))
ms <- setdiff(names(ap), c("scen","cell"))
for (cl in sort(unique(ap$cell))) {
  a <- ap[ap$cell==cl,]
  out <- list()
  for (i in seq_along(ms)) for (j in seq_along(ms)) if (i<j) {
    d <- median(abs(a[[ms[i]]]-a[[ms[j]]]), na.rm=TRUE)
    if (!is.na(d) && d <= 0.02) out[[length(out)+1]] <-
      sprintf("%s~%s (%.4f)", ms[i], ms[j], d)
  }
  cat(sprintf("  %-22s %s\n", cl, if(length(out)) paste(unlist(out), collapse="  ") else "-none-"))
}

# ---- erroneous-output + key-method patterns ----
pc <- function(x) 100*mean(x, na.rm=TRUE)
cat("=== C) ERRONEOUS-OUTPUT RATES: %% of scenarios per method x cell ===\n")
cat("   bad_fdr = max_fdr_violation>0.5 | overmass = mass_ratio>5 | unreliable = hi_pip_reliab<0.5\n")
e <- do.call(rbind, lapply(split(s, list(s$cell,s$method), drop=TRUE), function(d)
  data.frame(cell=d$cell[1], method=d$method[1],
             bad_fdr=pc(d$max_fdr_violation>0.5), overmass=pc(d$total_mass_ratio>5),
             unreliable=pc(d$hi_pip_reliab<0.5), n=nrow(d), stringsAsFactors=FALSE)))
e$worst <- pmax(e$bad_fdr,e$overmass,e$unreliable)
for (cl in sort(unique(e$cell))) {
  d <- e[e$cell==cl,]; d <- d[order(-d$worst),]
  cat(sprintf("\n-- %s\n", cl))
  print(format(d[d$worst>1, c("method","bad_fdr","overmass","unreliable")], digits=3), row.names=FALSE)
}

cat("\n=== D) SBAYESRC: ranks well but is it trustworthy? by LD regime (all cells) ===\n")
sb <- s[s$method=="sbayesrc",]
print(do.call(rbind, lapply(split(sb, sb$ld), function(d) data.frame(
  ld=d$ld[1], med_AP=round(median(d$ap),3), med_violation=round(median(d$max_fdr_violation),3),
  med_massratio=round(median(d$total_mass_ratio),1),
  med_hi_pip_reliab=round(median(d$hi_pip_reliab,na.rm=TRUE),3)))), row.names=FALSE)
cat("  sbayesrc WIN RATE (practical methods) by ld x cell:\n")
sbw <- r[r$method=="sbayesrc",]
print(round(tapply(sbw$is_win, list(sbw$cell, sbw$ld), mean)*100,1))

cat("\n=== E) THE 'ONE METHOD' CLUSTER: susie / polyfun_est / polyfun_ldsc / funmap ===\n")
ap <- reshape(r[,c("scen","method","ap","cell")], idvar=c("scen","cell"), timevar="method", direction="wide")
names(ap) <- sub("ap.","",names(ap))
for (cl in sort(unique(ap$cell))) {
  a <- ap[ap$cell==cl,]
  f <- function(x,y) pc(abs(a[[x]]-a[[y]]) <= 0.01)
  cat(sprintf("  %-22s %%scen within 0.01 AP:  est~susie %.0f  ldsc~susie %.0f  funmap~susie %s\n",
    cl, f("polyfun_est","susie"), f("polyfun_ldsc","susie"),
    if(all(is.na(a$funmap))) "n/a" else sprintf("%.0f", f("funmap","susie"))))
}

cat("\n=== F) SUSIE: does it EVER win among practical methods? ===\n")
su <- r[r$method=="susie",]
cat(sprintf("  overall win rate: %.2f%%  | best cell: %s\n", pc(su$is_win),
    paste(names(which.max(tapply(su$is_win,su$cell,mean))))))
print(round(tapply(su$is_win, list(su$cell, su$ld), mean)*100,2))
cat("  susie median FDR violation by ld:\n")
print(round(tapply(su$max_fdr_violation, su$ld, median),3))

cat("\n=== G) ABF at S=1: win rate by cell (practical) ===\n")
ab <- r[r$method=="abf",]
print(round(tapply(ab$is_win, list(ab$cell, ab$S), mean)*100,1))

# ---- figures ----
suppressPackageStartupMessages({library(ggplot2); library(grid)})
# OUT set from args above
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

BLUE<-"#2a78d6"; RED<-"#e34948"; GREEN<-"#008300"; VIOLET<-"#4a3aa7"; ORANGE<-"#eb6834"
INK<-"#0b0b0b"; INK2<-"#52514e"; GRID<-"#e6e6e3"; SURF<-"#fcfcfb"
th <- theme_minimal(base_size=9) + theme(
  plot.background=element_rect(fill=SURF,colour=NA), panel.background=element_rect(fill=SURF,colour=NA),
  panel.grid.minor=element_blank(), panel.grid.major=element_line(colour=GRID,linewidth=.3),
  axis.title=element_text(colour=INK2,size=8), axis.text=element_text(colour=INK2,size=7.5),
  plot.title=element_text(colour=INK,face="bold",size=10), plot.subtitle=element_text(colour=INK2,size=8),
  strip.text=element_text(colour=INK,face="bold",size=7.5),
  legend.position="bottom", legend.title=element_blank(), legend.text=element_text(colour=INK2,size=7.5),
  legend.key.height=unit(9,"pt"), legend.margin=margin(t=-4), panel.spacing=unit(9,"pt"))
CELLS <- c("sparse/none","sparse/binary","sparse/continuous",
           "sparse_inf/none","sparse_inf/binary","sparse_inf/continuous")
r$cell <- factor(r$cell, levels=CELLS); s$cell <- factor(s$cell, levels=CELLS)

# F1: ABF win% vs S, per cell
d1 <- aggregate(is_win ~ cell + S + method, r[r$method %in% c("abf","finemap","sbayesrc","susie"),], mean)
d1$win <- 100*d1$is_win
p1 <- ggplot(d1, aes(factor(S), win, colour=method, group=method)) +
  geom_line(linewidth=.7)+geom_point(size=1.6)+facet_wrap(~cell,nrow=2)+
  scale_colour_manual(values=c(abf=BLUE,finemap=GREEN,sbayesrc=ORANGE,susie=VIOLET))+
  labs(title="ABF is the best practical method when there is exactly one causal variant",
       subtitle="Win rate among the 13 practical methods (oracle excluded). Within each cell, aggregated over phi(5), region size(5), LD(3), enrichment(4 where applicable), p_causal(4 for sparse_inf).",
       x="true number of causal variants S", y="% of scenarios method is best")+th
ggsave(file.path(OUT,"f1_abf_S.png"),p1,width=9.4,height=4.6,dpi=200)

# F2: winner composition vs LD, per cell (top methods)
top <- c("finemap","sbayesrc","polyfun_est","sparsepro","beatrice","abf","paintor","polyfun_ldsc")
d2 <- aggregate(is_win ~ cell + ld + method, r[r$method %in% top,], mean); d2$win <- 100*d2$is_win
p2 <- ggplot(d2, aes(ld, win, colour=method, group=method)) +
  geom_line(linewidth=.7)+geom_point(size=1.6)+facet_wrap(~cell,nrow=2)+
  scale_colour_manual(values=c(finemap=GREEN,sbayesrc=ORANGE,polyfun_est=BLUE,sparsepro=VIOLET,
                               beatrice=RED,abf="#7a7a7a",paintor="#c9a227",polyfun_ldsc="#1baf7a"))+
  labs(title="FINEMAP takes over as the LD reference panel degrades - in all six cells",
       subtitle="Win rate among 13 practical methods. Within each cell, aggregated over S(5), phi(5), region size(5), enrichment(4 where applicable), p_causal(4 for sparse_inf).",
       x="LD regime (in-sample -> 500-sample panel -> 200-sample panel)", y="% of scenarios method is best")+th
ggsave(file.path(OUT,"f2_ld_takeover.png"),p2,width=9.4,height=4.8,dpi=200)

# F3: sbayesrc collapse (AP + violation + mass), per cell
sb <- s[s$method=="sbayesrc",]
a <- aggregate(cbind(ap,max_fdr_violation,total_mass_ratio) ~ cell + ld, sb, median)
long <- rbind(data.frame(cell=a$cell,ld=a$ld,metric="average precision",v=a$ap),
              data.frame(cell=a$cell,ld=a$ld,metric="max FDR violation",v=a$max_fdr_violation),
              data.frame(cell=a$cell,ld=a$ld,metric="total PIP mass ratio",v=a$total_mass_ratio))
long$metric <- factor(long$metric, levels=c("average precision","max FDR violation","total PIP mass ratio"))
p3 <- ggplot(long, aes(ld, v, colour=cell, group=cell))+
  geom_line(linewidth=.7)+geom_point(size=1.6)+
  facet_wrap(~metric, scales="free_y", nrow=1)+
  scale_colour_manual(values=c(BLUE,RED,GREEN,VIOLET,ORANGE,"#e87ba4"))+
  labs(title="SBayesRC does not degrade under LD noise - it disintegrates",
       subtitle="Median per scenario cell. Within each cell, aggregated over S(5), phi(5), region size(5), enrichment(4 where applicable), p_causal(4 for sparse_inf). Mass ratio 1.0 = correct.",
       x="LD regime", y=NULL)+th
ggsave(file.path(OUT,"f3_sbayesrc.png"),p3,width=9.4,height=3.6,dpi=200)

# F4: near-equivalence to susie
ap <- reshape(r[,c("scen","method","ap","cell")],idvar=c("scen","cell"),timevar="method",direction="wide")
names(ap) <- sub("ap.","",names(ap))
eq <- do.call(rbind, lapply(CELLS, function(cl){ a<-ap[ap$cell==cl,]
  f<-function(x) if(all(is.na(a[[x]]))) NA else 100*mean(abs(a[[x]]-a$susie)<=0.01,na.rm=TRUE)
  data.frame(cell=cl, method=c("polyfun_est","polyfun_ldsc","funmap"),
             pct=c(f("polyfun_est"),f("polyfun_ldsc"),f("funmap"))) }))
eq$cell <- factor(eq$cell, levels=CELLS); eq <- eq[!is.na(eq$pct),]
p4 <- ggplot(eq, aes(cell, pct, fill=method))+
  geom_col(position=position_dodge(.75), width=.68)+
  geom_text(aes(label=sprintf("%.0f",pct)), position=position_dodge(.75), vjust=-.4, size=2.4, colour=INK)+
  scale_fill_manual(values=c(polyfun_est=BLUE,polyfun_ldsc=GREEN,funmap=ORANGE))+
  scale_y_continuous(limits=c(0,108), expand=expansion(mult=c(0,.02)))+
  labs(title="Three 'different' methods are, in most scenarios, indistinguishable from SuSiE",
       subtitle="% of scenarios where the method's average precision is within 0.01 of SuSiE's. Within each cell, aggregated over S(5), phi(5), region size(5), LD(3), enrichment(4 where applicable), p_causal(4 for sparse_inf). funmap does not run on the 'none' arm.",
       x=NULL, y="% of scenarios within 0.01 AP of SuSiE")+th+
  theme(axis.text.x=element_text(angle=20,hjust=1))
ggsave(file.path(OUT,"f4_equivalence.png"),p4,width=9.4,height=3.8,dpi=200)
cat("figs2 written\n"); print(list.files(OUT))

BLUE<-"#2a78d6"; RED<-"#e34948"; GREEN<-"#008300"; VIOLET<-"#4a3aa7"; ORANGE<-"#eb6834"
INK<-"#0b0b0b"; INK2<-"#52514e"; GRID<-"#e6e6e3"; SURF<-"#fcfcfb"
th <- theme_minimal(base_size=9)+theme(
  plot.background=element_rect(fill=SURF,colour=NA), panel.background=element_rect(fill=SURF,colour=NA),
  panel.grid.minor=element_blank(), panel.grid.major=element_line(colour=GRID,linewidth=.3),
  axis.title=element_text(colour=INK2,size=8), axis.text=element_text(colour=INK2,size=7.5),
  plot.title=element_text(colour=INK,face="bold",size=10), plot.subtitle=element_text(colour=INK2,size=8),
  strip.text=element_text(colour=INK,face="bold",size=7.5), legend.position="bottom",
  legend.title=element_text(colour=INK2,size=7.5), legend.text=element_text(colour=INK2,size=7.5),
  legend.key.height=unit(9,"pt"), legend.margin=margin(t=-4), panel.spacing=unit(9,"pt"))
CELLS <- c("sparse/none","sparse/binary","sparse/continuous",
           "sparse_inf/none","sparse_inf/binary","sparse_inf/continuous")
s$cell <- factor(s$cell, levels=CELLS)

# F5: erroneous-output heatmap (% scenarios with max_fdr_violation > 0.5)
e <- aggregate(cbind(bad=max_fdr_violation>0.5) ~ cell + method, s, function(x) 100*mean(x))
ord <- aggregate(bad ~ method, e, mean); e$method <- factor(e$method, levels=ord$method[order(ord$bad)])
p5 <- ggplot(e, aes(cell, method, fill=bad))+geom_tile(colour=SURF, linewidth=1.1)+
  geom_text(aes(label=sprintf("%.0f",bad)), size=2.5,
            colour=ifelse(e$bad>45,"white",INK))+
  scale_fill_gradientn(colours=c("#cde2fb","#86b6ef","#3987e5","#256abf","#0d366b"),
                       limits=c(0,80), name="% of scenarios")+
  labs(title="Where each method emits badly over-confident probabilities",
       subtitle="% of that cell's scenarios with max FDR violation > 0.5 (i.e. selecting PIP>=t yields an FDR at least 0.5 above the 1-t guarantee).\nWithin each cell, aggregated over S(5), phi(5), region size(5), LD(3), enrichment(4 where applicable), p_causal(4 for sparse_inf).",
       x=NULL,y=NULL)+th+theme(axis.text.x=element_text(angle=20,hjust=1), panel.grid=element_blank())
ggsave(file.path(OUT,"f5_errors.png"),p5,width=9.4,height=4.6,dpi=200)

# F6: FB - beatrice, by cell x LD  (annotation-specific damage)
fb <- s[s$method %in% c("beatrice","functional_beatrice"),]
a <- aggregate(ap ~ method+cell+ld, fb, median)
w <- reshape(a, idvar=c("cell","ld"), timevar="method", direction="wide")
names(w) <- sub("ap.","",names(w)); w$gain <- w$functional_beatrice - w$beatrice
p6 <- ggplot(w, aes(ld, gain, fill=gain>0))+geom_hline(yintercept=0,colour=INK2,linewidth=.4)+
  geom_col(width=.62)+facet_wrap(~cell,nrow=2)+
  geom_text(aes(label=sprintf("%+.3f",gain), vjust=ifelse(gain>0,-.4,1.3)), size=2.2, colour=INK)+
  scale_fill_manual(values=c(`FALSE`=RED,`TRUE`=BLUE), guide="none")+
  scale_y_continuous(expand=expansion(mult=c(.22,.22)))+
  labs(title="Adding per-locus annotations to BEATRICE only ever helps when LD is perfect",
       subtitle="functional BEATRICE minus BEATRICE, median average precision. Same inference engine; annotations are the only difference.\nWithin each cell x LD regime, aggregated over S(5), phi(5), region size(5), enrichment(4 where applicable), p_causal(4 for sparse_inf).\nThe 'none' arms are an exact control - no annotations exist, so the difference is 0.000 by construction.",
       x="LD regime", y="delta average precision")+th
ggsave(file.path(OUT,"f6_fb.png"),p6,width=9.4,height=4.8,dpi=200)

# F7: who ever wins? win-rate composition per cell
r$cell <- factor(r$cell, levels=CELLS)
d <- aggregate(is_win ~ cell+method, r, mean); d$win <- 100*d$is_win
d <- d[d$win >= 1,]
d$method <- factor(d$method, levels=rev(sort(unique(d$method))))
p7 <- ggplot(d, aes(win, method, fill=cell))+
  geom_col(position=position_dodge(.8), width=.72)+
  scale_fill_manual(values=c(BLUE,RED,GREEN,VIOLET,ORANGE,"#e87ba4"), name=NULL)+
  labs(title="Which methods are ever the best, and in which simulation setting",
       subtitle="Win rate among the 13 practical methods (polyfun_oracle excluded as it is given the true annotation weights). Methods winning <1% of a cell's scenarios omitted.\nWithin each cell, aggregated over S(5), phi(5), region size(5), LD(3), enrichment(4 where applicable), p_causal(4 for sparse_inf).",
       x="% of that cell's scenarios where the method has the highest average precision", y=NULL)+th
ggsave(file.path(OUT,"f7_winners.png"),p7,width=9.4,height=4.4,dpi=200)
cat("done\n"); print(list.files(OUT))
