`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
OUT <- "results/iter004_sparse/report"
L2 <- readRDS("results/iter004_sparse/combined_replicate_metrics.rds"); L2 <- L2[!L2$method %in% EXCLUDE_METHODS,]
L3 <- readRDS("results/iter004_sparse/combined_scenario_metrics.rds"); L3 <- L3[!L3$method %in% EXCLUDE_METHODS,]
esc <- function(s) gsub("_","\\\\_",s)
con <- file(file.path(OUT,"tables.tex"), "w")
w <- function(...) cat(..., "\n", sep="", file=con)

STR <- c(none="No annotations", binary="Binary annotations", continuous="Continuous annotations")
for (st in names(STR)) {
  d <- L3[L3$annotation_type==st,]
  ms <- sort(unique(d$method))
  f <- function(m,v,fn=mean) fn(d[[v]][d$method==m], na.rm=TRUE)
  tab <- do.call(rbind, lapply(ms, function(m) data.frame(
    method=m, ap=f(m,"ap"), fdr=f(m,"fdr_at_90"), mass=f(m,"total_mass_ratio"),
    ece=f(m,"ece_hi"), ss=f(m,"set_size_95"), cov=f(m,"set_hit_95"),
    res=f(m,"res"), bss=f(m,"bss"),
    ncell=sum(!is.na(d$ap[d$method==m])))))
  tab <- tab[order(-tab$ap),]
  w("\\begin{table}[htbp]\\centering\\small")
  w("\\caption{\\textbf{", STR[st], "} --- all methods, sparse arm. Ranked by mean AP. ",
    "\\emph{cells} is the number of non-missing cells out of 1{,}125; anything below that is a genuine absence, explained in \\S\\ref{sec:completeness}.}")
  w("\\label{tab:main-", st, "}")
  w("\\begin{tabular}{lrrrrrrrr r}\\toprule")
  w("method & AP & FDR@.9 & mass & ECE$_{hi}$ & $|CS_{95}|$ & cov$_{95}$ & RES & BSS & cells\\\\\\midrule")
  for (i in seq_len(nrow(tab))) {
    r <- tab[i,]
    bold <- grepl("^(fb_|functional_beatrice|beatrice)", r$method)
    nm <- esc(r$method); if (bold) nm <- paste0("\\textbf{", nm, "}")
    w(sprintf("%s & %.4f & %.4f & %.2f & %.4f & %.1f & %.3f & %.4f & %.3f & %d\\\\",
              nm, r$ap, r$fdr, r$mass, r$ece, r$ss, r$cov, r$res, r$bss, r$ncell))
  }
  w("\\bottomrule\\end{tabular}\\end{table}")
}

# ablation table
key <- paste(L2$job_dir,L2$S,L2$phi,L2$region_size,L2$iter,sep="|")
w("\\begin{table}[htbp]\\centering\\small")
w("\\caption{Paired ablation on AP. Each row compares two methods on the \\emph{same} simulated replicate; $n$ is the number of paired replicates. $t=\\bar d/\\mathrm{SE}(\\bar d)$.}\\label{tab:ablation}")
w("\\begin{tabular}{llrrrrr}\\toprule")
w("stratum & contrast & $\\bar d$ & sd & $n$ & SE & $t$\\\\\\midrule")
for (st in names(STR)) {
  sel <- L2$annotation_type==st
  sub <- function(m){v<-L2$ap[sel&L2$method==m];names(v)<-key[sel&L2$method==m];v}
  for (p3 in list(c("functional_beatrice","beatrice","annotations"),
                  c("fb_xregion","functional_beatrice","region sharing"),
                  c("fb_xregion","beatrice","combined"))) {
    a<-sub(p3[1]); b<-sub(p3[2]); ix<-intersect(names(a),names(b))
    d<-a[ix]-b[ix]; d<-d[is.finite(d)]; if(!length(d)) next
    se<-sd(d)/sqrt(length(d)); tv<-if(se>0) mean(d)/se else NA
    w(sprintf("%s & %s & %.4f & %.4f & %d & %.5f & %s\\\\", STR[st], p3[3],
              mean(d), sd(d), length(d), se,
              if (is.na(tv)) "---" else sprintf("%.1f", tv)))
  }
}
w("\\bottomrule\\end{tabular}\\end{table}")
close(con)
cat("tables.tex written\n")
