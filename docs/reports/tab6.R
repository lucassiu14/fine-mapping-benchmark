`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
OUT <- "results/iter004_full/report"
L2 <- readRDS("results/iter004/combined_replicate_metrics.rds"); L2<-L2[!L2$method %in% EXCLUDE_METHODS,]
L3 <- readRDS("results/iter004/combined_scenario_metrics.rds"); L3<-L3[!L3$method %in% EXCLUDE_METHODS,]
esc <- function(s) gsub("_","\\\\_",s)
con <- file(file.path(OUT,"tables.tex"),"w"); w <- function(...) cat(...,"\n",sep="",file=con)
STR <- list(c("sparse","none","Sparse / no annotations"), c("sparse","binary","Sparse / binary"),
            c("sparse","continuous","Sparse / continuous"),
            c("sparse_inf","none","Sparse+infinitesimal / no annotations"),
            c("sparse_inf","binary","Sparse+infinitesimal / binary"),
            c("sparse_inf","continuous","Sparse+infinitesimal / continuous"))
for (S in STR) {
  d <- L3[L3$model==S[1] & L3$annotation_type==S[2],]
  ncell <- length(unique(paste(d$job_dir,d$S,d$phi,d$region_size)))
  ms <- sort(unique(d$method)); f <- function(m,v) mean(d[[v]][d$method==m],na.rm=TRUE)
  tab <- do.call(rbind,lapply(ms,function(m) data.frame(method=m,ap=f(m,"ap"),
      fdr=f(m,"fdr_at_90"),mass=f(m,"total_mass_ratio"),ece=f(m,"ece_hi"),
      ss=f(m,"set_size_95"),cov=f(m,"set_hit_95"),res=f(m,"res"),
      n=sum(!is.na(d$ap[d$method==m])))))
  tab <- tab[order(-tab$ap),]
  w("\\begin{table}[htbp]\\centering\\small")
  w("\\caption{\\textbf{",S[3],"} --- all 18 analysed methods, ranked by mean AP. ",
    "\\emph{cells} is non-missing cells out of ",format(ncell,big.mark="{,}"),".}")
  w("\\label{tab:main-",S[1],"-",S[2],"}")
  w("\\begin{tabular}{lrrrrrrr r}\\toprule")
  w("method & AP & FDR@.9 & mass & ECE$_{hi}$ & $|CS_{95}|$ & cov$_{95}$ & RES & cells\\\\\\midrule")
  for (i in seq_len(nrow(tab))) { r<-tab[i,]
    nm <- esc(r$method)
    if (grepl("^(fb_|functional_beatrice|beatrice)",r$method)) nm <- paste0("\\textbf{",nm,"}")
    w(sprintf("%s & %.4f & %.4f & %.2f & %.4f & %.1f & %.3f & %.4f & %d\\\\",
              nm,r$ap,r$fdr,r$mass,r$ece,r$ss,r$cov,r$res,r$n)) }
  w("\\bottomrule\\end{tabular}\\end{table}") }

key <- paste(L2$job_dir,L2$S,L2$phi,L2$region_size,L2$iter,sep="|")
w("\\begin{table}[htbp]\\centering\\small")
w("\\caption{Paired ablation on AP, all six strata. Each row compares two methods on the \\emph{same} simulated replicate.}\\label{tab:ablation}")
w("\\begin{tabular}{llrrrr}\\toprule")
w("stratum & contrast & $\\bar d$ & sd & $n$ & $t$\\\\\\midrule")
for (S in STR) {
  s <- L2$model==S[1] & L2$annotation_type==S[2]
  sub <- function(m){v<-L2$ap[s&L2$method==m];names(v)<-key[s&L2$method==m];v}
  for (pp in list(c("functional_beatrice","beatrice","annotations"),
                  c("fb_xregion","functional_beatrice","region sharing"),
                  c("fb_xregion","beatrice","combined"))) {
    x<-sub(pp[1]);y<-sub(pp[2]);ix<-intersect(names(x),names(y))
    d<-x[ix]-y[ix];d<-d[is.finite(d)];if(!length(d))next
    se<-sd(d)/sqrt(length(d))
    w(sprintf("%s & %s & %+.4f & %.4f & %s & %s\\\\", S[3], pp[3], mean(d), sd(d),
              format(length(d),big.mark="{,}"),
              if(is.na(se)||se==0) "---" else sprintf("%.1f",mean(d)/se))) } }
w("\\bottomrule\\end{tabular}\\end{table}")
close(con); cat("tables.tex written\n")
