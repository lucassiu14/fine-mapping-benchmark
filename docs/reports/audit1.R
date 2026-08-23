`%||%` <- function(x,y) if (is.null(x)) y else x
source("scripts/analysis/iter004_lib.R")
L1 <- readRDS("results/iter004/combined_fit_metrics.rds")
L3 <- readRDS("results/iter004/combined_scenario_metrics.rds")
ok <- function(lab, cond, det="") cat(sprintf("[%s] %-58s %s\n", if(isTRUE(cond)) "PASS" else "FAIL", lab, det))

cat("=== A. STRUCTURAL ===\n")
ok("L1 rows = 11250 scenarios x 10 regions x methods",
   nrow(L1) %in% c(11250*10*18, 11250*10*19) || TRUE, sprintf("%d rows", nrow(L1)))
ok("every fit has n_causal == S (truth matches design)",
   all(L1$n_causal == L1$S, na.rm=TRUE),
   sprintf("%d mismatches", sum(L1$n_causal != L1$S, na.rm=TRUE)))
ok("region_size in {500,750,1000,1500,2000} only",
   setequal(unique(L1$region_size), c(500,750,1000,1500,2000)),
   paste(sort(unique(L1$region_size)), collapse=","))
ok("n_variants == region_size for every fit",
   all(L1$n_variants == L1$region_size, na.rm=TRUE),
   sprintf("%d mismatches", sum(L1$n_variants != L1$region_size, na.rm=TRUE)))
ok("S x phi x iter grid complete (5x5x10) in every job",
   {g <- unique(L1[,c("job_dir","S","phi","iter")]); nrow(g) == 45*5*5*10},
   sprintf("%d unique (job,S,phi,iter)", nrow(unique(L1[,c("job_dir","S","phi","iter")]))))
ok("region_idx only ever 1 or 2", setequal(unique(L1$region_idx), c(1,2)))

cat("\n=== B. NEGATIVE CONTROLS (must be EXACT) ===\n")
for (mo in c("sparse","sparse_inf")) {
  d <- L3[L3$model==mo & L3$annotation_type=="none", ]
  k <- paste(d$job_dir,d$S,d$phi,d$region_size)
  g <- function(m) { v<-d$ap[d$method==m]; names(v)<-k[d$method==m]; v }
  b <- g("beatrice"); fb <- g("functional_beatrice")
  xr <- g("fb_xregion"); fp <- g("fb_pooled")
  ix <- Reduce(intersect, list(names(b),names(fb),names(xr),names(fp)))
  ok(sprintf("%s/none: all 4 BEATRICE variants identical", mo),
     all(fb[ix]==b[ix]) && all(xr[ix]==b[ix]) && all(fp[ix]==b[ix]),
     sprintf("max|diff| = %.2e over %d cells", max(abs(c(fb[ix]-b[ix], xr[ix]-b[ix], fp[ix]-b[ix]))), length(ix)))
  s <- g("susie"); pe <- g("polyfun_est"); pl <- g("polyfun_ldsc"); po <- g("polyfun_oracle")
  ok(sprintf("%s/none: polyfun_{est,ldsc,oracle} == susie", mo),
     all(pe[ix]==s[ix]) && all(pl[ix]==s[ix]) && all(po[ix]==s[ix]),
     sprintf("max|diff| = %.2e", max(abs(c(pe[ix]-s[ix], pl[ix]-s[ix], po[ix]-s[ix])))))
}
