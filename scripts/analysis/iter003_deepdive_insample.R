#!/usr/bin/env Rscript
# =============================================================================
# Part 2: (a) FIXED robustness table - each metric aggregated separately, because
#             cbind() in aggregate() uses complete cases across all responses and
#             hi_pip_reliab is NA in ~26% of cells (abstention). That silently
#             dropped the hardest cells and inflated AP, unequally by method
#             (beatrice 0.493 -> 0.678, finemap 0.482 -> 0.569).
#         (b) IN-SAMPLE REGIME deep dive - methods that fail under a reference
#             panel are NOT written off; they get a full profile in the regime
#             where they do work.
# =============================================================================
options(width=250, stringsAsFactors=FALSE)
TAB <- "/private/tmp/claude-502/-Users-ls1020/cdc231a7-d353-43cf-ac37-fd155cce1deb/scratchpad/tabs_ultra"
W <- function(x,f){ write.csv(x, file.path(TAB,f), row.names=FALSE); invisible(x) }
sc <- readRDS("results/iter003/combined_scenario_metrics.rds")
sc$ld  <- factor(ifelse(is.na(sc$n_ref),"in-sample",paste0("ref",sc$n_ref)),
                 levels=c("in-sample","ref500","ref200"))
sc$ann <- factor(sc$annotation_type, levels=c("none","binary","continuous"))
FOCUS <- c("susie","abf","sbayesrc","finemap","funmap","polyfun_est","polyfun_ldsc",
           "beatrice","functional_beatrice","fb_xregion","susie_inf")
CTX   <- c("polyfun_oracle","sparsepro","paintor","marginal_z","fb_pooled")
ALLM  <- c(FOCUS,CTX); DEPLOY <- setdiff(ALLM,"polyfun_oracle")
d <- sc[sc$method %in% ALLM,]
med1 <- function(df, metric, by) {                    # ONE metric at a time
  a <- aggregate(df[[metric]], df[by], function(x) median(x, na.rm=TRUE))
  names(a)[ncol(a)] <- metric; a }

# ---------------- (a) robustness, corrected -----------------------------------
cat("=== ROBUSTNESS (corrected: per-metric aggregation, no complete-case filter) ===\n")
ap <- med1(d,"ap",c("method","ld")); re <- med1(d,"hi_pip_reliab",c("method","ld"))
w1 <- reshape(ap, idvar="method", timevar="ld", direction="wide"); names(w1) <- gsub("ap.","AP_",names(w1),fixed=TRUE)
w2 <- reshape(re, idvar="method", timevar="ld", direction="wide"); names(w2) <- gsub("hi_pip_reliab.","REL_",names(w2),fixed=TRUE)
rob <- merge(w1,w2,by="method")
rob$AP_retained  <- round(rob$AP_ref200 / rob$`AP_in-sample`, 3)
rob$REL_retained <- round(rob$REL_ref200/ rob$`REL_in-sample`,3)
n <- sapply(rob,is.numeric); rob[n] <- lapply(rob[n], function(z) round(z,3))
rob <- rob[order(-rob$AP_retained),]
W(rob,"robustness.csv"); print(rob,row.names=FALSE)

# ---------------- (b) IN-SAMPLE REGIME: full standing -------------------------
cat("\n\n############### IN-SAMPLE LD REGIME (methods judged on their own terms) ###############\n")
ins <- d[d$ld=="in-sample",]
key <- c("job_dir","S","phi","region_size")
dep <- ins[ins$method %in% DEPLOY,]
cb <- aggregate(ap ~ job_dir+S+phi+region_size, dep, function(x) max(x,na.rm=TRUE))
names(cb)[5] <- "best_ap"; dep <- merge(dep, cb, by=key)
dep$rel <- dep$ap/pmax(dep$best_ap,1e-9); dep$best <- abs(dep$ap-dep$best_ap)<1e-12
dep$near <- dep$rel>=0.95
dep$trust <- !is.na(dep$max_fdr_violation_n20) & dep$max_fdr_violation_n20<=0.05
for (arm in c("none","binary","continuous")) {
  z <- dep[dep$ann==arm,]
  if (!nrow(z)) next
  a <- aggregate(cbind(ap,rel,best,near,trust) ~ method, z, function(x) mean(x,na.rm=TRUE))
  m <- med1(z,"ap","method"); names(m)[2] <- "median_AP"
  a <- merge(a[,c("method","rel","best","near","trust")], m, by="method")
  a$rel <- round(a$rel,3); a[c("best","near","trust")] <- lapply(a[c("best","near","trust")], function(x) round(100*x,1))
  a$median_AP <- round(a$median_AP,3)
  a <- a[order(-a$median_AP), c("method","median_AP","rel","best","near","trust")]
  names(a) <- c("method","median_AP","rel_AP","pct_best","pct_within5","pct_FDRok")
  cat("\n--- IN-SAMPLE,", arm, "annotations ---\n"); print(a,row.names=FALSE)
  W(a, paste0("insample_standing_",arm,".csv"))
}

# ---------------- (b2) in-sample: where does each method peak? ----------------
cat("\n\n=== IN-SAMPLE: each method's BEST niche (highest rel_AP cell-group) ===\n")
dep$Sf <- as.numeric(dep$S); dep$phif <- as.numeric(dep$phi); dep$rsf <- as.numeric(dep$region_size)
niche <- NULL
for (m in DEPLOY) {
  z <- dep[dep$method==m,]
  if (!nrow(z)) next
  g <- aggregate(rel ~ ann + Sf + phif, z, function(x) mean(x,na.rm=TRUE))
  g <- g[order(-g$rel),]
  gw <- g[nrow(g),]
  niche <- rbind(niche, data.frame(method=m,
      best = sprintf("%s, S=%d, phi=%.4g", g$ann[1], g$Sf[1], g$phif[1]),
      best_rel = round(g$rel[1],3),
      worst = sprintf("%s, S=%d, phi=%.4g", gw$ann, gw$Sf, gw$phif),
      worst_rel = round(gw$rel,3)))
}
niche <- niche[order(-niche$best_rel),]
W(niche,"insample_niches.csv"); print(niche,row.names=FALSE)

# ---------------- (b3) in-sample: sensitivity to S / phi / region -------------
cat("\n\n=== IN-SAMPLE: median AP by S, by phi, by region length ===\n")
for (v in c("S","phi","region_size")) {
  a <- med1(ins[ins$method %in% FOCUS,], "ap", c("method",v))
  a[[v]] <- as.numeric(a[[v]])
  w <- reshape(a, idvar="method", timevar=v, direction="wide")
  names(w) <- gsub("ap.","",names(w),fixed=TRUE)
  num <- sapply(w,is.numeric); w[num] <- lapply(w[num], function(z) round(z,3))
  w <- w[match(FOCUS,w$method),]
  cat("\n-- by", v, "--\n"); print(w,row.names=FALSE)
  W(w, paste0("insample_by_",v,".csv"))
}

# ---------------- (c) FDR at practical thresholds, per LD regime --------------
cat("\n\n=== FDR at practical selection thresholds (from the true curves) ===\n")
jm <- unique(sc[,c("job_dir","annotation_type","n_ref")])
fdr <- readRDS("results/iter003/combined_fdr_curves.rds")
fdr <- fdr[fdr$method %in% ALLM & (abs(fdr$threshold-0.90)<1e-9 | abs(fdr$threshold-0.95)<1e-9),]
fdr <- merge(fdr, jm, by="job_dir")
fdr$ld <- factor(ifelse(is.na(fdr$n_ref),"in-sample",paste0("ref",fdr$n_ref)),
                 levels=c("in-sample","ref500","ref200"))
res <- NULL
for (L in levels(fdr$ld)) for (t in c(0.90,0.95)) {
  g <- fdr[fdr$ld==L & abs(fdr$threshold-t)<1e-9,]
  a <- aggregate(cbind(tp,fp)~method, g, sum)
  a$n_sel <- a$tp+a$fp; a$emp_FDR <- round(a$fp/pmax(a$n_sel,1),4)
  a$promised <- 1-t; a$excess <- round(a$emp_FDR-(1-t),4); a$ld <- L; a$threshold <- t
  res <- rbind(res, a[,c("ld","threshold","method","n_sel","emp_FDR","promised","excess")])
}
W(res,"fdr_at_thresholds.csv")
for (L in c("in-sample","ref500")) {
  z <- res[res$ld==L & res$threshold==0.95,]; z <- z[order(z$emp_FDR),]
  cat("\n-- LD =",L,", PIP>=0.95 (promise FDR<=0.05) --\n"); print(z[,c("method","n_sel","emp_FDR","excess")],row.names=FALSE)
}
cat("\nDONE\n")
