options(width=230, stringsAsFactors=FALSE)
sc <- readRDS("results/iter003/combined_scenario_metrics.rds")
sc$ld <- factor(ifelse(is.na(sc$n_ref),"in-sample",paste0("ref",sc$n_ref)),
                levels=c("in-sample","ref500","ref200"))
sc$enr <- ifelse(is.na(sc$enrichment_fold),"none",as.character(sc$enrichment_fold))
# comparator set: the joint models, their direct ancestors, and the baselines that
# actually earned attention in Iteration 002. Excludes the vacuous abstainers
# (abf, marginal_z), the dominated (susie_inf), and the collapsed (sbayesrc).
CORE <- c("fb_pooled","fb_xregion","functional_beatrice","beatrice",
          "finemap","susie","polyfun_est","funmap","polyfun_oracle")
d <- sc[sc$method %in% CORE,]
d$method <- factor(d$method, levels=CORE)

cat("=== is `ap` macro in this collect? (PR #50 made it so) ===\n")
cat("  mean |ap - ap_macro| =", signif(mean(abs(d$ap - d$ap_macro), na.rm=TRUE),3),
    " | mean(ap - ap_micro) =", signif(mean(d$ap - d$ap_micro, na.rm=TRUE),4), "\n\n")

med <- function(df, metric, by){
  a <- aggregate(df[[metric]], df[by], function(x) median(x, na.rm=TRUE))
  names(a)[ncol(a)] <- metric; a }
wide <- function(df, metric, idv){
  a <- med(df, metric, c(idv,"method"))
  w <- reshape(a, idvar=idv, timevar="method", direction="wide")
  names(w) <- gsub(paste0(metric,"."),"",names(w),fixed=TRUE); w }

ANN <- d[d$annotation_type %in% c("binary","continuous"),]

cat("################ 1. HEADLINE: did the joint prior fix calibration? ################\n")
for (metric in c("total_mass_ratio","hi_pip_reliab","max_fdr_violation_n20","ap")) {
  w <- wide(ANN, metric, c("annotation_type","ld"))
  keep <- intersect(CORE, names(w))
  w[keep] <- lapply(w[keep], function(z) round(z,3))
  cat("\n--- ", metric, " (median per arm x LD) ---\n", sep="")
  print(w[order(w$annotation_type,w$ld), c("annotation_type","ld",keep)], row.names=FALSE)
}

cat("\n\n################ 2. PAIRED within-cell win rates ################\n")
cat("  % of cells where the joint model beats its comparator, in the SAME cell.\n")
key <- c("job_dir","S","phi","region_size")
pw <- reshape(ANN[,c(key,"method","ap","total_mass_ratio","hi_pip_reliab","max_fdr_violation_n20")],
              idvar=key, timevar="method", direction="wide")
meta <- unique(sc[,c("job_dir","annotation_type","n_ref","enrichment_fold","model")])
pw <- merge(pw, meta, by="job_dir")
pw$ld <- factor(ifelse(is.na(pw$n_ref),"in-sample",paste0("ref",pw$n_ref)),
                levels=c("in-sample","ref500","ref200"))
pct <- function(x) round(100*mean(x, na.rm=TRUE),1)
for (j in c("fb_pooled","fb_xregion")) {
  cat("\n===", j, "vs functional_beatrice ===\n")
  t1 <- aggregate(cbind(
      AUPRC_better        = pw[[paste0("ap.",j)]] > pw[["ap.functional_beatrice"]],
      mass_closer_to_1    = abs(pw[[paste0("total_mass_ratio.",j)]]-1) <
                            abs(pw[["total_mass_ratio.functional_beatrice"]]-1),
      hiPIP_better        = pw[[paste0("hi_pip_reliab.",j)]] > pw[["hi_pip_reliab.functional_beatrice"]],
      FDR_better          = pw[[paste0("max_fdr_violation_n20.",j)]] <
                            pw[["max_fdr_violation_n20.functional_beatrice"]]
    ) ~ annotation_type + ld, pw, pct)
  print(t1, row.names=FALSE)
}
cat("\n=== fb_xregion vs fb_pooled (which head wins?) ===\n")
t2 <- aggregate(cbind(
    AUPRC_x_better = pw[["ap.fb_xregion"]] > pw[["ap.fb_pooled"]],
    mass_x_better  = abs(pw[["total_mass_ratio.fb_xregion"]]-1) < abs(pw[["total_mass_ratio.fb_pooled"]]-1),
    FDR_x_better   = pw[["max_fdr_violation_n20.fb_xregion"]] < pw[["max_fdr_violation_n20.fb_pooled"]]
  ) ~ annotation_type + ld, pw, pct)
print(t2, row.names=FALSE)

cat("\n\n################ 3. vs the best baseline (finemap) ################\n")
for (j in c("fb_pooled","fb_xregion","functional_beatrice")) {
  t3 <- aggregate(cbind(
      beats_finemap_AUPRC = pw[[paste0("ap.",j)]] > pw[["ap.finemap"]],
      beats_finemap_FDR   = pw[[paste0("max_fdr_violation_n20.",j)]] < pw[["max_fdr_violation_n20.finemap"]]
    ) ~ ld, pw, pct)
  cat("\n ", j, ":\n"); print(t3, row.names=FALSE)
}

cat("\n\n################ 4. enrichment response (does the shared prior USE annotations?) ################\n")
E <- ANN[!is.na(ANN$enrichment_fold),]
a <- med(E, "ap", c("annotation_type","ld","enrichment_fold","method"))
sl <- do.call(rbind, lapply(split(a, list(a$method,a$annotation_type,a$ld), drop=TRUE), function(z){
  if (nrow(z) < 3) return(NULL)
  data.frame(method=as.character(z$method[1]), annotation=z$annotation_type[1], ld=z$ld[1],
             slope=round(unname(coef(lm(ap ~ enrichment_fold, z))[2]),4))
}))
w <- reshape(sl, idvar=c("method","annotation"), timevar="ld", direction="wide")
names(w) <- gsub("slope.","",names(w),fixed=TRUE)
cat("  AUPRC change per unit enrichment fold:\n")
print(w[order(w$annotation, -w$`in-sample`),], row.names=FALSE)

cat("\n\n################ 5. by number of causals S (annotated, ref500) ################\n")
S5 <- ANN[ANN$ld=="ref500",]
w <- wide(S5, "ap", c("S")); keep <- intersect(CORE, names(w))
w[keep] <- lapply(w[keep], function(z) round(z,3))
print(w[order(as.numeric(w$S)), c("S",keep)], row.names=FALSE)
cat("\n  and mass ratio by S:\n")
w <- wide(S5, "total_mass_ratio", c("S")); keep <- intersect(CORE, names(w))
w[keep] <- lapply(w[keep], function(z) round(z,2))
print(w[order(as.numeric(w$S)), c("S",keep)], row.names=FALSE)
