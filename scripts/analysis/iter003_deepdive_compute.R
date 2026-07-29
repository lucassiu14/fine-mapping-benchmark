#!/usr/bin/env Rscript
# =============================================================================
# Ultra deep-dive: per-method profiles + comparative analyses.
# Writes every table to tabs_ultra/ as CSV for the report builder.
# =============================================================================
options(width=250, stringsAsFactors=FALSE)
TAB <- "/private/tmp/claude-502/-Users-ls1020/cdc231a7-d353-43cf-ac37-fd155cce1deb/scratchpad/tabs_ultra"
dir.create(TAB, showWarnings=FALSE, recursive=TRUE)
W <- function(x,f){ write.csv(x, file.path(TAB,f), row.names=FALSE); invisible(x) }

sc <- readRDS("results/iter003/combined_scenario_metrics.rds")

# ---- recode every variable so NA-as-a-level is explicit (never silently dropped)
sc$ld    <- factor(ifelse(is.na(sc$n_ref),"in-sample",paste0("ref",sc$n_ref)),
                   levels=c("in-sample","ref500","ref200"))
sc$pc    <- factor(ifelse(is.na(sc$p_causal),"sparse-model",as.character(sc$p_causal)),
                   levels=c("sparse-model","0.5","0.7","0.9","1"))
sc$enr   <- factor(ifelse(is.na(sc$enrichment_fold),"none",as.character(sc$enrichment_fold)),
                   levels=c("none","2.7","5.4","8.1","10.8"))
sc$ann   <- factor(sc$annotation_type, levels=c("none","binary","continuous"))
sc$mdl   <- factor(sc$model, levels=c("sparse","sparse_inf"))
sc$Sf    <- factor(as.numeric(sc$S), levels=c(1,2,3,5,10))
sc$phif  <- factor(as.numeric(sc$phi), levels=c(0.0075,0.05,0.1,0.2,0.4))
sc$rsf   <- factor(as.numeric(sc$region_size), levels=c(100,200,400,500,1000))

# The 10 methods the user asked for, plus context comparators.
FOCUS <- c("susie","abf","sbayesrc","finemap","funmap","polyfun_est","polyfun_ldsc",
           "beatrice","functional_beatrice","fb_xregion","susie_inf")
CONTEXT <- c("polyfun_oracle","sparsepro","paintor","marginal_z","fb_pooled")
ALLM <- c(FOCUS, CONTEXT)
# "deployable" = excludes the oracle (given true weights) - used for win rates.
DEPLOY <- setdiff(ALLM, "polyfun_oracle")
d <- sc[sc$method %in% ALLM, ]
cat("cells:", nrow(d), " methods:", length(unique(d$method)), "\n")

VARS <- c(mdl="effect model", ann="annotation type", enr="annotation enrichment",
          ld="LD condition", Sf="number of causals S", phif="heritability phi",
          rsf="region length", pc="p_causal")
METRICS <- c(ap="macro-AUPRC", total_mass_ratio="total mass ratio",
             hi_pip_reliab="high-PIP reliability", max_fdr_violation_n20="guarded FDR violation")

# =============================================================================
# 1. VARIANCE DECOMPOSITION - which variables actually move each method?
#    eta^2 = SS_term / SS_total from a main-effects ANOVA fitted PER METHOD.
#    This is a decomposition of variability, not a pooled performance claim:
#    all performance numbers elsewhere stay stratified.
# =============================================================================
cat("\n=== 1. variance decomposition (eta^2, %) ===\n")
eta_for <- function(dd, metric) {
  dd <- dd[is.finite(dd[[metric]]), ]
  # drop variables that do not vary for this method/metric subset
  use <- names(VARS)[vapply(names(VARS), function(v) nlevels(droplevels(dd[[v]])) > 1, logical(1))]
  if (!length(use) || nrow(dd) < 50) return(NULL)
  f <- stats::as.formula(paste(metric, "~", paste(use, collapse=" + ")))
  a <- tryCatch(stats::anova(stats::lm(f, data=dd)), error=function(e) NULL)
  if (is.null(a)) return(NULL)
  ss <- a[["Sum Sq"]]; names(ss) <- rownames(a)
  tot <- sum(ss)
  out <- data.frame(term=names(ss), eta2=round(100*ss/tot,2), row.names=NULL)
  out[out$term != "Residuals", ]
}
VD <- NULL
for (m in ALLM) for (mt in names(METRICS)) {
  e <- eta_for(d[d$method==m, ], mt)
  if (is.null(e)) next
  e$method <- m; e$metric <- mt
  VD <- rbind(VD, e)
}
VD$variable <- VARS[VD$term]
W(VD, "variance_decomposition_long.csv")
# wide view for AUPRC (the headline)
vw <- VD[VD$metric=="ap", c("method","variable","eta2")]
vw <- reshape(vw, idvar="method", timevar="variable", direction="wide")
names(vw) <- gsub("eta2.","",names(vw), fixed=TRUE)
vw <- vw[match(FOCUS, vw$method), ]
W(vw, "variance_decomposition_auprc.csv")
print(vw, row.names=FALSE)

cat("\n=== 1b. same for high-PIP reliability (trustworthiness) ===\n")
vr <- VD[VD$metric=="hi_pip_reliab", c("method","variable","eta2")]
vr <- reshape(vr, idvar="method", timevar="variable", direction="wide")
names(vr) <- gsub("eta2.","",names(vr), fixed=TRUE)
vr <- vr[match(FOCUS, vr$method), ]
W(vr, "variance_decomposition_reliab.csv")
print(vr, row.names=FALSE)

# =============================================================================
# 2. MARGINAL PROFILES - median metric at each level of each variable, per method
# =============================================================================
cat("\n=== 2. marginal profiles written ===\n")
prof <- NULL
for (v in names(VARS)) for (mt in names(METRICS)) {
  a <- aggregate(d[[mt]], list(method=d$method, level=d[[v]]),
                 function(x) median(x, na.rm=TRUE))
  names(a)[3] <- "value"; a$variable <- VARS[[v]]; a$metric <- mt
  prof <- rbind(prof, a)
}
prof$value <- round(prof$value, 4)
W(prof, "marginal_profiles.csv")
cat("  rows:", nrow(prof), "\n")

# =============================================================================
# 3. RELATIVE STANDING per stratum: AP as a fraction of the best DEPLOYABLE
#    method in the same cell, plus win / near-win rates.
# =============================================================================
cat("\n=== 3. per-stratum standing ===\n")
key <- c("job_dir","S","phi","region_size")
dd <- d[d$method %in% DEPLOY, ]
cellbest <- aggregate(ap ~ job_dir + S + phi + region_size, dd,
                      function(x) max(x, na.rm=TRUE))
names(cellbest)[5] <- "best_ap"
dd <- merge(dd, cellbest, by=key)
dd$rel_ap  <- dd$ap / pmax(dd$best_ap, 1e-9)
dd$is_best <- abs(dd$ap - dd$best_ap) < 1e-12
dd$near    <- dd$rel_ap >= 0.95
# calibration gate: FDR violation <= 0.05 AND at least some confident calls
dd$trust   <- (!is.na(dd$max_fdr_violation_n20) & dd$max_fdr_violation_n20 <= 0.05)
dd$good    <- dd$near & dd$trust

standing <- aggregate(cbind(rel_ap, is_best, near, trust, good) ~ method + ann + ld, dd,
                      function(x) mean(x, na.rm=TRUE))
standing[,4:8] <- lapply(standing[,4:8], function(z) round(100*z,1))
standing$rel_ap <- round(standing$rel_ap/100, 3)   # rel_ap back to a fraction
names(standing) <- c("method","annotation","ld","rel_AP","pct_best","pct_within5","pct_FDRok","pct_both")
W(standing, "standing_by_stratum.csv")
cat("  written; example (binary, ref500):\n")
ex <- standing[standing$annotation=="binary" & standing$ld=="ref500", ]
print(ex[order(-ex$pct_both), ], row.names=FALSE)

# overall standing (for the per-method headline) - kept stratified then averaged
ov <- aggregate(cbind(rel_AP, pct_best, pct_within5, pct_FDRok, pct_both) ~ method,
                standing, function(x) round(mean(x),2))
ov <- ov[order(-ov$pct_both), ]
W(ov, "standing_overall.csv")
cat("\n  overall (mean over the 9 annotation x LD strata):\n"); print(ov, row.names=FALSE)

# =============================================================================
# 4. BEST / WORST conditions per method
# =============================================================================
cat("\n=== 4. best/worst strata per method ===\n")
bw <- NULL
for (m in DEPLOY) {
  z <- standing[standing$method==m, ]
  if (!nrow(z)) next
  z <- z[order(-z$rel_AP), ]
  bw <- rbind(bw, data.frame(method=m,
    best_stratum  = paste0(z$annotation[1], " / ", z$ld[1]),
    best_relAP    = z$rel_AP[1],
    worst_stratum = paste0(z$annotation[nrow(z)], " / ", z$ld[nrow(z)]),
    worst_relAP   = z$rel_AP[nrow(z)],
    spread        = round(z$rel_AP[1] - z$rel_AP[nrow(z)], 3)))
}
bw <- bw[order(-bw$best_relAP), ]
W(bw, "best_worst_strata.csv"); print(bw, row.names=FALSE)

# =============================================================================
# 5. ROBUSTNESS: degradation from in-sample to ref200
# =============================================================================
cat("\n=== 5. robustness to LD mis-specification ===\n")
rb <- aggregate(cbind(ap, hi_pip_reliab, total_mass_ratio) ~ method + ld, d,
                function(x) median(x, na.rm=TRUE))
r1 <- reshape(rb[,c("method","ld","ap")], idvar="method", timevar="ld", direction="wide")
names(r1) <- gsub("ap.","AP_",names(r1), fixed=TRUE)
r2 <- reshape(rb[,c("method","ld","hi_pip_reliab")], idvar="method", timevar="ld", direction="wide")
names(r2) <- gsub("hi_pip_reliab.","REL_",names(r2), fixed=TRUE)
rob <- merge(r1, r2, by="method")
rob$AP_retained  <- round(rob$`AP_ref200` / rob$`AP_in-sample`, 3)
rob$REL_retained <- round(rob$`REL_ref200`/ rob$`REL_in-sample`, 3)
num <- sapply(rob, is.numeric); rob[num] <- lapply(rob[num], function(z) round(z,3))
rob <- rob[order(-rob$AP_retained), ]
W(rob, "robustness.csv"); print(rob, row.names=FALSE)

# =============================================================================
# 6. INTERCHANGEABILITY: pairwise median |dAP| and % of cells within 0.01
# =============================================================================
cat("\n=== 6. interchangeability (pairwise) ===\n")
wide <- reshape(d[d$method %in% FOCUS, c(key,"method","ap")],
                idvar=key, timevar="method", direction="wide")
names(wide) <- gsub("ap.","",names(wide), fixed=TRUE)
ms <- intersect(FOCUS, names(wide))
pairs <- NULL
for (i in seq_along(ms)) for (j in seq_along(ms)) if (i < j) {
  a <- wide[[ms[i]]]; b <- wide[[ms[j]]]
  ok <- is.finite(a) & is.finite(b)
  pairs <- rbind(pairs, data.frame(m1=ms[i], m2=ms[j],
      med_abs_diff = round(median(abs(a[ok]-b[ok])), 4),
      pct_within_01 = round(100*mean(abs(a[ok]-b[ok]) < 0.01), 1)))
}
pairs <- pairs[order(-pairs$pct_within_01), ]
W(pairs, "interchangeability.csv")
print(utils::head(pairs, 12), row.names=FALSE)

cat("\nDONE - tables in", TAB, "\n"); print(list.files(TAB))
