#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter002_macro_auprc.R
#
# Corrects the AUPRC aggregation from MICRO (pooled) to MACRO (per-replicate
# average), and writes combined_scenario_metrics_macroAP.{rds,csv}.
#
# WHY. combined_scenario_metrics.rds computes ap by POOLING every fit in a
# scenario cell (10 iterations x 2 same-length regions) into one ranked list and
# taking a single AUPRC. AUPRC is a ranking metric: the standard, more defensible
# estimator of expected per-dataset performance is to compute AUPRC PER REPLICATE
# and average (macro), which also yields a standard error. Pooling instead
# estimates "AUPRC of all loci merged into one ranked list" - not a real
# experiment - and systematically UNDER-states AUPRC (a locus's false positives
# outrank another locus's true positives), most at low S. Measured on Iteration
# 002: micro < macro by ~0.05 at S=1, ~0.015 at S=10.
#
# WHAT changes / stays. Only AUPRC is affected. FDR-violation and PIP-calibration
# are RATE metrics estimated from counts; per-replicate they are too sparse to
# estimate (a single locus rarely has enough PIP>=0.9 calls or high-threshold
# selections), so pooling is correct for them and they are carried through
# unchanged. (They are also not reconstructable locally - combined_evaluation.rds
# stores only ranking scalars per task, not the curves.)
#
# GRANULARITY. The replicate unit here is the ITERATION (10 per cell); each
# iteration's AUPRC itself pools its 2 same-length regions, because that is the
# finest granularity stored in combined_evaluation.rds (the by_region_size
# stratum). A future collect run could add a per-region stratum for a true
# per-locus macro; the dominant correction - not pooling across iterations - is
# captured here.
#
# Usage: Rscript scripts/analysis/iter002_macro_auprc.R [eval_rds] [metrics_rds] [out_rds]
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
EV  <- if (length(args) >= 1) args[1] else "results/iter002_fixed/combined_evaluation.rds"
MET <- if (length(args) >= 2) args[2] else "results/iter002_fixed/combined_scenario_metrics.rds"
OUT <- if (length(args) >= 3) args[3] else "results/iter002_fixed/combined_scenario_metrics_macroAP.rds"

ev <- readRDS(EV); sc <- readRDS(MET)

# Per-iteration AUPRC: each `scenario` (scenario_XXX) is one (S, phi, iteration);
# its by_region_size rows give that iteration's AUPRC per region length.
rs <- ev[ev$stratum == "by_region_size", c("job_dir", "scenario", "method", "stratum_value", "ap")]
names(rs)[4] <- "region_size"
# S / phi for each (job, scenario) come from its single-valued by_S / by_phi strata.
Smap <- unique(ev[ev$stratum == "by_S",   c("job_dir", "scenario", "stratum_value")]); names(Smap)[3] <- "S"
Pmap <- unique(ev[ev$stratum == "by_phi", c("job_dir", "scenario", "stratum_value")]); names(Pmap)[3] <- "phi"
rs <- merge(merge(rs, Smap, by = c("job_dir", "scenario")), Pmap, by = c("job_dir", "scenario"))

# Macro AUPRC = mean of the per-iteration AUPRCs; SE across iterations.
macro <- aggregate(ap ~ job_dir + S + phi + region_size + method, data = rs,
                   FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                       se   = stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))),
                                       n    = sum(!is.na(x))))
macro <- do.call(data.frame, macro)
names(macro)[6:8] <- c("ap_macro", "ap_se", "ap_n")

key <- c("job_dir", "S", "phi", "region_size", "method")
for (k in c("S", "phi", "region_size")) { sc[[k]] <- as.character(sc[[k]]); macro[[k]] <- as.character(macro[[k]]) }

out <- merge(sc, macro[, c(key, "ap_macro", "ap_se", "ap_n")], by = key, all.x = TRUE)
out$ap_micro <- out$ap          # keep the pooled value for reference
out$ap <- out$ap_macro          # `ap` now = the corrected macro AUPRC

saveRDS(out, OUT)
write.csv(out, sub("\\.rds$", ".csv", OUT), row.names = FALSE)

cat(sprintf("wrote %s\n  matched %d / %d scenario cells; iterations/cell: %s\n",
            OUT, sum(!is.na(out$ap_macro)), nrow(out),
            paste(range(out$ap_n, na.rm = TRUE), collapse = "-")))
d <- out$ap_micro - out$ap_macro
cat(sprintf("  micro - macro AUPRC: mean %.4f  (pooling under-states AUPRC); |max| %.4f\n",
            mean(d, na.rm = TRUE), max(abs(d), na.rm = TRUE)))
