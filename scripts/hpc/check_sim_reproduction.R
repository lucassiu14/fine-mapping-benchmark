#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/check_sim_reproduction.R
#
# Canary for a comparator re-run: did a fresh output root reproduce Iteration
# 004's simulations exactly?
#
# A re-run that must leave Iteration 004 untouched writes to its own output root,
# so there is no cached sim.rds to inherit and every row is simulated again from
# seed = 1000 + row. That reproduces Iteration 004 only if the simulator's RNG
# stream and the VCF directory are unchanged since it ran. R/simulate_phenotypes.R
# was edited after Iteration 004 (the Iteration 006 relationships), so this has to
# be checked rather than assumed.
#
# COMPARE AT FIT LEVEL, NOT CELL LEVEL. The first version of this script compared
# scripts/hpc/collect_results.R's cell table against Iteration 004's stored one and
# reported differences that had nothing to do with the simulations. The two tables
# are built by different pipelines and their `ap` columns are different estimators:
# collect_results.R reads AP off the 0.005 FDR grid and ranks a scenario's two
# same-size regions POOLED, whereas Iteration 004's numbers come from the per-fit
# recomputation in scripts/analysis/iter004_collect.R (its header sets out why).
# Their units differ too - ap_n counts replicates, n_fits counts fits. Nothing can
# be concluded from comparing them.
#
# So this compares the L1 fit-level tables, which are the same object on both
# sides, on two things:
#
#   n_variants   the realised region size after MAF filtering. It comes from the
#                genotype simulation alone - no method, no metric - so it is the
#                cleanest fingerprint of the simulation there is. If the VCF draw
#                or the RNG moved, this moves.
#   ap           per-fit average precision for an UNCHANGED, deterministic method
#                (susie). It depends on the genotypes, the annotations, the causal
#                set and the phenotype, so it catches anything n_variants misses.
#
# Build the new side on the cluster first, one task per row (no refitting - it
# reads the stored PIPs):
#
#   Rscript scripts/analysis/iter004_collect.R 1 \
#     $EPHEMERAL/fmbench/results/iter008_fidelity /path/to/L1 supp
#
# then run this against the directory of L1_*.rds files it wrote:
#
#   Rscript scripts/hpc/check_sim_reproduction.R /path/to/L1
#
# Options:
#   argument 2        reference L1 (default results/iter004/combined_fit_metrics.rds)
#   FMB_CANARY_METHOD method to compare (default susie)
#   FMB_CANARY_TOL    absolute AP tolerance (default 0, i.e. bit-identical)
# =============================================================================

args   <- commandArgs(trailingOnly = TRUE)
NEW    <- if (length(args) >= 1) args[1] else
  stop("usage: check_sim_reproduction.R <dir of L1_*.rds> [reference L1.rds]",
       call. = FALSE)
REF    <- if (length(args) >= 2) args[2] else
  "results/iter004/combined_fit_metrics.rds"
METHOD <- Sys.getenv("FMB_CANARY_METHOD", unset = "susie")
TOL    <- as.numeric(Sys.getenv("FMB_CANARY_TOL", unset = "0"))

if (!file.exists(REF)) stop("no reference L1 at ", REF, call. = FALSE)

new_files <- if (dir.exists(NEW)) Sys.glob(file.path(NEW, "L1_*.rds")) else NEW
if (length(new_files) == 0L)
  stop("no L1_*.rds under ", NEW,
       "\n  Build them with scripts/analysis/iter004_collect.R <row> <root> <out> supp",
       call. = FALSE)

new <- do.call(rbind, lapply(new_files, readRDS))
ref <- readRDS(REF)

KEY <- c("job_dir", "scenario_id", "region_id")
need <- c(KEY, "method", "ap", "n_variants")
for (nm in c("new", "ref")) {
  d <- get(nm)
  miss <- setdiff(need, names(d))
  if (length(miss))
    stop(nm, " L1 is missing: ", paste(miss, collapse = ", "), call. = FALSE)
}

new <- new[new$method == METHOD, ]
ref <- ref[ref$method == METHOD, ]
if (nrow(new) == 0L)
  stop("the re-run has no ", METHOD, " fits. Include it in FMB_METHODS - it is ",
       "the only thing that can show the simulations were reproduced.", call. = FALSE)

k <- function(d) do.call(paste, c(d[KEY], sep = "|"))
new$key <- k(new); ref$key <- k(ref)
common <- intersect(new$key, ref$key)

cat(sprintf("Canary method : %s\n", METHOD))
cat(sprintf("Re-run fits   : %d, over %d row(s)\n", nrow(new),
            length(unique(new$job_dir))))
cat(sprintf("Reference     : %s\n", REF))
cat(sprintf("Fits in both  : %d\n\n", length(common)))
if (length(common) == 0L)
  stop("no fits in common - the job labels or the scenario numbering differ.",
       call. = FALSE)

a <- new[match(common, new$key), ]
b <- ref[match(common, ref$key), ]

# 1. the simulation itself
vd <- a$n_variants != b$n_variants
cat(sprintf("n_variants  : %d of %d fits differ\n", sum(vd), length(common)))

# 2. the deterministic method on top of it
finite <- is.finite(a$ap) & is.finite(b$ap)
ad <- abs(a$ap - b$ap)
ad[!finite] <- NA_real_
bad_ap <- which(!is.na(ad) & ad > TOL)
cat(sprintf("%-11s : %d of %d comparable fits differ by more than %g (max %.3e)\n",
            paste0("ap(", METHOD, ")"), length(bad_ap), sum(finite), TOL,
            if (any(finite)) max(ad, na.rm = TRUE) else 0))
if (sum(!finite)) cat(sprintf("              (%d fits non-finite on one side)\n",
                              sum(!finite)))

if (sum(vd) == 0L && length(bad_ap) == 0L) {
  cat("\nPASS - the regenerated simulations are Iteration 004's, fit for fit.\n")
  cat("The re-run comparator results are exactly paired with the stored ones.\n")
  quit(status = 0)
}

cat("\nFAIL - the regenerated simulations are NOT Iteration 004's.\n")
if (sum(vd) > 0L) {
  cat("  n_variants differs, so the GENOTYPES differ: the VCF draw or the\n",
      "  simulator's RNG stream moved. Nothing downstream can be paired.\n", sep = "")
} else {
  cat("  Region sizes match, so the genotypes are the same, but a deterministic\n",
      "  method scores them differently: the phenotypes, causal sets or\n",
      "  annotations moved, or the method itself changed.\n", sep = "")
}
show <- head(order(-ad, na.last = NA), 10)
print(data.frame(job_dir  = a$job_dir[show],
                 scenario = a$scenario_id[show],
                 region   = a$region_id[show],
                 nvar_new = a$n_variants[show],
                 nvar_ref = b$n_variants[show],
                 ap_new   = round(a$ap[show], 6),
                 ap_ref   = round(b$ap[show], 6)),
      row.names = FALSE)
quit(status = 1)
