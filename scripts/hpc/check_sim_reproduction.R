#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/check_sim_reproduction.R
#
# Canary for the comparator-fidelity re-run: does a fresh output root reproduce
# Iteration 004's simulations exactly?
#
# The re-run writes to its own output root (Iteration 004 is the standard and is
# never touched), so there is no cached sim.rds to inherit and every row is
# simulated again from seed = 1000 + row. That reproduces Iteration 004 only if
# the simulator's RNG stream and the VCF directory are unchanged since it ran.
# R/simulate_phenotypes.R was edited after Iteration 004 (the Iteration 006
# relationships), so this cannot be assumed - it has to be checked.
#
# The check: run an UNCHANGED, deterministic method alongside the three being
# re-run. susie is the obvious one. If the simulations are identical, susie's
# per-cell average precision must equal Iteration 004's to the last digit. If it
# does not, the re-run is not paired with Iteration 004 and the comparison the
# LSR makes would be across different data.
#
# Usage, after the canary rows have been collected into the new root:
#   FMB_OUTPUT_ROOT=$EPHEMERAL/fmbench/results/iter008_fidelity \
#     Rscript scripts/hpc/check_sim_reproduction.R
#
# Options:
#   FMB_REFERENCE   stored Iteration 004 cell table
#                   (default results/iter004/combined_scenario_metrics_with_power.rds)
#   FMB_CANARY_METHOD   method to compare (default susie)
#   FMB_CANARY_TOL      absolute AP tolerance (default 0, i.e. bit-identical)
# =============================================================================

OUTPUT_ROOT <- Sys.getenv("FMB_OUTPUT_ROOT", unset = "results/benchmark")
REFERENCE   <- Sys.getenv("FMB_REFERENCE",
                          unset = "results/iter004/combined_scenario_metrics_with_power.rds")
METHOD      <- Sys.getenv("FMB_CANARY_METHOD", unset = "susie")
TOL         <- as.numeric(Sys.getenv("FMB_CANARY_TOL", unset = "0"))

new_file <- file.path(OUTPUT_ROOT, "combined_scenario_metrics.rds")
if (!file.exists(new_file))
  stop("no collected table at ", new_file,
       "\n  Run scripts/hpc/collect_results.R against this output root first.",
       call. = FALSE)
if (!file.exists(REFERENCE))
  stop("no reference table at ", REFERENCE, call. = FALSE)

new <- readRDS(new_file)
ref <- readRDS(REFERENCE)

KEY <- c("job_dir", "S", "phi", "region_size")
stopifnot(all(c(KEY, "method", "ap") %in% names(new)),
          all(c(KEY, "method", "ap") %in% names(ref)))

new <- new[new$method == METHOD, c(KEY, "ap", "n_fits", "n_failed")]
ref <- ref[ref$method == METHOD, c(KEY, "ap", "n_fits", "n_failed")]
if (nrow(new) == 0L)
  stop("the re-run has no ", METHOD, " cells. Include it in FMB_METHODS - it is ",
       "the only thing that can show the simulations were reproduced.", call. = FALSE)

k <- function(d) do.call(paste, c(d[KEY], sep = "|"))
new$key <- k(new); ref$key <- k(ref)

cat(sprintf("Canary method : %s\n", METHOD))
cat(sprintf("Re-run cells  : %d (%s)\n", nrow(new), OUTPUT_ROOT))
cat(sprintf("Reference     : %d cells (%s)\n", nrow(ref), REFERENCE))

common <- intersect(new$key, ref$key)
cat(sprintf("Cells in both : %d\n\n", length(common)))
if (length(common) == 0L)
  stop("no cells in common - the job labels or the grid differ.", call. = FALSE)

a <- new[match(common, new$key), ]
b <- ref[match(common, ref$key), ]
d <- abs(a$ap - b$ap)

bad <- which(d > TOL | (a$n_fits != b$n_fits))
if (length(bad) == 0L) {
  cat(sprintf("PASS - every one of the %d shared cells matches to within %g.\n",
              length(common), TOL))
  cat("The re-run's simulations are the same data Iteration 004 used, so the\n")
  cat("re-run comparator results are exactly paired with the stored ones.\n")
  quit(status = 0)
}

cat(sprintf("FAIL - %d of %d cells differ (max |dAP| = %.3e).\n",
            length(bad), length(common), max(d)))
cat("The regenerated simulations are NOT Iteration 004's. Do not read the\n")
cat("re-run against the stored results until this is resolved; the likely\n")
cat("causes are a changed simulator RNG path or a changed VCF directory.\n\n")
show <- head(bad, 10)
print(data.frame(cell      = common[show],
                 ap_rerun  = round(a$ap[show], 6),
                 ap_iter004 = round(b$ap[show], 6),
                 diff      = signif(d[show], 3),
                 fits_rerun = a$n_fits[show],
                 fits_ref   = b$n_fits[show]),
      row.names = FALSE)
quit(status = 1)
