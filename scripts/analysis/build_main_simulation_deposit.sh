#!/bin/bash
# =============================================================================
# scripts/analysis/build_main_simulation_deposit.sh
#
# Stage the main simulation's results for deposit (Zenodo), under names a reader
# can navigate without the code.
#
# Keeps the posterior inclusion probabilities, the credible sets and the
# evaluation metrics. Drops the figures, the sensitivity analyses and the
# scheduler logs, which are either reproducible from what remains or specific to
# the cluster this ran on.
#
# Files are APFS clones where available and hard links otherwise, so staging
# 1.7 GB costs no additional disk. Nothing under results/iter004 is modified.
#
#   bash scripts/analysis/build_main_simulation_deposit.sh [SRC] [DEST]
#
# SRC defaults to results/iter004, DEST to zenodo/main-simulation.
# =============================================================================
set -euo pipefail

SRC="${1:-results/iter004}"
DEST="${2:-zenodo/main-simulation}"
[[ -d "$SRC" ]] || { echo "ERROR: no $SRC" >&2; exit 1; }
rm -rf "$DEST"
mkdir -p "$DEST"/{evaluation-metrics,credible-sets/per-fit,posterior-inclusion-probabilities,runtime}

# Clone if the filesystem supports it (APFS), else hard link, else copy.
place() {
  cp -c "$1" "$2" 2>/dev/null || ln "$1" "$2" 2>/dev/null || cp "$1" "$2"
}

echo "staging from $SRC -> $DEST"

# --- evaluation metrics -----------------------------------------------------
place "$SRC/combined_scenario_metrics_with_power.rds" "$DEST/evaluation-metrics/per-cell.rds"
place "$SRC/combined_replicate_metrics.rds"           "$DEST/evaluation-metrics/per-replicate.rds"
place "$SRC/combined_fit_metrics.rds"                 "$DEST/evaluation-metrics/per-fit.rds"
place "$SRC/fdr_thresholds.rds"      "$DEST/evaluation-metrics/false-discovery-rate-by-threshold.rds"
place "$SRC/calibration_bands9.rds"  "$DEST/evaluation-metrics/calibration-by-pip-band.rds"
place "$SRC/calib_fdr_periter.rds"   "$DEST/evaluation-metrics/calibration-and-fdr-per-iteration.rds"

# --- credible sets ----------------------------------------------------------
place "$SRC/native_credible_sets.rds" "$DEST/credible-sets/summary-by-method.rds"
for f in "$SRC"/credsets/cs_*.rds; do
  place "$f" "$DEST/credible-sets/per-fit/$(basename "$f" | sed 's/^cs_//')"
done

# --- posterior inclusion probabilities --------------------------------------
for f in "$SRC"/piptail/piptail_*.rds; do
  place "$f" "$DEST/posterior-inclusion-probabilities/$(basename "$f" | sed 's/^piptail_//')"
done

# --- runtime ----------------------------------------------------------------
# The aux store holds two record lists: `runtime`, which is a reported result,
# and `importance`, the LassoNet annotation weights, which are a model output
# rather than an evaluation metric and are not deposited. Written fresh rather
# than cloned, because the clone shares storage with results/iter004 and must
# not be edited in place.
Rscript -e '
  src <- sort(Sys.glob(file.path("'"$SRC"'", "aux", "aux_*.rds")))
  dst <- "'"$DEST"'/runtime"
  n <- 0L
  for (f in src) {
    a <- readRDS(f)
    a$importance <- NULL
    out <- file.path(dst, sub("^aux_", "", basename(f)))
    tmp <- paste0(out, ".tmp"); saveRDS(a, tmp); file.rename(tmp, out)
    n <- n + length(a$runtime)
  }
  cat("runtime:", length(src), "files,", format(n, big.mark = ","),
      "records, annotation importance dropped\n")
' 2>&1 | grep -v "out-of-sync\|renv::status"

# --- the documentation that makes the folder navigable ----------------------
# Kept in docs/deposit so it is version controlled and survives the rm -rf above.
for d in README.md DATA-DICTIONARY.md; do
  [[ -f "docs/deposit/$d" ]] || { echo "ERROR: docs/deposit/$d missing" >&2; exit 1; }
  cp "docs/deposit/$d" "$DEST/$d"
done

# --- the design, with the join key spelled out ------------------------------
Rscript -e '
  g <- read.csv("scripts/hpc/params_grid.csv", stringsAsFactors = FALSE)
  g$job_dir <- sprintf("job_%03d_%s", g$job_id, g$label)
  g <- g[, c("job_dir", setdiff(names(g), "job_dir"))]
  write.csv(g, file.path("'"$DEST"'", "simulation-grid.csv"), row.names = FALSE)
  cat("simulation-grid.csv:", nrow(g), "rows\n")
  # The per-cell table as CSV too, as the interoperable entry point.
  L3 <- readRDS(file.path("'"$DEST"'", "evaluation-metrics", "per-cell.rds"))
  f  <- file.path("'"$DEST"'", "evaluation-metrics", "per-cell.csv")
  write.csv(L3, f, row.names = FALSE)
  cat("per-cell.csv:", nrow(L3), "rows,", round(file.size(f)/1024^2), "MB\n")
' 2>&1 | grep -v "out-of-sync\|renv::status"

echo
echo "staged:"
du -sh "$DEST"/* | sed 's/^/  /'
echo "total: $(du -sh "$DEST" | cut -f1)"
echo
echo "file counts: pips $(ls "$DEST"/posterior-inclusion-probabilities | wc -l | tr -d ' '), credible sets $(ls "$DEST"/credible-sets/per-fit | wc -l | tr -d ' '), runtime $(ls "$DEST"/runtime | wc -l | tr -d ' ')"
