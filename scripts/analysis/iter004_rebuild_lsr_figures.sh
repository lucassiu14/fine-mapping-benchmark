#!/bin/bash
# =============================================================================
# scripts/analysis/iter004_rebuild_lsr_figures.sh
#
# Everything that has to happen after the pip-tail rescue lands, in order.
# Run from the repo root, locally, once results/iter004/piptail/ is populated.
#
#   bash scripts/analysis/iter004_rebuild_lsr_figures.sh [figure_out_dir]
#
# The validator gates the rebuild deliberately: a failed truth join inside
# extract_pip_tail.R produces tails in which nothing is causal, which would
# render as "every method is catastrophically overconfident" rather than as an
# error. Do not skip it.
# =============================================================================
set -euo pipefail

FIGDIR="${1:-${HOME}/Downloads/LSR_joint__1_/figures}"
RES="results/iter004"

command -v Rscript >/dev/null || { echo "Rscript not on PATH" >&2; exit 1; }
[[ -d "$RES/piptail" ]] || { echo "no $RES/piptail - rsync it down from the cluster first:
  rsync -av login-b.cx3.hpc.ic.ac.uk:fine-mapping-benchmark/$RES/piptail/ $RES/piptail/" >&2; exit 1; }

n=$(ls "$RES"/piptail/piptail_*.rds 2>/dev/null | wc -l | tr -d ' ')
echo "== piptail files: $n (expect 45)"
(( n == 45 )) || echo "   WARNING: incomplete; figures will be built from a subset"

echo "== 1/4  validating the rescue against the frozen band counts"
Rscript scripts/analysis/iter004_check_piptail.R "$RES/piptail" "$RES/combined_scenario_metrics.rds"

echo "== 2/4  attaching credible-set power to L3"
Rscript scripts/analysis/iter004_add_power.R "$RES"

echo "== 3/4  rebuilding calibration at the nine bands"
Rscript scripts/analysis/iter004_calibration_bands.R \
  "$RES/piptail" "$RES/combined_scenario_metrics_with_power.rds" "$RES/calibration_bands9.rds"

echo "== 4/4  regenerating the five figures into $FIGDIR"
mkdir -p "$FIGDIR"
Rscript scripts/analysis/iter004_lsr_figures.R "$FIGDIR" \
  "$RES/combined_scenario_metrics_with_power.rds" "$RES/calibration_bands9.rds"

echo
echo "Done. Recompile the LSR; only fig_calibration.pdf should differ from the"
echo "version now in the manuscript, and it should carry ten points per method"
echo "rather than four."
