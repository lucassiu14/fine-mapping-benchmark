#!/bin/bash
# =============================================================================
# scripts/hpc/check_lambda_sweep.sh
#
# One status report for the end of the lambda_l1 sweep, read from job logs and
# output folders only. It never walks the scenario tree, so it is fine on a
# login node.
#
#   1. the ten reruns of array 4005367's walltime-killed annotated tasks
#   2. the Iteration 004 full-PIP rescue (FLOOR=0 -> results/iter004_fullpips)
#   3. the sweep extraction: aux, pip tail, and the L1/L2/L3 collect
#
#   cd ~/fine-mapping-benchmark && bash scripts/hpc/check_lambda_sweep.sh
# =============================================================================
set -o pipefail
shopt -s nullglob

REPO="${REPO:-$PWD}"
BENCH_LOGS="${BENCH_LOGS:-${EPHEMERAL:-/rds/general/user/$USER/ephemeral}/fmbench_iter004/logs/benchmark}"
FULLPIPS="${FULLPIPS:-$REPO/results/iter004_fullpips}"
SWEEP="${SWEEP:-$REPO/results/fb_lambda_sweep}"
RERUNS="4028028 4028030 4028032 4028033 4028034 4028036 4028038 4028042 4028123 4028268"

# Expected values and their sources:
#   2031230  fits in results/iter004/piptail, made from the same results.rds at floor 0.01
#   11250    scenarios in the grid (45 rows x 250), all present when those tails were made
#   ~10500   scenarios with a sweep result: 40 annotated rows x 250, plus the 500
#            finished on the skipped no-annotation rows (4005367 logs)
#   ~840000  sweep fits: ~10500 scenarios x 8 arms x 10 regions
EXP_FULL_FITS=2031230
EXP_SCEN=11250

n_killed() { grep -l 'job killed\|Execution halted' "$@" 2>/dev/null | wc -l | tr -d ' '; }

echo "== 1. Reruns of the walltime-killed sweep tasks"
n_ok=0; n_seen=0
for j in $RERUNS; do
  files=( "$BENCH_LOGS"/$j\[*\].pbs-7.OU )
  if (( ${#files[@]} == 0 )); then echo "   $j: no log yet"; continue; fi
  for f in "${files[@]}"; do
    n_seen=$((n_seen + 1))
    if grep -q 'finished at' "$f"; then st=finished; n_ok=$((n_ok + 1)); else st=NOT-DONE; fi
    printf '   %-22s %-9s %s\n' "$(basename "$f")" "$st" \
      "$(grep -o '\[task complete\] [0-9]* scenario\|Nothing to do' "$f" | head -1)"
  done
done
echo "   -> $n_ok of $n_seen logs finished (expect 20 of 20)"

echo
echo "== 2. Iteration 004 full-PIP rescue"
if [[ ! -d "$FULLPIPS" ]]; then
  echo "   not run ($FULLPIPS absent)"
else
  files=( "$FULLPIPS"/piptail_*.rds ); ers=( "$FULLPIPS"/logs/*.ER )
  echo "   files: ${#files[@]} (expect 45)"
  if (( ${#ers[@]} )); then
    grep -h 'wrote piptail' "${ers[@]}" | awk -v e="$EXP_FULL_FITS" '
      {gsub(/\(/, "", $3); n++; f += $3; mb += $(NF-1)}
      END {printf "   rows written: %d, fits: %d (expect %d), size: %.1f GB\n", n, f, e, mb / 1024}'
    echo "   killed or failed tasks: $(n_killed "${ers[@]}")"
  fi
fi

echo
echo "== 3. Sweep extraction"
if [[ ! -d "$SWEEP" ]]; then
  echo "   not run ($SWEEP absent)"
else
  # aux: runtimes and annotation importances (results.rds with results_supp.rds overlaid)
  files=( "$SWEEP"/aux/aux_*.rds ); ers=( "$SWEEP"/aux/logs/*.ER )
  echo "   aux      files: ${#files[@]} (expect 45)"
  if (( ${#ers[@]} )); then
    t=$(grep -h 'wrote aux_' "${ers[@]}" | awk '{gsub(/\(/, "", $3); t += $3} END {print t + 0}')
    s=$(grep -h 'overlaid results_supp.rds in' "${ers[@]}" | awk '{s += $4} END {print s + 0}')
    echo "            scenarios read: $t (expect $EXP_SCEN), with a sweep result: $s (expect ~10500)"
    echo "            killed or failed tasks: $(n_killed "${ers[@]}")"
  fi

  # pip tail of the sweep's own results
  files=( "$SWEEP"/piptail/piptail_*.rds ); ers=( "$SWEEP"/piptail/logs/*.ER )
  echo "   piptail  files: ${#files[@]} (expect 45)"
  if (( ${#ers[@]} )); then
    grep -h 'floor = ' "${ers[@]}" | awk '{c["floor=" $6 " source=" $9]++}
      END {for (k in c) printf "            %d rows at %s  (expect 45 at floor=0 source=supp)\n", c[k], k}'
    grep -h 'wrote piptail' "${ers[@]}" | awk '
      {gsub(/\(/, "", $3); f += $3; mb += $(NF-1)}
      END {printf "            fits: %d (expect ~840000), size: %.1f GB\n", f, mb / 1024}'
    echo "            killed or failed tasks: $(n_killed "${ers[@]}")"
  fi

  # collect: Stage A array (L1 partials) and Stage B (L2/L3, validity, senses)
  l1=( "$SWEEP"/L1/L1_*.rds ); ersA=( "$SWEEP"/logs/*\[*\]*.ER ); ersB=()
  for f in "$SWEEP"/logs/*.ER; do [[ "$f" == *\[* ]] || ersB+=( "$f" ); done
  echo "   collect  L1 partials: ${#l1[@]} (expect 45)"
  if (( ${#ersA[@]} )); then
    grep -h 'wrote L1_' "${ersA[@]}" | awk '
      {f = $2; s = $6
       if (f ~ /anNone/) {nn++; ns += s}
       else {na++; as += s; if (s > 0) printf "            INCOMPLETE %s: %d scenarios missing\n", f, s}}
      END {printf "            annotated rows: %d collected, %d scenarios missing (expect 40, 0)\n", na, as
           printf "            no-annotation rows: %d collected, %d scenarios missing (skipped by design)\n", nn, ns}'
    echo "            killed or failed Stage A tasks: $(n_killed "${ersA[@]}")"
  fi
  if (( ${#ersB[@]} )); then
    # The sense scripts may exit non-zero on a sweep-only table (Sense G wants susie);
    # report.R only messages on that, so look for its own milestones instead.
    grep -h '^  L1: \|^  L3: \|^validity: \|^Done\. ' "${ersB[@]}" | sed 's/^ */            Stage B: /'
  fi
  for f in combined_fit_metrics.rds combined_replicate_metrics.rds combined_scenario_metrics.rds validity_checks.txt; do
    if [[ -e "$SWEEP/$f" ]]; then echo "   present  $f"; else echo "   MISSING  $f"; fi
  done
  [[ -e "$SWEEP/validity_checks.txt" ]] && sed 's/^/            /' "$SWEEP/validity_checks.txt"
fi
