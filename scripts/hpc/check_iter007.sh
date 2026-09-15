#!/bin/bash
# =============================================================================
# scripts/hpc/check_iter007.sh
#
# >>> ITERATION 007 ONLY - TEMPORARY. See docs/autoresearch/iteration-007-REVERT.md
#
# One status report for Iteration 007, read from job logs and output folders
# only. It never walks the scenario tree, so it is fine on a login node.
#
#   1. the fitting array: 120 tasks of 5 scenarios, 12 methods, five LD levels
#   2. the extraction into results/iter007: aux (runtimes, importances), the
#      full PIPs, and the L1 collect
#
#   cd ~/fine-mapping-benchmark && bash scripts/hpc/check_iter007.sh
# =============================================================================
set -o pipefail
shopt -s nullglob

REPO="${REPO:-$PWD}"
SCRATCH="${FMB_SCRATCH:-${EPHEMERAL:-/rds/general/user/$USER/ephemeral}/fmbench_iter007}"
BENCH="${BENCH_ROOT:-$SCRATCH/results/benchmark}"
LOGS="${BENCH_LOGS:-$SCRATCH/logs/benchmark}"
OUT="${OUT_DIR:-$REPO/results/iter007}"
ARRAY="${ARRAY:-4045723}"

# Expected values. From generate_params_grid_iter007.R: 10 rows (5 LD levels x 2
# annotation types) x 60 scenarios. From the submission: 5 scenarios per task, so
# 120 tasks and 24 per LD level. FMB_ITER007_METHODS selects 12 methods, each
# fitting 10 regions per scenario, so 6,000 fits per method and 72,000 in all.
EXP_TASKS=120; EXP_PER_TASK=5; EXP_ROWS=10; EXP_SCEN=600; EXP_METHODS=12
EXP_FITS_METHOD=6000; EXP_FITS=72000; EXP_TASKS_PER_LD=24

n_killed() { grep -l 'job killed\|Execution halted' "$@" 2>/dev/null | wc -l | tr -d ' '; }

echo "== 1. Fitting array $ARRAY"
ou=( "$LOGS"/$ARRAY\[*\].pbs-7.OU )
er=( "$LOGS"/$ARRAY\[*\].pbs-7.ER )
echo "   task logs: ${#ou[@]} (expect $EXP_TASKS)"
if (( ${#ou[@]} )); then
  nf=$(grep -L 'finished at' "${ou[@]}" | wc -l | tr -d ' ')
  echo "   tasks that did not finish cleanly: $nf (expect 0)"
  if (( nf )); then grep -L 'finished at' "${ou[@]}" | xargs -n1 basename | head -20 | sed 's/^/      /'; fi
  echo "   killed or R error: $(n_killed "${er[@]}") (expect 0)"
  echo "   scenarios completed per task - count of tasks, scenarios (expect $EXP_TASKS tasks with $EXP_PER_TASK):"
  grep -h '^\[task complete\]' "${ou[@]}" | awk '{print $3}' | sort | uniq -c | sed 's/^ */      /'
  echo "   method set (expect $EXP_TASKS tasks with $EXP_METHODS methods):"
  grep -h 'reduced method set' "${ou[@]}" | sort | uniq -c | sed 's/^ */      /'
  echo "   fits and failures per method (expect $EXP_FITS_METHOD fits each):"
  grep -h 'n_fits=.*failed=' "${ou[@]}" | awk '
    {sub("n_fits=", "", $2); sub("failed=", "", $3); n[$1] += $2; f[$1] += $3}
    END {for (m in n) printf "      %-22s fits=%d failed=%d\n", m, n[m], f[m]}' | sort
  echo "   LD in each task's banner (expect $EXP_TASKS_PER_LD tasks at each of the five levels):"
  grep -ho 'in-sample (perfect)\|reference panel n_ref=[0-9]*' "${ou[@]}" | sort | uniq -c | sed 's/^ */      /'
fi
sims=( "$BENCH"/job_*/sim.rds )
echo "   row simulations (sim.rds): ${#sims[@]} (expect $EXP_ROWS)"

echo
echo "== 2. Extraction into $OUT"
if [[ ! -d "$OUT" ]]; then
  echo "   not run ($OUT absent)"
else
  files=( "$OUT"/aux/aux_*.rds ); ers=( "$OUT"/aux/logs/*.ER )
  echo "   aux      files: ${#files[@]} (expect $EXP_ROWS)"
  if (( ${#ers[@]} )); then
    t=$(grep -h 'wrote aux_' "${ers[@]}" | awk '{gsub(/\(/, "", $3); v[$2] = $3} END {for (k in v) t += v[k]; print t + 0}')
    echo "            scenarios read: $t (expect $EXP_SCEN)"
    echo "            killed or failed tasks: $(n_killed "${ers[@]}")"
  fi

  files=( "$OUT"/piptail/piptail_*.rds ); ers=( "$OUT"/piptail/logs/*.ER )
  echo "   piptail  files: ${#files[@]} (expect $EXP_ROWS)"
  if (( ${#ers[@]} )); then
    grep -h 'floor = ' "${ers[@]}" | awk -v e="$EXP_ROWS" '{row[$3] = "floor=" $6 " source=" $9}
      END {for (r in row) c[row[r]]++; for (k in c) printf "            %d rows at %s  (expect %d at floor=0 source=results)\n", c[k], k, e}'
    grep -h 'wrote piptail' "${ers[@]}" | awk -v e="$EXP_FITS" '{gsub(/\(/, "", $3); f[$2] = $3; mb[$2] = $(NF-1)}
      END {for (k in f) {F += f[k]; MB += mb[k]}; printf "            fits: %d (expect %d), size: %.1f GB\n", F, e, MB / 1024}'
    echo "            killed or failed tasks: $(n_killed "${ers[@]}")"
  fi

  l1=( "$OUT"/L1/L1_*.rds ); ers=( "$OUT"/logs/*.ER )
  echo "   collect  L1 partials: ${#l1[@]} (expect $EXP_ROWS)"
  if (( ${#ers[@]} )); then
    grep -h 'wrote L1_' "${ers[@]}" | awk -v e="$EXP_ROWS" -v ef="$EXP_FITS" '{gsub(/\(/, "", $3); n[$2] = $3; m[$2] = $6}
      END {for (k in n) {r++; F += n[k]; s += m[k]; if (m[k] > 0) printf "            INCOMPLETE %s: %d scenarios skipped\n", k, m[k]}
           printf "            rows collected: %d, fit rows: %d, scenarios skipped: %d (expect %d, %d, 0)\n", r, F, s, e, ef}'
    echo "            rows reporting non-finite PIPs: $(grep -h 'NON-FINITE PIPs' "${ers[@]}" | wc -l | tr -d ' ') (expect 0)"
    echo "            killed or failed tasks: $(n_killed "${ers[@]}")"
  fi
fi
