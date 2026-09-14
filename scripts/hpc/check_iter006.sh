#!/bin/bash
# =============================================================================
# scripts/hpc/check_iter006.sh
#
# >>> ITERATION 006 ONLY - TEMPORARY. See docs/autoresearch/iteration-006-REVERT.md
#
# One status report for Iteration 006, read from job logs and output folders
# only. It never walks the scenario tree, so it is fine on a login node.
#
#   1. the fitting array: 160 tasks of 5 scenarios, 9 methods
#   2. the extraction into results/iter006: aux (runtimes, importances), the
#      full PIPs, and the L1/L2/L3 collect with its validity checks
#
#   cd ~/fine-mapping-benchmark && bash scripts/hpc/check_iter006.sh
# =============================================================================
set -o pipefail
shopt -s nullglob

REPO="${REPO:-$PWD}"
SCRATCH="${FMB_SCRATCH:-${EPHEMERAL:-/rds/general/user/$USER/ephemeral}/fmbench_iter006}"
BENCH="${BENCH_ROOT:-$SCRATCH/results/benchmark}"
LOGS="${BENCH_LOGS:-$SCRATCH/logs/benchmark}"
OUT="${OUT_DIR:-$REPO/results/iter006}"
ARRAY="${ARRAY:-4031997}"

# Expected values. From generate_params_grid_iter006.R: 8 rows x 100 scenarios.
# From the submission: 5 scenarios per task, so 160 tasks; FMB_ITER006_METHODS
# selects 9 methods, each fitting 10 regions per scenario, so 8,000 fits per
# method and 72,000 in all.
EXP_TASKS=160; EXP_PER_TASK=5; EXP_ROWS=8; EXP_SCEN=800; EXP_FITS_METHOD=8000; EXP_FITS=72000

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
  echo "   method set (expect $EXP_TASKS tasks with 9 methods):"
  grep -h 'reduced method set' "${ou[@]}" | sort | uniq -c | sed 's/^ */      /'
  echo "   fits and failures per method (expect $EXP_FITS_METHOD fits each):"
  grep -h 'n_fits=.*failed=' "${ou[@]}" | awk '
    {sub("n_fits=", "", $2); sub("failed=", "", $3); n[$1] += $2; f[$1] += $3}
    END {for (m in n) printf "      %-22s fits=%d failed=%d\n", m, n[m], f[m]}' | sort
  echo "   relationships simulated (each row logs it once per task that built the row):"
  grep -h '\[iter005\] relationship=' "${ou[@]}" | awk '{sub("relationship=", "", $2); print $2}' | sort | uniq -c | sed 's/^ */      /'
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
    t=$(grep -h 'wrote aux_' "${ers[@]}" | awk '{gsub(/\(/, "", $3); t += $3} END {print t + 0}')
    echo "            scenarios read: $t (expect $EXP_SCEN)"
    echo "            killed or failed tasks: $(n_killed "${ers[@]}")"
  fi

  files=( "$OUT"/piptail/piptail_*.rds ); ers=( "$OUT"/piptail/logs/*.ER )
  echo "   piptail  files: ${#files[@]} (expect $EXP_ROWS)"
  if (( ${#ers[@]} )); then
    grep -h 'floor = ' "${ers[@]}" | awk -v e="$EXP_ROWS" '{c["floor=" $6 " source=" $9]++}
      END {for (k in c) printf "            %d rows at %s  (expect %d at floor=0 source=results)\n", c[k], k, e}'
    grep -h 'wrote piptail' "${ers[@]}" | awk -v e="$EXP_FITS" '{gsub(/\(/, "", $3); f += $3; mb += $(NF-1)}
      END {printf "            fits: %d (expect %d), size: %.1f GB\n", f, e, mb / 1024}'
    echo "            killed or failed tasks: $(n_killed "${ers[@]}")"
  fi

  l1=( "$OUT"/L1/L1_*.rds ); ersA=( "$OUT"/logs/*\[*\]*.ER ); ersB=()
  for f in "$OUT"/logs/*.ER; do [[ "$f" == *\[* ]] || ersB+=( "$f" ); done
  echo "   collect  L1 partials: ${#l1[@]} (expect $EXP_ROWS)"
  if (( ${#ersA[@]} )); then
    grep -h 'wrote L1_' "${ersA[@]}" | awk -v e="$EXP_ROWS" '
      {n++; s += $6; if ($6 > 0) printf "            INCOMPLETE %s: %d scenarios missing\n", $2, $6}
      END {printf "            rows collected: %d, scenarios missing: %d (expect %d, 0)\n", n, s, e}'
    echo "            killed or failed Stage A tasks: $(n_killed "${ersA[@]}")"
  fi
  if (( ${#ersB[@]} )); then
    grep -h '^  L1: \|^  L3: \|^validity: \|^Done\. ' "${ersB[@]}" | sed 's/^ */            Stage B: /'
  fi
  for f in combined_fit_metrics.rds combined_replicate_metrics.rds combined_scenario_metrics.rds validity_checks.txt; do
    if [[ -e "$OUT/$f" ]]; then echo "   present  $f"; else echo "   MISSING  $f"; fi
  done
  if [[ -e "$OUT/validity_checks.txt" ]]; then sed 's/^/            /' "$OUT/validity_checks.txt"; fi
fi
