#!/bin/bash
# =============================================================================
# scripts/analysis/submit_iter007_collect.sh
#
# >>> ITERATION 007 ONLY - TEMPORARY. See docs/autoresearch/iteration-007-REVERT.md
#
# Stage A only: re-collect Iteration 007's stored PIPs into L1 partials, one task
# per grid row, with scripts/analysis/iter004_collect.R. That collector reads
# only sim.rds (the truth) and results.rds (the PIPs), so the reference panel
# does not touch it.
#
# There is deliberately no Stage B. iter004_report.R and iter005_report.R both
# assume in-sample LD: they leave n_ref out of the aggregation keys, overwrite it
# with NA, and check that every row is in-sample. Iteration 007's strata are LD
# level x annotation type, so L1 is bound and analysed locally with n_ref kept.
#
#   SKIP_COUNT=1 bash scripts/analysis/submit_iter007_collect.sh
# =============================================================================
set -euo pipefail

PROJECT_ROOT="$(pwd)"
SCRATCH="${FMB_SCRATCH:-${EPHEMERAL:-/rds/general/user/$USER/ephemeral}/fmbench_iter007}"
BENCH_ROOT="${BENCH_ROOT:-${SCRATCH}/results/benchmark}"
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/results/iter007}"
L1_DIR="${L1_DIR:-${OUT_DIR}/L1}"
LOG_DIR="${LOG_DIR:-${OUT_DIR}/logs}"
R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
RSCRIPT="${RSCRIPT:-Rscript}"
# v1_small72a: v1_small24 rejects arrays (see submit_iter006_analysis.sh).
QUEUE="${QUEUE:-v1_small72a}"
SELECT="${SELECT:-1:ncpus=1:mem=16gb}"
WALLTIME="${WALLTIME:-06:00:00}"
EXP_ROWS=10; EXP_SCEN=600       # generate_params_grid_iter007.R: 10 rows x 60 scenarios

[[ -d "$BENCH_ROOT" ]] || { echo "ERROR: benchmark output not found at $BENCH_ROOT" >&2; exit 1; }
mkdir -p "$L1_DIR" "$LOG_DIR" "$OUT_DIR"

# The grid, written into OUT_DIR for the local analysis to join on job_dir. Not
# scripts/hpc/params_grid.csv, which holds whichever design was submitted last.
GRID_CSV="${OUT_DIR}/params_grid_iter007.csv"
"$RSCRIPT" "${PROJECT_ROOT}/scripts/hpc/generate_params_grid_iter007.R" "$GRID_CSV" >/dev/null
"$RSCRIPT" "${PROJECT_ROOT}/scripts/hpc/check_grid_columns.R" "$GRID_CSV" \
    "${PROJECT_ROOT}/scripts/hpc/generate_params_grid.R" \
  || { echo "ABORTING: grid is missing columns the pipeline reads." >&2; exit 1; }

N_ROWS=$(find "$BENCH_ROOT" -maxdepth 1 -type d -name 'job_*' | wc -l | tr -d ' ')
if (( N_ROWS != EXP_ROWS )); then
  echo "ERROR: expected $EXP_ROWS job_* rows under $BENCH_ROOT, found $N_ROWS" >&2; exit 1
fi
if [[ "${SKIP_COUNT:-0}" == "1" ]]; then
  # The count reads one directory per scenario. Skip it once check_iter007.sh
  # shows every fitting task finished; the Stage A logs report skipped scenarios.
  N_SCEN=$EXP_SCEN
else
  N_SCEN=$(find "$BENCH_ROOT" -mindepth 3 -maxdepth 3 -path '*/job_*/scenario_*/results.rds' | wc -l | tr -d ' ')
fi
echo "Benchmark root : $BENCH_ROOT"
echo "Grid           : $GRID_CSV"
echo "Rows           : $N_ROWS"
if [[ "${SKIP_COUNT:-0}" == "1" ]]; then echo "Scenarios      : count skipped (SKIP_COUNT=1)"
else echo "Scenarios      : $N_SCEN of $EXP_SCEN expected"; fi
if (( N_SCEN < EXP_SCEN )); then
  echo "ERROR: the run is incomplete - run scripts/hpc/check_iter007.sh first." >&2; exit 1
fi

A_SCRIPT="$(mktemp -t fmb7a_XXXXXX.sh)"
cat > "$A_SCRIPT" <<EOF
#!/bin/bash
#PBS -N fmb7collect
#PBS -q ${QUEUE}
#PBS -l select=${SELECT}
#PBS -l walltime=${WALLTIME}
#PBS -J 1-${N_ROWS}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
echo "[collect \${PBS_ARRAY_INDEX} on \$(hostname)] start \$(date)"
Rscript scripts/analysis/iter004_collect.R "\${PBS_ARRAY_INDEX}" "${BENCH_ROOT}" "${L1_DIR}" results
echo "[collect \${PBS_ARRAY_INDEX}] done \$(date)"
EOF
A_ID="$(qsub "$A_SCRIPT")"
echo
echo "Stage A (collect, ${N_ROWS} tasks): ${A_ID}"
cat <<EOF

Output : ${L1_DIR}/L1_<row>.rds, one per row, and ${GRID_CSV}
Check  : bash scripts/hpc/check_iter007.sh   (section 2, "collect")
EOF
