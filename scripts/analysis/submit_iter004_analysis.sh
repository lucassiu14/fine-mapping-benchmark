#!/bin/bash
# =============================================================================
# scripts/analysis/submit_iter004_analysis.sh
#
# ONE COMMAND: from the benchmark's stored PIPs to the finished analysis.
#
# Submits two chained PBS jobs.
#
#   STAGE A  array 1-45, one task per params_grid row. Re-collects from
#            results.rds + sim.rds into the L1 fit-level table. No refitting -
#            the expensive part is already done. Parallel because sim.rds is
#            cached per row and holds the truth for that row's 250 scenarios.
#
#   STAGE B  single job, -W depend=afterok on A. Binds the L1 partials, builds
#            L2 and L3, runs the §11 validity checks, then drives Sense V, D
#            and G.
#
# afterok means Stage B only runs if EVERY Stage A task succeeded. A partial
# collect would silently analyse a subset of the grid, which is worse than not
# running at all.
#
# Usage (from the project root):
#   bash scripts/analysis/submit_iter004_analysis.sh
#   BENCH_ROOT=/path/to/results/benchmark bash scripts/analysis/submit_iter004_analysis.sh
#
# Run it only once the fitting array has finished. Check first:
#   Rscript scripts/hpc/qc_run.R $BENCH_ROOT --sample 300
# =============================================================================
set -euo pipefail

PROJECT_ROOT="$(pwd)"
SCRATCH="${FMB_SCRATCH:-${EPHEMERAL:-/rds/general/user/$USER/ephemeral}/fmbench_iter004}"
BENCH_ROOT="${BENCH_ROOT:-${SCRATCH}/results/benchmark}"

# Analysis output goes to PROJECT SPACE, not ephemeral. The L1 table and the
# decomposition are small and are the thing we actually want to keep; ephemeral
# auto-purges and would take the raw PIPs with it.
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/results/iter004}"
L1_DIR="${L1_DIR:-${OUT_DIR}/L1}"
LOG_DIR="${LOG_DIR:-${OUT_DIR}/logs}"
GRID_CSV="${GRID_CSV:-${PROJECT_ROOT}/scripts/hpc/params_grid.csv}"

R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
QUEUE_A="${QUEUE_A:-v1_small24}"
QUEUE_B="${QUEUE_B:-v1_small24}"
# Stage A: 250 scenarios x ~190 fits per row, each an O(p log p) sort at p <= 2000.
SELECT_A="${SELECT_A:-1:ncpus=1:mem=32gb}"
WALL_A="${WALL_A:-08:00:00}"
# Stage B holds the whole L1 table (~2.1M rows) plus the aov per stratum x method
# x response. Memory-hungry, not CPU-hungry.
SELECT_B="${SELECT_B:-1:ncpus=1:mem=96gb}"
WALL_B="${WALL_B:-12:00:00}"

if [[ ! -d "$BENCH_ROOT" ]]; then
  echo "ERROR: benchmark output not found at $BENCH_ROOT" >&2
  echo "       Set BENCH_ROOT or FMB_SCRATCH." >&2
  exit 1
fi
if [[ ! -f "$GRID_CSV" ]]; then
  echo "ERROR: $GRID_CSV not found. Regenerate it with" >&2
  echo "       Rscript scripts/hpc/generate_params_grid.R $GRID_CSV" >&2
  exit 1
fi

N_ROWS=$(find "$BENCH_ROOT" -maxdepth 1 -type d -name 'job_*' | wc -l | tr -d ' ')
if (( N_ROWS < 1 )); then
  echo "ERROR: no job_* directories under $BENCH_ROOT" >&2
  exit 1
fi

# Refuse to analyse an unfinished run by accident. A partial grid is not a
# smaller version of the answer - the design is balanced by construction and the
# §7 closed-form correction assumes it.
N_SCEN=$(find "$BENCH_ROOT" -name 'results.rds' | wc -l | tr -d ' ')
EXPECTED=$(( N_ROWS * 250 ))
echo "Benchmark root : $BENCH_ROOT"
echo "Rows           : $N_ROWS"
echo "Scenarios      : $N_SCEN of $EXPECTED expected"
if (( N_SCEN < EXPECTED )); then
  echo
  echo "WARNING: the run is INCOMPLETE ($(( 100 * N_SCEN / EXPECTED ))%)."
  echo "         The design is balanced by construction and the variance"
  echo "         decomposition's closed-form noise correction assumes 20 fits in"
  echo "         every cell. Analysing a partial grid gives unbalanced strata and"
  echo "         a residual line the correction cannot interpret."
  echo
  if [[ "${ALLOW_PARTIAL:-0}" != "1" ]]; then
    echo "         Re-run with ALLOW_PARTIAL=1 to proceed anyway." >&2
    exit 1
  fi
  echo "         ALLOW_PARTIAL=1 set - proceeding."
fi

mkdir -p "$L1_DIR" "$LOG_DIR" "$OUT_DIR"

# --- Stage A -----------------------------------------------------------------
A_SCRIPT="$(mktemp -t fmb4a_XXXXXX.sh)"
cat > "$A_SCRIPT" <<EOF
#!/bin/bash
#PBS -N fmb4collect
#PBS -q ${QUEUE_A}
#PBS -l select=${SELECT_A}
#PBS -l walltime=${WALL_A}
#PBS -J 1-${N_ROWS}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
echo "[collect \${PBS_ARRAY_INDEX} on \$(hostname)] start \$(date)"
Rscript scripts/analysis/iter004_collect.R \\
    "\${PBS_ARRAY_INDEX}" "${BENCH_ROOT}" "${L1_DIR}"
echo "[collect \${PBS_ARRAY_INDEX}] done \$(date)"
EOF
A_ID="$(qsub "$A_SCRIPT")"
echo
echo "Stage A (collect, ${N_ROWS} tasks): ${A_ID}"

# --- Stage B -----------------------------------------------------------------
B_SCRIPT="$(mktemp -t fmb4b_XXXXXX.sh)"
cat > "$B_SCRIPT" <<EOF
#!/bin/bash
#PBS -N fmb4analyse
#PBS -q ${QUEUE_B}
#PBS -l select=${SELECT_B}
#PBS -l walltime=${WALL_B}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
echo "[analyse on \$(hostname)] start \$(date)"
Rscript scripts/analysis/iter004_report.R "${L1_DIR}" "${GRID_CSV}" "${OUT_DIR}"
echo "[analyse] done \$(date)"
EOF
B_ID="$(qsub -W depend=afterok:"${A_ID}" "$B_SCRIPT")"
echo "Stage B (analyse, held on afterok): ${B_ID}"

cat <<EOF

Output    : ${OUT_DIR}
  combined_fit_metrics.rds        L1 - the pivot; every statistic re-derivable from it
  combined_replicate_metrics.rds  L2 - paired differences for the Sense D gate
  combined_scenario_metrics.rds   L3 - cell means, the ANOVA input
  validity_checks.txt             §11 - read this FIRST
  iter004_sense_v/sense_v.txt     which variables move the metric
  iter004_sense_d/sense_d_ap.txt  which variables change the winner
  iter004_sense_g/sense_g_ap.txt  which variables hurt THIS method specifically

Track     : qstat -tan -u \$USER | grep -E 'fmb4collect|fmb4analyse'

Stage B is held until every Stage A task succeeds (afterok). If A partly fails,
B stays held: fix the cause, resubmit A for the failed indices, then release B
with  qrls ${B_ID}  -- or qdel it and re-run this script, which is safe because
the L1 partials are per-row and simply overwritten.
EOF
