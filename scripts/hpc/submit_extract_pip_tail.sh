#!/bin/bash
# =============================================================================
# scripts/hpc/submit_extract_pip_tail.sh
#
# Rescue the PIP tail so threshold-based metrics can be recomputed at ANY
# threshold after $EPHEMERAL purges.
#
# The collected tables hold counts at five fixed thresholds only. A power-FDR
# curve, or a rate at any other threshold, needs the PIP vectors themselves,
# which exist nowhere but results.rds on ephemeral. Storing the full vectors is
# ~17 GB; storing only the tail above a floor is a few percent of that and
# reconstructs every threshold above the floor exactly. Causal variants are
# retained regardless of their PIP, so recall is never inflated by truncation.
#
# Third of three rescue jobs, alongside submit_extract_truth.sh (causal indices
# from sim.rds) and submit_extract_aux.sh (runtime and annotation importance).
#
#   BENCH_ROOT=$EPHEMERAL/fmbench_iter004/results/benchmark \
#   OUT_DIR=$HOME/fine-mapping-benchmark/results/iter004/piptail \
#   bash scripts/hpc/submit_extract_pip_tail.sh
# =============================================================================
set -euo pipefail

PROJECT_ROOT="$(pwd)"
BENCH_ROOT="${BENCH_ROOT:-${EPHEMERAL}/fmbench_iter004/results/benchmark}"
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/results/iter004/piptail}"
LOG_DIR="${LOG_DIR:-${OUT_DIR}/logs}"
FLOOR="${FLOOR:-0.01}"
R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
QUEUE="${QUEUE:-v1_small72a}"
SELECT="${SELECT:-1:ncpus=1:mem=32gb}"
WALLTIME="${WALLTIME:-06:00:00}"

[[ -d "$BENCH_ROOT" ]] || { echo "ERROR: no such tree: $BENCH_ROOT" >&2; exit 1; }
EXTRACT="${PROJECT_ROOT}/scripts/hpc/extract_pip_tail.R"
[[ -f "$EXTRACT" ]] || { echo "ERROR: extract_pip_tail.R missing" >&2; exit 1; }

N_ROWS=$(find "$BENCH_ROOT" -maxdepth 1 -type d -name 'job_*' | wc -l | tr -d ' ')
(( N_ROWS >= 1 )) || { echo "ERROR: no job_* directories" >&2; exit 1; }
N_SIM=$(find "$BENCH_ROOT" -maxdepth 2 -name sim.rds | wc -l | tr -d ' ')

mkdir -p "$OUT_DIR" "$LOG_DIR"
echo "Benchmark root : $BENCH_ROOT"
echo "Rows           : $N_ROWS   (sim.rds present: $N_SIM)"
echo "PIP floor      : $FLOOR"
echo "Output         : $OUT_DIR"
if (( N_SIM < N_ROWS )); then
  echo
  echo "WARNING: $(( N_ROWS - N_SIM )) row(s) have no sim.rds, so their tails"
  echo "         will be stored WITHOUT causal flags and cannot be scored."
  echo "         Recover truth for those rows first if it still exists."
fi

JOB="$(mktemp -t fmbtail_XXXXXX.sh)"
cat > "$JOB" <<EOF
#!/bin/bash
#PBS -N fmbtail
#PBS -q ${QUEUE}
#PBS -l select=${SELECT}
#PBS -l walltime=${WALLTIME}
#PBS -J 1-${N_ROWS}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
echo "[tail \${PBS_ARRAY_INDEX} on \$(hostname)] start \$(date)"
Rscript scripts/hpc/extract_pip_tail.R "\${PBS_ARRAY_INDEX}" "${BENCH_ROOT}" "${OUT_DIR}" "${FLOOR}"
echo "[tail \${PBS_ARRAY_INDEX}] done \$(date)"
EOF
echo
echo "Submitting array 1-${N_ROWS} to ${QUEUE} ..."
qsub "$JOB"
cat <<EOF

When it finishes:
  ls ${OUT_DIR}/piptail_*.rds | wc -l     # should be ${N_ROWS}
  du -sh ${OUT_DIR}

Each log line reports the percentage of PIPs retained. In Iteration 004 this
should be a few percent; if it is far higher the floor is too low for the
available space, and the job can be re-run at a larger FLOOR.
EOF
