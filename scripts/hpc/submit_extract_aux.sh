#!/bin/bash
# =============================================================================
# scripts/hpc/submit_extract_truth.sh
#
# Rescue the ground truth from a benchmark tree on $EPHEMERAL before it purges.
#
# $EPHEMERAL is deleted 30 days after write and is not backed up. sim.rds is the
# only copy of the causal variant indices, and it is the OLDEST object in the
# tree: written once per row at the start of the run and never rewritten, even
# when results are recomputed. It therefore expires FIRST, and once it goes the
# surviving results.rds cannot be scored against anything.
#
# One array task per grid row. Each reads that row's sim.rds (large - it carries
# the genotype matrices) and writes a few-KB truth table, so the output is a
# rounding error next to the input.
#
#   BENCH_ROOT=$EPHEMERAL/fmbench_iter004/results/benchmark \
#   OUT_DIR=$HOME/fine-mapping-benchmark/results/iter004/aux \
#   bash scripts/hpc/submit_extract_truth.sh
# =============================================================================
set -euo pipefail

PROJECT_ROOT="$(pwd)"
BENCH_ROOT="${BENCH_ROOT:-${EPHEMERAL}/fmbench_iter004/results/benchmark}"
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/results/iter004/aux}"
LOG_DIR="${LOG_DIR:-${OUT_DIR}/logs}"
R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
QUEUE="${QUEUE:-v1_small72a}"
SELECT="${SELECT:-1:ncpus=1:mem=16gb}"
WALLTIME="${WALLTIME:-04:00:00}"

[[ -d "$BENCH_ROOT" ]] || { echo "ERROR: no such tree: $BENCH_ROOT" >&2; exit 1; }
EXTRACT="${PROJECT_ROOT}/scripts/hpc/extract_aux.R"
[[ -f "$EXTRACT" ]] || { echo "ERROR: extract_aux.R missing" >&2; exit 1; }

N_ROWS=$(find "$BENCH_ROOT" -maxdepth 1 -type d -name 'job_*' | wc -l | tr -d ' ')
(( N_ROWS >= 1 )) || { echo "ERROR: no job_* directories under $BENCH_ROOT" >&2; exit 1; }

# How many rows still HAVE a sim.rds. If this is below N_ROWS, truth has already
# been lost for the difference and no job can recover it.
N_SIM=$(find "$BENCH_ROOT" -maxdepth 2 -name sim.rds | wc -l | tr -d ' ')

mkdir -p "$OUT_DIR" "$LOG_DIR"
echo "Benchmark root : $BENCH_ROOT"
echo "Rows           : $N_ROWS   (sim.rds still present: $N_SIM)"
echo "Output         : $OUT_DIR"
if (( N_SIM < N_ROWS )); then
  echo
  echo "WARNING: $(( N_ROWS - N_SIM )) row(s) have already lost their sim.rds."
  echo "         Those rows are unrecoverable; the rest will still be extracted."
fi

JOB="$(mktemp -t fmbaux_XXXXXX.sh)"
cat > "$JOB" <<EOF
#!/bin/bash
#PBS -N fmbaux
#PBS -q ${QUEUE}
#PBS -l select=${SELECT}
#PBS -l walltime=${WALLTIME}
#PBS -J 1-${N_ROWS}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
echo "[truth \${PBS_ARRAY_INDEX} on \$(hostname)] start \$(date)"
Rscript scripts/hpc/extract_aux.R "\${PBS_ARRAY_INDEX}" "${BENCH_ROOT}" "${OUT_DIR}"
echo "[truth \${PBS_ARRAY_INDEX}] done \$(date)"
EOF
echo
echo "Submitting array 1-${N_ROWS} to ${QUEUE} ..."
qsub "$JOB"
cat <<EOF

When it finishes, check every row was rescued:
  ls ${OUT_DIR}/aux_*.rds | wc -l      # should be ${N_SIM}

The output is small and lives in project space, which is not purged. It is
enough to re-score any PIP file that still exists, and to recompute every
truth-dependent metric in the analysis.
EOF
