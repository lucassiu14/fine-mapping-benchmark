#!/bin/bash
# =============================================================================
# scripts/hpc/submit_extract_credible_sets.sh
#
# Rescue the methods' OWN credible sets before $EPHEMERAL purges results.rds.
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
#   OUT_DIR=$HOME/fine-mapping-benchmark/results/iter004/credsets \
#   bash scripts/hpc/submit_extract_credible_sets.sh
# =============================================================================
set -euo pipefail

PROJECT_ROOT="$(pwd)"
BENCH_ROOT="${BENCH_ROOT:-${EPHEMERAL}/fmbench_iter004/results/benchmark}"
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/results/iter004/credsets}"
LOG_DIR="${LOG_DIR:-${OUT_DIR}/logs}"
# results | supp | overlay - see extract_credible_sets.R. Expanded at submit time
# into the job script, so no #PBS -V is needed for it to reach the node. The
# guard mirrors submit_extract_pip_tail.sh: a non-default source must not write
# into Iteration 004's credsets, which the LSR was submitted with.
SOURCE="${SOURCE:-results}"
if [[ "$SOURCE" != "results" && "$OUT_DIR" == "${PROJECT_ROOT}/results/iter004/credsets" ]]; then
  echo "ERROR: SOURCE=$SOURCE would write into $OUT_DIR, Iteration 004's credible sets." >&2
  echo "       Set a separate OUT_DIR, e.g. OUT_DIR=\$PWD/results/iter009_fb_xregion/credsets" >&2
  exit 1
fi
R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
QUEUE="${QUEUE:-v1_small72a}"
SELECT="${SELECT:-1:ncpus=1:mem=32gb}"
WALLTIME="${WALLTIME:-06:00:00}"

[[ -d "$BENCH_ROOT" ]] || { echo "ERROR: no such tree: $BENCH_ROOT" >&2; exit 1; }
EXTRACT="${PROJECT_ROOT}/scripts/hpc/extract_credible_sets.R"
[[ -f "$EXTRACT" ]] || { echo "ERROR: extract_credible_sets.R missing" >&2; exit 1; }

N_ROWS=$(find "$BENCH_ROOT" -maxdepth 1 -type d -name 'job_*' | wc -l | tr -d ' ')
(( N_ROWS >= 1 )) || { echo "ERROR: no job_* directories" >&2; exit 1; }
mkdir -p "$OUT_DIR" "$LOG_DIR"
echo "Benchmark root : $BENCH_ROOT"
echo "Rows           : $N_ROWS"
echo "Output         : $OUT_DIR"
echo "Source         : $SOURCE"

JOB="$(mktemp -t fmbcs_XXXXXX.sh)"
cat > "$JOB" <<EOF
#!/bin/bash
#PBS -N fmbcs
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
Rscript scripts/hpc/extract_credible_sets.R "\${PBS_ARRAY_INDEX}" "${BENCH_ROOT}" "${OUT_DIR}" "${SOURCE}"
echo "[tail \${PBS_ARRAY_INDEX}] done \$(date)"
EOF
echo
echo "Submitting array 1-${N_ROWS} to ${QUEUE} ..."
qsub "$JOB"
cat <<EOF

When it finishes:
  ls ${OUT_DIR}/cs_*.rds | wc -l     # should be ${N_ROWS}
  du -sh ${OUT_DIR}

Each log line reports how many fits carried at least one reported set.
EOF
