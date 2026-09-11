#!/bin/bash
# =============================================================================
# scripts/hpc/submit_extract_aux.sh
#
# Rescue RUNTIME and ANNOTATION-IMPORTANCE from a benchmark tree on $EPHEMERAL.
#
# iter004_collect.R keeps only the accuracy metrics, so two quantities the paper
# needs survive nowhere but results.rds on ephemeral: per-fit runtime, and the
# LassoNet annotation importances returned by the Functional BEATRICE family.
# Neither can be recomputed without re-running the benchmark.
#
# Companion to submit_extract_truth.sh, which rescues the causal indices from
# sim.rds. Run both: they read different files and neither subsumes the other.
#
#   BENCH_ROOT=$EPHEMERAL/fmbench_iter004/results/benchmark \
#   OUT_DIR=$HOME/fine-mapping-benchmark/results/iter004/aux \
#   bash scripts/hpc/submit_extract_aux.sh
#
# extract_aux.R overlays results_supp.rds onto results.rds per method, so after
# a supplemental re-run the extract is no longer Iteration 004's: point OUT_DIR
# somewhere new. The script refuses an OUT_DIR that already holds aux_*.rds
# (OVERWRITE=1 to allow). SKIP_COUNT=1 skips the results.rds count, which reads
# one directory per scenario.
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

# How many scenario results survive. results.rds is rewritten whenever a row is
# recomputed, so it is younger than sim.rds and expires later - but not by much.
# The count reads one directory per scenario (~11,000 on RDS) - not for a login
# node. SKIP_COUNT=1 skips it; extract_aux.R logs scenarios read per row.
if [[ "${SKIP_COUNT:-0}" == "1" ]]; then
  N_RES="count skipped"
else
  N_RES=$(find "$BENCH_ROOT" -maxdepth 3 -name results.rds | wc -l | tr -d ' ')
fi

# An existing extract is not overwritten by accident: after a supplemental
# re-run the overlay would silently replace the re-run methods' records.
if compgen -G "${OUT_DIR}/aux_*.rds" > /dev/null && [[ "${OVERWRITE:-0}" != "1" ]]; then
  echo "ERROR: $OUT_DIR already holds aux_*.rds. Choose a new OUT_DIR, or set OVERWRITE=1." >&2
  exit 1
fi

mkdir -p "$OUT_DIR" "$LOG_DIR"
echo "Benchmark root : $BENCH_ROOT"
echo "Rows           : $N_ROWS   (results.rds present: $N_RES)"
echo "Output         : $OUT_DIR"
if [[ "$N_RES" == "0" ]]; then
  echo "ERROR: no results.rds found under $BENCH_ROOT - nothing to extract." >&2
  exit 1
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
  ls ${OUT_DIR}/aux_*.rds | wc -l      # should be ${N_ROWS}

The output lives in project space, which is not purged, and carries the
per-fit runtimes and annotation importances that the collected metric tables
do not retain.
EOF
