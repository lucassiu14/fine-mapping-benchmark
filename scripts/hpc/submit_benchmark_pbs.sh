#!/bin/bash
# scripts/hpc/submit_benchmark_pbs.sh
# PBS Pro driver for the fine-mapping benchmark array (Imperial HPC).
# SLURM users: see scripts/hpc/submit_benchmark.sh.

set -euo pipefail

PROJECT_ROOT="$(pwd)"
RSCRIPT="${RSCRIPT:-Rscript}"
PARAMS_CSV="${PROJECT_ROOT}/scripts/hpc/params_grid.csv"

# Output + logs go on scratch, NOT home. Each task writes ~100 MB+ of RDS
# files; a full array overflows a 1 TB home quota and saveRDS() dies with
# "error writing to connection". Default to personal ephemeral (multi-TB).
# Override with FMB_SCRATCH, or FMB_OUTPUT_ROOT / FMB_LOG_DIR individually.
FMB_SCRATCH="${FMB_SCRATCH:-${EPHEMERAL:-/rds/general/user/$USER/ephemeral}/fmbench}"
OUTPUT_ROOT="${FMB_OUTPUT_ROOT:-${FMB_SCRATCH}/results/benchmark}"
LOG_DIR="${FMB_LOG_DIR:-${FMB_SCRATCH}/logs/benchmark}"
export FMB_OUTPUT_ROOT="${OUTPUT_ROOT}"

# Config (override with env vars: PBS_QUEUE=... bash submit_benchmark_pbs.sh)
PBS_QUEUE="${PBS_QUEUE:-v1_small72a}"
PBS_WALLTIME="${PBS_WALLTIME:-72:00:00}"       # queue max; BEATRICE-family
                                                # is compute-heavy so start
                                                # generous, cut later if fine.
# 64gb: the 8gb default OOM-killed tasks. Each task holds the full
# 20-region sim (p up to 1000, so ~1000x1000 LD matrices), accumulates
# all 14 methods' per-fit results across 2500 fits in memory, AND spawns
# a torch subprocess for BEATRICE/FB - the R process + Python child both
# count against the cgroup limit. 64gb is comfortably under the
# v1_small72a 128gb node cap.
PBS_SELECT="${PBS_SELECT:-1:ncpus=1:mem=64gb}"
ARRAY_RANGE="${ARRAY_RANGE:-}"                  # e.g. "1-2" canary, "" full

R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
GSL_MODULE="${GSL_MODULE:-GSL/2.8-GCC-14.3.0}"
PYTHON_MODULE="${PYTHON_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
PY_VENV_ACTIVATE="${PY_VENV_ACTIVATE:-$HOME/tools/fmpy-venv/bin/activate}"

# Generate grid if missing
if [[ ! -f "$PARAMS_CSV" ]]; then
  echo "Generating params_grid.csv ..."
  "$RSCRIPT" "${PROJECT_ROOT}/scripts/hpc/generate_params_grid.R" "$PARAMS_CSV"
fi
N_JOBS=$(( $(wc -l < "$PARAMS_CSV") - 1 ))
if (( N_JOBS < 1 )); then
  echo "ERROR: params_grid.csv has no jobs" >&2
  exit 1
fi
if [[ -z "$ARRAY_RANGE" ]]; then ARRAY_RANGE="1-${N_JOBS}"; fi

mkdir -p "$LOG_DIR" "$OUTPUT_ROOT"
echo "Output root: $OUTPUT_ROOT"
echo "Log dir:     $LOG_DIR"

JOB_SCRIPT="$(mktemp -t fmbench_pbs_XXXXXX.sh)"
cat > "$JOB_SCRIPT" <<PBS_EOF
#!/bin/bash
#PBS -N fmbench
#PBS -q ${PBS_QUEUE}
#PBS -l select=${PBS_SELECT}
#PBS -l walltime=${PBS_WALLTIME}
#PBS -J ${ARRAY_RANGE}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/

set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
module load ${GSL_MODULE}
module load ${PYTHON_MODULE}
source ${PY_VENV_ACTIVATE}
echo "[node:\$(hostname)] task \${PBS_ARRAY_INDEX} of ${ARRAY_RANGE} starting at \$(date)"
${RSCRIPT} scripts/hpc/run_benchmark_job.R "\${PBS_ARRAY_INDEX}"
echo "[node:\$(hostname)] task \${PBS_ARRAY_INDEX} finished at \$(date)"
PBS_EOF

echo "Job script: $JOB_SCRIPT"
echo "Submitting array ${ARRAY_RANGE} to queue ${PBS_QUEUE} ..."
qsub "$JOB_SCRIPT"

echo
echo "Track with:"
echo "  qstat -tan \$USER          # per-array-element status"
echo "  ls -lh ${LOG_DIR}          # log files as they land"
