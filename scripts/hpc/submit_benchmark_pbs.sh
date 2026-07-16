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
# 32gb: each task handles ONE scenario but ALL 20 regions (the region
# panel must stay together for cross-region pooling in scenario_setup).
# It holds the 20-region sim (~1 GB in memory) + one scenario's results
# (small) + a torch subprocess for BEATRICE/FB. The earlier 8gb OOM was
# the OLD full-row task accumulating all 125 scenarios; a single-scenario
# task needs far less. 32gb is a safe margin that still packs ~4 tasks
# onto a 128gb node - important for parallelism across 3125 tasks.
PBS_SELECT="${PBS_SELECT:-1:ncpus=1:mem=32gb}"
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
# The array is subdivided by scenario: one task per (grid row, scenario),
# with all regions kept together in each task. Total tasks =
# N_ROWS * SCENARIOS_PER_ROW, where SCENARIOS_PER_ROW = |S| * |phi| * n_iter.
# Columns are located by header name so column order can change.
N_ROWS=$(( $(wc -l < "$PARAMS_CSV") - 1 ))
if (( N_ROWS < 1 )); then
  echo "ERROR: params_grid.csv has no rows" >&2
  exit 1
fi
_col() { head -1 "$PARAMS_CSV" | tr ',' '\n' | grep -n "\"$1\"" | cut -d: -f1; }
S_COL=$(_col S_values); PHI_COL=$(_col phi_values); NITER_COL=$(_col n_iter)
if [[ -z "$S_COL" || -z "$PHI_COL" || -z "$NITER_COL" ]]; then
  echo "ERROR: params_grid.csv missing S_values/phi_values/n_iter columns" >&2
  exit 1
fi
SCENARIOS_PER_ROW=$(awk -F, -v s="$S_COL" -v ph="$PHI_COL" -v ni="$NITER_COL" '
  NR==2 {
    gsub(/"/,"",$s); gsub(/"/,"",$ph); gsub(/"/,"",$ni);
    ns = split($s, a, "[|]"); np = split($ph, b, "[|]");
    print ns * np * ($ni + 0);
  }' "$PARAMS_CSV")
N_JOBS=$(( N_ROWS * SCENARIOS_PER_ROW ))
echo "Grid: ${N_ROWS} rows x ${SCENARIOS_PER_ROW} scenarios = ${N_JOBS} array tasks"
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
# CRITICAL: export the scratch output root INTO the compute-node job.
# The login-side export does not propagate to the node, so without this
# the worker falls back to results/ under home and overflows the quota.
# \${OUTPUT_ROOT} is expanded here at submit time to the literal path.
export FMB_OUTPUT_ROOT="${OUTPUT_ROOT}"
# Supplemental re-run selector, forwarded the same way (a login-side export
# does NOT reach the node). Empty string = normal full-method run.
export FMB_METHODS="${FMB_METHODS:-}"
module load ${R_MODULE}
module load ${GSL_MODULE}
module load ${PYTHON_MODULE}
source ${PY_VENV_ACTIVATE}
echo "[node:\$(hostname)] task \${PBS_ARRAY_INDEX} of ${ARRAY_RANGE} starting at \$(date)"
echo "[node:\$(hostname)] FMB_OUTPUT_ROOT=\${FMB_OUTPUT_ROOT}"
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
