#!/bin/bash
# =============================================================================
# scripts/hpc/submit_benchmark.sh
#
# Submit the fine-mapping benchmark as a SLURM job array.
#
# Usage (from project root, after `renv::restore()` and installing the package):
#   bash scripts/hpc/submit_benchmark.sh
#
# Generates the parameter grid if it doesn't exist, then submits one task per
# row. Each task is independent; failures don't take down the rest of the array.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration - edit to match your cluster
# -----------------------------------------------------------------------------
PROJECT_ROOT="$(pwd)"
RSCRIPT="${RSCRIPT:-Rscript}"
LOG_DIR="${PROJECT_ROOT}/logs/benchmark"

# SLURM resource requests (defaults: edit for your jobs and partition)
SLURM_TIME="06:00:00"        # wall-clock limit per task
SLURM_MEM="8G"               # memory per task
SLURM_CPUS=1                 # CPUs per task (methods are single-threaded)
SLURM_PARTITION=""           # e.g. "short", "compute" - SET THIS
SLURM_ACCOUNT=""             # set if your cluster requires it

# Optional: load a conda env that has Python deps for Tier 3 methods. Leave
# empty to skip. The env activation is sourced before each Rscript call.
SLURM_CONDA_ENV=""

# -----------------------------------------------------------------------------
# Generate the parameter grid if it doesn't already exist
# -----------------------------------------------------------------------------
PARAMS_CSV="${PROJECT_ROOT}/scripts/hpc/params_grid.csv"
if [[ ! -f "$PARAMS_CSV" ]]; then
  echo "Generating params_grid.csv ..."
  "$RSCRIPT" "${PROJECT_ROOT}/scripts/hpc/generate_params_grid.R" "$PARAMS_CSV"
fi
N_JOBS=$(( $(wc -l < "$PARAMS_CSV") - 1 ))   # minus header row
if (( N_JOBS < 1 )); then
  echo "ERROR: params_grid.csv has no jobs" >&2
  exit 1
fi
echo "Submitting array of ${N_JOBS} jobs."

# -----------------------------------------------------------------------------
# Prepare log directory
# -----------------------------------------------------------------------------
mkdir -p "$LOG_DIR"

# -----------------------------------------------------------------------------
# Build SBATCH directives
# -----------------------------------------------------------------------------
SBATCH_ARGS=(
  --job-name=fmbench
  --array=1-${N_JOBS}
  --time=${SLURM_TIME}
  --mem=${SLURM_MEM}
  --cpus-per-task=${SLURM_CPUS}
  --output=${LOG_DIR}/slurm-%A_%a.out
  --error=${LOG_DIR}/slurm-%A_%a.err
)
if [[ -n "${SLURM_PARTITION}" ]]; then
  SBATCH_ARGS+=(--partition="${SLURM_PARTITION}")
fi
if [[ -n "${SLURM_ACCOUNT}" ]]; then
  SBATCH_ARGS+=(--account="${SLURM_ACCOUNT}")
fi

# -----------------------------------------------------------------------------
# Submit
# -----------------------------------------------------------------------------
sbatch "${SBATCH_ARGS[@]}" --wrap "
  set -e
  cd '${PROJECT_ROOT}'
  if [[ -n '${SLURM_CONDA_ENV}' ]]; then
    # shellcheck disable=SC1091
    source \"\$(conda info --base)/etc/profile.d/conda.sh\"
    conda activate '${SLURM_CONDA_ENV}'
  fi
  '${RSCRIPT}' scripts/hpc/run_benchmark_job.R \"\${SLURM_ARRAY_TASK_ID}\"
"

echo "Submitted. Logs will appear under: ${LOG_DIR}"
echo "After the array finishes, run:"
echo "  ${RSCRIPT} scripts/hpc/collect_results.R"
