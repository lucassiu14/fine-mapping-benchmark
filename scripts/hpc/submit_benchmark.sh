#!/bin/bash
# =============================================================================
# scripts/hpc/submit_benchmark.sh
#
# Submit the fine-mapping benchmark as a SLURM job array.
#
# Usage (from project root):
#   bash scripts/hpc/submit_benchmark.sh
#
# The script generates the parameter grid (if not already present), then
# submits a 40-task array job.  Each task runs one (model, p, annotation)
# combination, sweeping all S × phi (× p_causal) values internally.
#
# Edit the SBATCH directives and configuration block below before submitting.
# =============================================================================

set -euo pipefail

# =============================================================================
# Configuration — edit to match your HPC environment
# =============================================================================

PROJECT_ROOT="$(pwd)"          # must be run from the project root
RSCRIPT="Rscript"              # or full path: /path/to/Rscript
LOG_DIR="${PROJECT_ROOT}/logs/benchmark"

# SLURM resource requests
#   Runtime: sparse jobs ~1-2 h, sparse_inf jobs ~4-6 h; 8 h is safe
#   Memory : ~8 GB per job (mostly Python/torch for BEATRICE)
#   CPUs   : 1 per job (all methods are single-threaded here)
#   Account/partition: edit to match your cluster
SLURM_TIME="08:00:00"
SLURM_MEM="8G"
SLURM_CPUS=1
SLURM_PARTITION="short"        # change to your partition name
SLURM_ACCOUNT=""               # set to your account/project code if required
SLURM_CONDA_ENV="beatrice"     # conda environment with Python + torch deps

# =============================================================================
# Generate parameter grid if not present
# =============================================================================

GRID_PATH="${PROJECT_ROOT}/scripts/hpc/params_grid.csv"

if [[ ! -f "$GRID_PATH" ]]; then
  echo "Generating parameter grid..."
  "$RSCRIPT" "${PROJECT_ROOT}/scripts/hpc/generate_params_grid.R"
fi

N_JOBS=$(tail -n +2 "$GRID_PATH" | wc -l | tr -d ' ')
echo "Parameter grid: ${N_JOBS} jobs"

# =============================================================================
# Create log directory
# =============================================================================

mkdir -p "$LOG_DIR"

# =============================================================================
# Build optional SLURM account flag
# =============================================================================

ACCOUNT_FLAG=""
if [[ -n "$SLURM_ACCOUNT" ]]; then
  ACCOUNT_FLAG="--account=${SLURM_ACCOUNT}"
fi

# =============================================================================
# Submit array job
# =============================================================================

JOB_SCRIPT=$(cat <<'SBATCH_SCRIPT'
#!/bin/bash
#SBATCH --job-name=fm_benchmark
#SBATCH --output=LOGDIR/job_%A_%a.out
#SBATCH --error=LOGDIR/job_%A_%a.err
#SBATCH --time=SLURM_TIME_PLACEHOLDER
#SBATCH --mem=SLURM_MEM_PLACEHOLDER
#SBATCH --cpus-per-task=SLURM_CPUS_PLACEHOLDER
#SBATCH --partition=SLURM_PARTITION_PLACEHOLDER
ACCOUNT_LINE

# Activate conda environment (adjust to your cluster's module system)
# Option 1: conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate CONDA_ENV_PLACEHOLDER

# Option 2: module system — uncomment if your HPC uses modules instead:
# module load miniconda3
# conda activate CONDA_ENV_PLACEHOLDER

# Move to project root
cd PROJECT_ROOT_PLACEHOLDER

# Run the worker script for this array task
Rscript scripts/hpc/run_benchmark_job.R ${SLURM_ARRAY_TASK_ID}
SBATCH_SCRIPT
)

# Substitute placeholders
ACCOUNT_LINE=""
if [[ -n "$SLURM_ACCOUNT" ]]; then
  ACCOUNT_LINE="#SBATCH --account=${SLURM_ACCOUNT}"
fi

JOB_SCRIPT="${JOB_SCRIPT/LOGDIR/$LOG_DIR}"
JOB_SCRIPT="${JOB_SCRIPT/SLURM_TIME_PLACEHOLDER/$SLURM_TIME}"
JOB_SCRIPT="${JOB_SCRIPT/SLURM_MEM_PLACEHOLDER/$SLURM_MEM}"
JOB_SCRIPT="${JOB_SCRIPT/SLURM_CPUS_PLACEHOLDER/$SLURM_CPUS}"
JOB_SCRIPT="${JOB_SCRIPT/SLURM_PARTITION_PLACEHOLDER/$SLURM_PARTITION}"
JOB_SCRIPT="${JOB_SCRIPT/ACCOUNT_LINE/$ACCOUNT_LINE}"
JOB_SCRIPT="${JOB_SCRIPT/CONDA_ENV_PLACEHOLDER/$SLURM_CONDA_ENV}"
JOB_SCRIPT="${JOB_SCRIPT/PROJECT_ROOT_PLACEHOLDER/$PROJECT_ROOT}"

# Write to temp file and submit
TMPSCRIPT=$(mktemp /tmp/fm_benchmark_XXXXXX.sh)
echo "$JOB_SCRIPT" > "$TMPSCRIPT"

echo ""
echo "Submitting array job (1-${N_JOBS})..."
sbatch --array="1-${N_JOBS}" "$TMPSCRIPT"
rm -f "$TMPSCRIPT"

echo ""
echo "Logs will appear in: $LOG_DIR"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  tail -f ${LOG_DIR}/job_<ARRAY_ID>_<TASK>.out"
echo ""
echo "Once all jobs complete, collect results with:"
echo "  Rscript scripts/hpc/collect_results.R"
