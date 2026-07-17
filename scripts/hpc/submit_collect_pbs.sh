#!/bin/bash
# scripts/hpc/submit_collect_pbs.sh
# Submit collect_results.R as a single PBS job. Collect reads ~2 x 3125
# evaluation files (tens of GB of I/O) and pools the FDR/calibration curves,
# which is far too heavy for a login node.
#
# Usage (from project root):
#   bash scripts/hpc/submit_collect_pbs.sh

set -euo pipefail

PROJECT_ROOT="$(pwd)"
RSCRIPT="${RSCRIPT:-Rscript}"

FMB_SCRATCH="${FMB_SCRATCH:-${EPHEMERAL:-/rds/general/user/$USER/ephemeral}/fmbench}"
OUTPUT_ROOT="${FMB_OUTPUT_ROOT:-${FMB_SCRATCH}/results/benchmark}"
LOG_DIR="${FMB_LOG_DIR:-${FMB_SCRATCH}/logs/collect}"

PBS_QUEUE="${PBS_QUEUE:-v1_small24}"
PBS_WALLTIME="${PBS_WALLTIME:-04:00:00}"
# Memory: the pooled curve accumulators + the ~300k-row scalar frame are
# modest, but rbind of many frames peaks. 64gb is comfortable headroom.
PBS_SELECT="${PBS_SELECT:-1:ncpus=1:mem=64gb}"

R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"

mkdir -p "$LOG_DIR"
echo "Output root: $OUTPUT_ROOT"
echo "Log dir:     $LOG_DIR"

JOB_SCRIPT="$(mktemp -t fmbcollect_XXXXXX.sh)"
cat > "$JOB_SCRIPT" <<PBS_EOF
#!/bin/bash
#PBS -N fmbcollect
#PBS -q ${PBS_QUEUE}
#PBS -l select=${PBS_SELECT}
#PBS -l walltime=${PBS_WALLTIME}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/

set -euo pipefail
cd "${PROJECT_ROOT}"
export FMB_OUTPUT_ROOT="${OUTPUT_ROOT}"
module load ${R_MODULE}
echo "[node:\$(hostname)] collect starting at \$(date)"
echo "[node:\$(hostname)] FMB_OUTPUT_ROOT=\${FMB_OUTPUT_ROOT}"
${RSCRIPT} scripts/hpc/collect_results.R
echo "[node:\$(hostname)] collect finished at \$(date)"
PBS_EOF

echo "Job script: $JOB_SCRIPT"
qsub "$JOB_SCRIPT"
echo
echo "Track with:  qstat -an | grep \$USER"
echo "Log:         ls -lh ${LOG_DIR}"
