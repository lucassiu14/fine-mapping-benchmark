#!/bin/bash
#PBS -N fmb4reset
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=72:00:00
#PBS -j oe
# =============================================================================
# scripts/hpc/reset_and_resubmit.sh
#
# Invalidate a ROW RANGE's checkpoints and resubmit just those rows - as a JOB,
# because none of this belongs on a login node. Deleting checkpoints across 36
# job directories and counting 11,250 files over RDS takes minutes there and
# competes with everyone else's shell.
#
#   qsub -v PROJECT=$HOME/fine-mapping-benchmark scripts/hpc/reset_and_resubmit.sh
#
# Override the range or chunk with -v:
#   qsub -v PROJECT=...,FIRST_ROW=10,LAST_ROW=45,CHUNK=5 scripts/hpc/reset_and_resubmit.sh
#
# Set DRY=1 to report what WOULD be deleted and submitted, without doing either.
#
# WHY A ROW RANGE. Resume keys on evaluation.rds, so the worker redoes exactly
# the scenarios whose checkpoint is gone. Rows 1-9 (the sparse arm) were fitted
# against a healthy toolchain and must NOT be touched; rows 10-45 (sparse_inf)
# were fitted against a decayed one and are void.
# =============================================================================
set -uo pipefail

PROJECT="${PROJECT:-$HOME/fine-mapping-benchmark}"
cd "$PROJECT" || { echo "cannot cd to $PROJECT"; exit 1; }
SCRATCH="${FMB_SCRATCH:-${EPHEMERAL}/fmbench_iter004}"
BENCH="${SCRATCH}/results/benchmark"
FIRST_ROW="${FIRST_ROW:-10}"
LAST_ROW="${LAST_ROW:-45}"
CHUNK="${CHUNK:-5}"
SCEN_PER_ROW="${SCEN_PER_ROW:-250}"

# Task range derived, never hardcoded. run_benchmark_job.R maps
#   row = ((task-1) / TASKS_PER_ROW) + 1
# so getting this wrong silently refits the WRONG rows - including the clean
# sparse arm we are trying to preserve.
TASKS_PER_ROW=$(( (SCEN_PER_ROW + CHUNK - 1) / CHUNK ))
FIRST_TASK=$(( (FIRST_ROW - 1) * TASKS_PER_ROW + 1 ))
LAST_TASK=$((  LAST_ROW * TASKS_PER_ROW ))

echo "=== plan ==="
echo "rows            : $FIRST_ROW-$LAST_ROW"
echo "scenarios/row   : $SCEN_PER_ROW   chunk: $CHUNK   -> $TASKS_PER_ROW tasks/row"
echo "array range     : ${FIRST_TASK}-${LAST_TASK}"
echo "benchmark root  : $BENCH"

if [[ ! -d "$BENCH" ]]; then echo "ERROR: $BENCH not found"; exit 1; fi

echo
echo "=== invalidating checkpoints ==="
n_del=0
for (( i = FIRST_ROW; i <= LAST_ROW; i++ )); do
  d=$(printf '%s/job_%03d_*' "$BENCH" "$i")
  for jd in $d; do
    [[ -d "$jd" ]] || continue
    c=$(find "$jd" -name 'evaluation.rds' | wc -l | tr -d ' ')
    n_del=$(( n_del + c ))
    if [[ "${DRY:-0}" != "1" ]]; then
      find "$jd" -name 'evaluation.rds' -delete
    fi
  done
done
echo "evaluation.rds removed from rows $FIRST_ROW-$LAST_ROW : $n_del"

remaining=$(find "$BENCH" -name 'evaluation.rds' | wc -l | tr -d ' ')
echo "evaluation.rds remaining (rows 1-$((FIRST_ROW-1)))    : $remaining"

if [[ "${DRY:-0}" == "1" ]]; then
  echo
  echo "DRY=1 - nothing deleted, nothing submitted."
  exit 0
fi

echo
echo "=== submitting rows $FIRST_ROW-$LAST_ROW ==="
if ARRAY_RANGE="${FIRST_TASK}-${LAST_TASK}" FMB_SCENARIOS_PER_TASK="$CHUNK" \
   PBS_QUEUE="${PBS_QUEUE:-v1_small72a}" PBS_WALLTIME="${PBS_WALLTIME:-72:00:00}" \
   FMB_SCRATCH="$SCRATCH" bash scripts/hpc/submit_benchmark_pbs.sh; then
  echo "submitted."
else
  echo "COULD NOT SUBMIT from the compute node. Run on the login node:"
  echo "  cd $PROJECT && ARRAY_RANGE=${FIRST_TASK}-${LAST_TASK} FMB_SCENARIOS_PER_TASK=$CHUNK PBS_QUEUE=v1_small72a PBS_WALLTIME=72:00:00 FMB_SCRATCH=$SCRATCH bash scripts/hpc/submit_benchmark_pbs.sh"
  exit 1
fi
