#!/bin/bash
#PBS -N fmb4verify
#PBS -l select=1:ncpus=1:mem=16gb
#PBS -l walltime=72:00:00
#PBS -j oe
# =============================================================================
# scripts/analysis/verify_and_analyse.sh
#
# NOTE: the #PBS block sits directly under the shebang deliberately. Directives
# must precede the first executable line, and some PBS builds stop scanning at
# the first blank line - burying them below a long comment header is a silent
# way to have every resource request ignored.
#
# ONE SUBMISSION: verify the finished benchmark, and only if it is sound, launch
# the analysis chain.
#
# WHY THIS EXISTS. Every check in here is one that was skipped or done too late
# in Iteration 004, and each omission cost real time:
#
#   * results.rds completeness  - the run looked done at 10,861 evaluation.rds
#     when the file the analysis actually reads was already complete at 11,250.
#   * preflight                 - a decayed Python environment let 7 of 19
#     methods record 40-50% failed fits for hours; nothing stopped the job.
#   * per-method failure rates  - only discovered AFTER building a full variance
#     decomposition on top of an arm that was never measured.
#
# Running them on the login node is slow (11,250 stat calls over RDS, plus
# reading a sample of results.rds), which is why this is a job.
#
#   qsub -v PROJECT=$HOME/fine-mapping-benchmark scripts/analysis/verify_and_analyse.sh
#
# Set DRY=1 to verify and report WITHOUT submitting the analysis.
# =============================================================================
set -uo pipefail

PROJECT="${PROJECT:-$HOME/fine-mapping-benchmark}"
cd "$PROJECT" || { echo "cannot cd to $PROJECT"; exit 1; }
SCRATCH="${FMB_SCRATCH:-${EPHEMERAL}/fmbench_iter004}"
BENCH="${SCRATCH}/results/benchmark"
LOGS="${SCRATCH}/logs/benchmark"
EXPECT="${EXPECT:-11250}"
MAX_FAIL="${MAX_FAIL:-5}"      # % declared failures tolerated per method/arm
MAX_SILENT="${MAX_SILENT:-10}" # % silent all-NA tolerated (finemap_inf exempt)

module load "${R_MODULE:-R/4.5.2-gfbf-2025b}" 2>/dev/null
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

FAIL=0
note() { printf "\n=== %s ===\n" "$1"; }

note "1. completeness"
N_RES=$(find "$BENCH" -name 'results.rds'    | wc -l | tr -d ' ')
N_EVA=$(find "$BENCH" -name 'evaluation.rds' | wc -l | tr -d ' ')
echo "results.rds    : $N_RES of $EXPECT   <- the file the analysis reads"
echo "evaluation.rds : $N_EVA of $EXPECT   <- informational only; finemap_inf's"
echo "                 phi=0.6 all-NA output crashes the OLD evaluator, so a"
echo "                 shortfall here is expected and does not block anything."
if (( N_RES < EXPECT )); then
  echo "FAIL: results.rds incomplete - $((EXPECT - N_RES)) missing"; FAIL=1
else
  echo "ok"
fi

note "2. python preflight across every task"
PF_OK=$(grep -lh "preflight ok"     "$LOGS"/* 2>/dev/null | wc -l | tr -d ' ')
PF_NO=$(grep -lh "PREFLIGHT FAILED" "$LOGS"/* 2>/dev/null | wc -l | tr -d ' ')
echo "tasks reporting preflight ok     : $PF_OK"
echo "tasks reporting PREFLIGHT FAILED : $PF_NO"
if (( PF_NO > 0 )); then
  echo "FAIL: the toolchain was broken on at least one node"; FAIL=1
else
  echo "ok"
fi

note "3. per-method failure rates"
GATE_OUT="$(Rscript scripts/hpc/qc_gate.R "$BENCH" --sample "${QC_SAMPLE:-300}" 2>&1)"
echo "$GATE_OUT"
MAXF=$(printf '%s\n' "$GATE_OUT" | awk '/^GATE_MAX_FAIL_PCT/   {print int($2+0.5)}')
MAXS=$(printf '%s\n' "$GATE_OUT" | awk '/^GATE_MAX_SILENT_PCT/ {print int($2+0.5)}')
MAXF=${MAXF:-999}; MAXS=${MAXS:-999}
if (( MAXF > MAX_FAIL )); then
  echo "FAIL: a method declares ${MAXF}% failed fits (limit ${MAX_FAIL}%)"; FAIL=1
elif (( MAXS > MAX_SILENT )); then
  # Checked SEPARATELY from declared failures on purpose: a method returning an
  # all-NA posterior without setting an error scores 0% failed, which is how
  # CARMA went unnoticed for three iterations.
  echo "FAIL: a method silently returns ${MAXS}% all-NA fits (limit ${MAX_SILENT}%)"; FAIL=1
else
  echo "ok - worst declared ${MAXF}%, worst silent ${MAXS}%"
fi

note "verdict"
if (( FAIL )); then
  echo "NOT SUBMITTING THE ANALYSIS. Fix the above and resubmit this job."
  exit 1
fi
echo "All checks passed."
if [[ "${DRY:-0}" == "1" ]]; then
  echo "DRY=1 set - stopping before submission."
  exit 0
fi

note "4. submitting the analysis chain"
# qsub from a compute node is permitted on this cluster, but if that ever
# changes the failure must be actionable rather than silent - print the exact
# command instead of exiting 0 as though the analysis were running.
if QUEUE_A=v1_small72a QUEUE_B=v1_small72a WALL_A=72:00:00 WALL_B=72:00:00 \
   SELECT_A=1:ncpus=1:mem=12gb SELECT_B=1:ncpus=1:mem=16gb \
   FMB_SCRATCH="$SCRATCH" bash scripts/analysis/submit_iter004_analysis.sh; then
  echo "analysis submitted."
else
  echo "COULD NOT SUBMIT from the compute node. Run this on the login node:"
  echo "  cd $PROJECT && QUEUE_A=v1_small72a QUEUE_B=v1_small72a WALL_A=72:00:00 WALL_B=72:00:00 SELECT_A=1:ncpus=1:mem=12gb SELECT_B=1:ncpus=1:mem=16gb FMB_SCRATCH=$SCRATCH bash scripts/analysis/submit_iter004_analysis.sh"
  exit 1
fi
