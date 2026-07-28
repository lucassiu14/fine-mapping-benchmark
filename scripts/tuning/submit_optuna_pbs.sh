#!/bin/bash
# =============================================================================
# scripts/tuning/submit_optuna_pbs.sh
#
# Submit an array of Optuna workers that all tune Functional BEATRICE against
# ONE shared study. Re-running this script CONTINUES the same study: every
# worker opens the storage with load_if_exists=True, so completed trials are
# remembered and the sampler keeps its history. Nothing is recomputed.
#
# Each worker is given a --timeout of (walltime - TRIAL_MARGIN) so it stops
# starting new trials and exits cleanly before PBS would kill it mid-write.
#
# Usage (from project root):
#   bash scripts/tuning/submit_optuna_pbs.sh                     # defaults
#   WORKERS=20 PBS_WALLTIME=24:00:00 bash scripts/tuning/submit_optuna_pbs.sh
#   OBJECTIVE=scalar bash scripts/tuning/submit_optuna_pbs.sh    # enables pruning
#
# Resume after walltime: just run the same command again.
# =============================================================================
set -euo pipefail

PROJECT_ROOT="$(pwd)"

# --- study identity ----------------------------------------------------------
# One study per stratum (never-pool rule): a prior tuned on binary/ref500 has no
# reason to be optimal on continuous/in-sample. STUDY names the stratum.
STUDY="${STUDY:-fb_binary_ref500}"
SIM="${SIM:-tuning/sim_binary_ref500.rds}"

# Storage must OUTLIVE ephemeral (which auto-purges ~30 days) - the study is the
# accumulated knowledge of the search and is only a few MB. Default: project dir.
STORAGE="${STORAGE:-${PROJECT_ROOT}/tuning/${STUDY}.journal}"

# --- search configuration ----------------------------------------------------
OBJECTIVE="${OBJECTIVE:-multi}"        # multi | scalar
PENALTY="${PENALTY:-0.5}"              # scalar mode only
SCENARIOS="${SCENARIOS:-4}"            # scenarios per trial
REGIONS="${REGIONS:-4}"                # regions per scenario
MAX_ITER="${MAX_ITER:-2000}"           # FB max_iter (fixed, never tuned)

# --- stratum for the tuning sim (used only when SIM has to be built) ---------
# These must describe the SAME stratum that STUDY names.
ANNOTATIONS="${ANNOTATIONS:-binary}"   # none | binary | continuous
ENRICHMENT="${ENRICHMENT:-10.8}"
N_REF="${N_REF:-500}"                  # "NA" for in-sample LD
MODEL="${MODEL:-sparse}"               # sparse | sparse_inf
P_CAUSAL="${P_CAUSAL:-}"               # sparse_inf only
SIM_WALLTIME="${SIM_WALLTIME:-04:00:00}"
SIM_SELECT="${SIM_SELECT:-1:ncpus=1:mem=32gb}"
SIM_QUEUE="${SIM_QUEUE:-v1_small24}"

# --- cluster resources -------------------------------------------------------
WORKERS="${WORKERS:-10}"               # parallel array elements sharing the study
PBS_QUEUE="${PBS_QUEUE:-v1_small72a}"
PBS_WALLTIME="${PBS_WALLTIME:-24:00:00}"
PBS_SELECT="${PBS_SELECT:-1:ncpus=1:mem=16gb}"
R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
# The venv python is dynamically linked against the module's libpython, so the
# Python module MUST be loaded before calling it - otherwise it dies with
# "error while loading shared libraries: libpython3.12.so.1.0".
PY_MODULE="${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
# Python that has optuna installed (see the doc for the one-off pip install).
PY="${PY:-$HOME/tools/fmpy-venv/bin/python}"
# Python that runs BEATRICE itself (torch); the R wrapper shells out to this.
FB_PY="${FB_PY:-$HOME/tools/py-venv-runner.sh}"

LOG_DIR="${FMB_LOG_DIR:-${PROJECT_ROOT}/tuning/logs}"
mkdir -p "$LOG_DIR" "$(dirname "$STORAGE")"

# --- walltime -> worker timeout ---------------------------------------------
# Leave one typical trial's worth of margin so a worker never gets killed
# mid-trial: it finishes what it has, exits, and the next submit resumes.
hh=${PBS_WALLTIME%%:*}; rest=${PBS_WALLTIME#*:}; mm=${rest%%:*}; ss=${rest##*:}
WALL_SEC=$((10#$hh*3600 + 10#$mm*60 + 10#$ss))
TRIAL_MARGIN="${TRIAL_MARGIN:-1800}"                 # 30 min default
TIMEOUT=$((WALL_SEC - TRIAL_MARGIN))
if (( TIMEOUT < 300 )); then
  echo "ERROR: walltime ${PBS_WALLTIME} is too short for TRIAL_MARGIN=${TRIAL_MARGIN}s" >&2
  exit 1
fi

# --- bootstrap: build the tuning sim if it is missing ------------------------
# Submitted as its own job, with the worker array made dependent on it
# (-W depend=afterok), so ONE invocation of this script does the whole thing:
# build-if-needed, then tune. If the sim already exists the build is skipped and
# the array starts immediately - which is also what happens on every resume.
DEPEND=""
if [[ -f "$SIM" ]]; then
  echo "Tuning sim: $SIM (exists - reusing)"
else
  echo "Tuning sim: $SIM (missing - submitting a build job first)"
  SIM_SCRIPT="$(mktemp -t fmbtunesim_XXXXXX.sh)"
  PC_ARG=""
  [[ -n "$P_CAUSAL" ]] && PC_ARG="--p_causal ${P_CAUSAL}"
  cat > "$SIM_SCRIPT" <<SIM_EOF
#!/bin/bash
#PBS -N fbtunesim
#PBS -q ${SIM_QUEUE}
#PBS -l select=${SIM_SELECT}
#PBS -l walltime=${SIM_WALLTIME}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
echo "[sim build on \$(hostname)] start \$(date)"
Rscript scripts/tuning/make_tuning_sim.R \
    --out "${SIM}" \
    --annotations "${ANNOTATIONS}" \
    --enrichment ${ENRICHMENT} \
    --n_ref "${N_REF}" \
    --model "${MODEL}" \
    --regions ${REGIONS} ${PC_ARG}
echo "[sim build] done \$(date)"
SIM_EOF
  SIM_JOBID="$(qsub "$SIM_SCRIPT")"
  echo "  sim build job: ${SIM_JOBID}"
  # afterok: the workers only start if the build SUCCEEDS. If it fails they stay
  # held and can be deleted - they never run against a missing/partial sim.
  DEPEND="-W depend=afterok:${SIM_JOBID}"
fi

echo "Study:     $STUDY   (objective=$OBJECTIVE)"
echo "Sim:       $SIM"
echo "Storage:   $STORAGE"
echo "Workers:   $WORKERS x ${PBS_WALLTIME}  (worker timeout ${TIMEOUT}s)"
echo "Log dir:   $LOG_DIR"

JOB_SCRIPT="$(mktemp -t fmbtune_XXXXXX.sh)"
cat > "$JOB_SCRIPT" <<PBS_EOF
#!/bin/bash
#PBS -N fbtune
#PBS -q ${PBS_QUEUE}
#PBS -l select=${PBS_SELECT}
#PBS -l walltime=${PBS_WALLTIME}
#PBS -J 1-${WORKERS}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/

set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
module load ${PY_MODULE}

echo "[worker \${PBS_ARRAY_INDEX} on \$(hostname)] start \$(date)"
${PY} scripts/tuning/optuna_fb.py \
    --sim "${SIM}" \
    --study "${STUDY}" \
    --storage "${STORAGE}" \
    --objective "${OBJECTIVE}" \
    --penalty ${PENALTY} \
    --scenarios ${SCENARIOS} \
    --regions ${REGIONS} \
    --max_iter ${MAX_ITER} \
    --timeout ${TIMEOUT} \
    --python "${FB_PY}"
echo "[worker \${PBS_ARRAY_INDEX}] done \$(date)"
PBS_EOF

echo "Job script: $JOB_SCRIPT"
# ${DEPEND} is intentionally unquoted: it is either empty or "-W depend=afterok:ID".
# shellcheck disable=SC2086
qsub ${DEPEND} "$JOB_SCRIPT"
echo
if [[ -n "$DEPEND" ]]; then
  echo "NOTE: the worker array is HELD until the sim build finishes successfully"
  echo "      (shows as state H in qstat, then switches to Q/R by itself)."
fi
echo "Track with:   qstat -tan -u \$USER | grep -E 'fbtune|fbtunesim'"
echo "Progress:     ${PY} scripts/tuning/optuna_fb.py --sim ${SIM} --study ${STUDY} \\"
echo "                  --storage ${STORAGE} --objective ${OBJECTIVE} --report"
echo "Resume later: re-run this exact script (continues the same study)"
