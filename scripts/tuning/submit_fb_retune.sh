#!/bin/bash
# =============================================================================
# scripts/tuning/submit_fb_retune.sh
#
# Re-tune fb_xregion after the reformat_memo fix, over BOTH generative models.
#
# WHY THIS EXISTS, rather than reusing submit_optuna_pbs.sh:
#
#   1. The shipped configuration was tuned on 2026-08-06 (31b4794) against a
#      saturated memo - reformat_memo subtracted the MEAN log-BF, so the top
#      half of every memo clipped to one value and the posterior was flattened.
#      The fix (max offset) reached the benchmark on 2026-11-04 (119f6c6), so
#      every hyperparameter was selected against an objective computed on a
#      posterior that no longer exists.
#
#   2. fb_objective.R scored "functional_beatrice" - the PER-REGION method.
#      run_benchmark_job.R then assigns that configuration to fb_xregion
#      (METHOD_ARGS$fb_xregion <- .FBv). lambda_l1, prior_regularisation and
#      hierarchy_M are the SHARED prior head's knobs, so the cross-region
#      method has never been tuned at all - it inherited another method's
#      settings. This study tunes fb_xregion itself (--method default).
#
#   3. The old grid was sparse-only. Measured on Iteration 009, 95.6% of
#      fb_xregion's false calls at PIP >= 0.95 are in sparse_inf, so a
#      sparse-only grid contains 4.4% of the failure being tuned against -
#      the v1 mistake docs/autoresearch/hyperparameter-tuning.md records.
#
# GRID (measured placement of the failure, Iteration 009):
#   model      sparse AND sparse_inf, binary annotations
#   S          1, 3, 10         FDR 0.866 at S=1 -> 0.572 at S=10
#   phi        0.05, 0.2, 0.6   FDR 0.318 -> 0.816; phi=0.6 holds 43% of calls
#   rho        0.2 and 0.6      FDR 0.879 -> 0.686 (sparse_inf only)
#   regions    10 per scenario, 2 each of 500/750/1000/1500/2000
#   enrichment 5.4 (flat axis: FDR 0.319-0.354 across all four folds)
#   n = 5000, in-sample LD
# => 9 + 9 + 9 = 27 scenarios x 10 regions = 270 region-fits per trial.
#
# REGIONS IS NOT THE AXIS TO ECONOMISE ON. fb_xregion trains ONE prior head
# over all of a scenario's regions, so regions-per-scenario IS the amount of
# cross-region pooling, and an L1 penalty tuned on 4 regions is not the same
# quantity deployed on 10. generate_params_grid.R uses 10; so does this. If a
# trial is too slow, cut S or phi levels, never regions.
#
# NEVER POOLED. sparse and sparse_inf are separate strata: each yields its own
# AP and its own FDR (pooled from its own counts), and the objective takes the
# WORST stratum's FDR. The two rho sims share the sparse_inf label because the
# benchmark pools over rho within a stratum; give them distinct labels to
# target the worst rho instead.
#
#   CANARY=1 bash scripts/tuning/submit_fb_retune.sh   # 1 trial, for timing
#   WORKERS=20 bash scripts/tuning/submit_fb_retune.sh # the study
# =============================================================================
set -euo pipefail

PROJECT_ROOT="$(pwd)"
[[ -f scripts/tuning/optuna_fb.py ]] || { echo "ERROR: run from the project root" >&2; exit 1; }

STUDY="${STUDY:-fb_retune_cells_v4}"
TUNE_DIR="${TUNE_DIR:-${PROJECT_ROOT}/tuning}"
STORAGE="${STORAGE:-${TUNE_DIR}/${STUDY}.journal}"
LOG_DIR="${FMB_LOG_DIR:-${TUNE_DIR}/logs}"
METHOD="${METHOD:-fb_xregion}"
# v3 used OBJECTIVE=triple, which worst-cased FDR over strata but pooled over
# the true causal count S. Its winning point was a two-variant detector:
# sparse_concrete=2, five times below the published floor of 10, +26% AP at
# S=1 and -48% at S=10. v4 worst-cases over (stratum x S) cells and makes
# posterior mass a hard constraint instead of a printed warning.
OBJECTIVE="${OBJECTIVE:-cells}"
# The published BEATRICE declares sparse_concrete with lower_bound=10 and
# refuses to start below it. Do not lower SC_MIN.
SC_MIN="${SC_MIN:-10}"
SC_MAX="${SC_MAX:-300}"
MASS_LO="${MASS_LO:-0.5}"
MASS_HI="${MASS_HI:-2.0}"
# prior_mass: the per-region anchor on sum(p_0), added 2026-11-14. Until it
# existed, mass was the one constrained quantity that no search parameter
# controlled directly, which is why v4 found 0/76 feasible trials. PM_MAX=0
# keeps the anchor off and out of the search space (the v1-v4 behaviour).
PM_MIN="${PM_MIN:-0.5}"
PM_MAX="${PM_MAX:-0}"
CANARY="${CANARY:-0}"

# --- grid -------------------------------------------------------------------
SIM_S="${SIM_S:-1,3,10}"
SIM_PHI="${SIM_PHI:-0.05,0.2,0.6}"
SIM_N="${SIM_N:-5000}"
SIM_REGION_SIZES="${SIM_REGION_SIZES:-500,750,1000,1500,2000}"
REGIONS="${REGIONS:-10}"
ENRICHMENT="${ENRICHMENT:-5.4}"
RHOS="${RHOS:-0.2 0.6}"
SEED="${SEED:-20261110}"           # one seed for every sim: the strata then
                                   # differ only in the generative model
# Scenario count per sim MUST equal |S| x |phi|, or trials silently score only
# a prefix of the grid - the v1 blindness.
N_S=$(( $(awk -F, '{print NF}' <<< "$SIM_S") ))
N_PHI=$(( $(awk -F, '{print NF}' <<< "$SIM_PHI") ))
PER_SIM=$(( N_S * N_PHI ))

# --- queues: EMPTY by default. v1_small72a/v1_small24 deny access, and naming
# one kills the submission. Let the cx routing queue place the job. -----------
PBS_QUEUE="${PBS_QUEUE:-}"
SIM_QUEUE="${SIM_QUEUE:-}"
qline() { if [[ -n "${1:-}" ]]; then echo "#PBS -q $1"; else
  echo "# no queue named: the cx routing queue places this job"; fi; }

PBS_SELECT="${PBS_SELECT:-1:ncpus=1:mem=16gb}"
PBS_WALLTIME="${PBS_WALLTIME:-24:00:00}"
SIM_SELECT="${SIM_SELECT:-1:ncpus=1:mem=48gb}"
SIM_WALLTIME="${SIM_WALLTIME:-06:00:00}"
R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
PY_MODULE="${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
PY="${PY:-$HOME/tools/fmpy-venv/bin/python}"
FB_PY="${FB_PY:-$HOME/tools/py-venv-runner.sh}"
MAX_ITER="${MAX_ITER:-2000}"
WORKERS="${WORKERS:-10}"

# Trial cost scales with the number of REGION FITS, which is
# scenarios-per-sim x sims x regions, so the measurement is kept per fit and
# multiplied up below once the strata are known. The 2026-11-10 canary (job
# 4596630) ran 27 scenarios x 10 regions = 270 fits in 175.8 minutes: 10548s a
# trial, ~390s a scenario, ~39s a region fit.
#
# Keeping it per fit matters now that STRATA can select one stratum: a
# sparse-only study is 9 scenarios x 10 regions = 90 fits, a third of the work,
# and a hardcoded 10548 would reserve three and a half hours of margin for a
# one-hour trial and launch a third as many trials as the walltime allows.
MEASURED_FIT_SEC="${MEASURED_FIT_SEC:-39}"
mkdir -p "$TUNE_DIR" "$LOG_DIR"

# The venv python is dynamically linked against the module's libpython, so the
# module has to be loaded before it will even start (submit_optuna_pbs.sh:72).
# Checked here because a missing optuna is only discovered inside the job, which
# on 2026-11-10 cost a 90-minute sim build and a queue round trip to learn.
if ! ( module load "${PY_MODULE}" 2>/dev/null; "$PY" -c 'import optuna' ) 2>/dev/null; then
  cat >&2 <<MSG
ERROR: $PY cannot import optuna, so no trial can run. Install it:

  module load ${PY_MODULE} && ${PY%/python}/pip install optuna scikit-learn

scikit-learn is optional - without it --report just omits the fANOVA parameter
importance table (optuna_fb.py guards that call).
MSG
  exit 1
fi

# --- which strata go into THIS study ----------------------------------------
# STRATA selects them. make_tuning_sim.R's own header says to run ONE study per
# stratum, because a configuration optimal under `sparse` has no reason to be
# optimal under `sparse_inf` and the two are never pooled. v4 put all three sims
# in one study; its binding constraint (posterior mass) then pooled them, and no
# single point could satisfy a window on the average - 0/76 feasible trials.
#   STRATA="sparse"              one study on the sparse stratum
#   STRATA="sparse_inf"          one study on sparse_inf (all RHOS)
#   STRATA="sparse sparse_inf"   both in one study (the v4 layout); metrics stay
#                                per stratum and mass is now constrained on the
#                                WORST stratum, but the search must then serve
#                                two generative models with one configuration.
STRATA="${STRATA:-sparse}"
SIM_SPARSE="${TUNE_DIR}/sim_sparse_binary_e${ENRICHMENT}_n${SIM_N}.rds"
SIM_ARGS=()
BUILD=()
for st in $STRATA; do
  case "$st" in
    sparse)
      SIM_ARGS+=( "sparse=${SIM_SPARSE}" )
      [[ -f "$SIM_SPARSE" ]] || BUILD+=( "sparse::${SIM_SPARSE}" )
      ;;
    sparse_inf)
      for rho in $RHOS; do
        f="${TUNE_DIR}/sim_sparseinf_rho${rho}_binary_e${ENRICHMENT}_n${SIM_N}.rds"
        SIM_ARGS+=( "sparse_inf=${f}" )
        [[ -f "$f" ]] || BUILD+=( "sparse_inf:${rho}:${f}" )
      done
      ;;
    *)
      echo "ERROR: unknown stratum '$st' in STRATA (want sparse and/or sparse_inf)" >&2
      exit 1
      ;;
  esac
done
if [[ ${#SIM_ARGS[@]} -eq 0 ]]; then
  echo "ERROR: STRATA selected no sims" >&2; exit 1
fi

N_SIMS=${#SIM_ARGS[@]}
echo "Study      : $STUDY   (method=$METHOD, objective=$OBJECTIVE)"
echo "Grid       : S={$SIM_S} x phi={$SIM_PHI} = ${PER_SIM} scenarios per sim"
echo "Strata     : ${N_SIMS} sims -> $(( PER_SIM * N_SIMS )) scenarios x ${REGIONS} regions"
echo "             = $(( PER_SIM * N_SIMS * REGIONS )) region-fits per trial"
for s in "${SIM_ARGS[@]}"; do echo "             $s"; done
echo "Storage    : $STORAGE"
echo "Objective  : ${OBJECTIVE}  sparse_concrete in [${SC_MIN}, ${SC_MAX}]"
echo "             mass constraint: sum(PIP)/#causal in [${MASS_LO}, ${MASS_HI}]"
if [[ "${PM_MAX}" != "0" ]]; then
  echo "             prior_mass in [${PM_MIN}, ${PM_MAX}] - anchors sum(p_0) per region"
else
  echo "             prior_mass: OFF (learned prior scale left unanchored)"
fi
FITS_PER_TRIAL=$(( PER_SIM * N_SIMS * REGIONS ))
MEASURED_TRIAL_SEC="${MEASURED_TRIAL_SEC:-$(( MEASURED_FIT_SEC * FITS_PER_TRIAL ))}"
# --timeout goes to study.optimize(), which stops LAUNCHING trials once it
# elapses - a trial already running is not interrupted. So the margin has to
# exceed one whole trial, or every worker's last trial is launched near the
# deadline, killed mid-flight by PBS, and left RUNNING in the journal forever.
# 1.2x gives a trial that runs slower than measured somewhere to land.
TRIAL_MARGIN="${TRIAL_MARGIN:-$(( MEASURED_TRIAL_SEC * 12 / 10 ))}"
WALL_S=$(awk -F: '{print $1*3600+$2*60+$3}' <<< "$PBS_WALLTIME")
if (( TRIAL_MARGIN < MEASURED_TRIAL_SEC )); then
  echo "ERROR: TRIAL_MARGIN=${TRIAL_MARGIN}s is less than one trial" >&2
  echo "       (${MEASURED_TRIAL_SEC}s measured). Every worker would lose its" >&2
  echo "       last trial to the walltime kill. Raise it." >&2
  exit 1
fi
PER_WORKER=$(( (WALL_S - TRIAL_MARGIN) / MEASURED_TRIAL_SEC ))
echo "Budget     : trial ~$(( MEASURED_TRIAL_SEC / 60 )) min measured; walltime ${PBS_WALLTIME}"
echo "             ${PER_WORKER} trials/worker x ${WORKERS} workers = $(( PER_WORKER * WORKERS )) trials"
echo "             NSGA-II population is 50, so that is about \
$(( PER_WORKER * WORKERS / 50 )) generations of selection"

DEPEND=""
if (( ${#BUILD[@]} > 0 )); then
  echo "Sims to build: ${#BUILD[@]}"
  S_SCRIPT="$(mktemp -t fbretsim_XXXXXX.sh)"
  { cat <<SIM_EOF
#!/bin/bash
#PBS -N fbretsim
$(qline "$SIM_QUEUE")
#PBS -l select=${SIM_SELECT}
#PBS -l walltime=${SIM_WALLTIME}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
echo "[sim build on \$(hostname)] start \$(date)"
SIM_EOF
    for spec in "${BUILD[@]}"; do
      model="${spec%%:*}"; rest="${spec#*:}"; rho="${rest%%:*}"; out="${rest#*:}"
      pc=""; [[ -n "$rho" ]] && pc="--p_causal ${rho}"
      cat <<SIM_EOF
Rscript scripts/tuning/make_tuning_sim.R \\
    --out "${out}" --model "${model}" ${pc} \\
    --annotations binary --enrichment ${ENRICHMENT} --n_annotations 10 \\
    --n_ref NA --n ${SIM_N} --regions ${REGIONS} \\
    --region_sizes "${SIM_REGION_SIZES}" \\
    --S "${SIM_S}" --phi "${SIM_PHI}" --seed ${SEED}
SIM_EOF
    done
    echo 'echo "[sim build] done $(date)"'
  } > "$S_SCRIPT"
  SIM_JOBID="$(qsub "$S_SCRIPT")"
  echo "  sim build job: ${SIM_JOBID}"
  DEPEND="-W depend=afterok:${SIM_JOBID}"
else
  echo "All sims present - reusing."
fi

SIM_FLAGS=""; for s in "${SIM_ARGS[@]}"; do SIM_FLAGS+=" --sim \"$s\""; done

if [[ "$CANARY" == "1" ]]; then
  # One trial at the shipped configuration, to measure a trial's wall-clock
  # before any budget is committed. --fix pins every parameter, so the trial
  # measures cost, not search.
  WALL="${CANARY_WALLTIME:-08:00:00}"
  J="$(mktemp -t fbretcan_XXXXXX.sh)"
  cat > "$J" <<PBS_EOF
#!/bin/bash
#PBS -N fbretcan
$(qline "$PBS_QUEUE")
#PBS -l select=${PBS_SELECT}
#PBS -l walltime=${WALL}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
module load ${PY_MODULE}
echo "[canary on \$(hostname)] start \$(date)"
${PY} scripts/tuning/optuna_fb.py ${SIM_FLAGS} \\
    --study "${STUDY}_canary" --storage "${TUNE_DIR}/${STUDY}_canary.journal" \\
    --method "${METHOD}" --objective "${OBJECTIVE}" \\
    --scenarios ${PER_SIM} --regions ${REGIONS} --max_iter ${MAX_ITER} \\
    --sc-min ${SC_MIN} --sc-max ${SC_MAX} \\
    --mass-lo ${MASS_LO} --mass-hi ${MASS_HI} \\
    --pm-min ${PM_MIN} --pm-max ${PM_MAX} \\
    --n-trials 1 --python "${FB_PY}" \\
    --fix sigma_sq=0.0899 n_caus=10 sparse_concrete=51 \\
          prior_regularisation=0.1255 lambda_l1=0.1223 hierarchy_M=10.15 \\
          gamma=0.1 gamma_key=0.2
echo "[canary] done \$(date)"
PBS_EOF
  echo "Canary job script: $J"
  # shellcheck disable=SC2086
  JID="$(qsub ${DEPEND} "$J")"
  echo "Canary job: $JID"
  echo
  echo "The canary pins every parameter at the SHIPPED values, so its per-scenario"
  echo "seconds= lines are the cost of one trial, not of the search. This cluster"
  echo "writes <jobid>.OU / .ER into the log dir, so read them with:"
  echo "  grep 'scenario' ${LOG_DIR}/${JID}.OU | tail -30"
  echo "  tail -20 ${LOG_DIR}/${JID}.ER      # on failure"
  exit 0
fi

J="$(mktemp -t fbretune_XXXXXX.sh)"
cat > "$J" <<PBS_EOF
#!/bin/bash
#PBS -N fbretune
$(qline "$PBS_QUEUE")
#PBS -l select=${PBS_SELECT}
#PBS -l walltime=${PBS_WALLTIME}
#PBS -J 1-${WORKERS}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
module load ${PY_MODULE}
WALL_SEC=\$(( $(awk -F: '{print $1*3600+$2*60+$3}' <<< "$PBS_WALLTIME") ))
echo "[worker \${PBS_ARRAY_INDEX} on \$(hostname)] start \$(date)"
${PY} scripts/tuning/optuna_fb.py ${SIM_FLAGS} \\
    --study "${STUDY}" --storage "${STORAGE}" \\
    --method "${METHOD}" --objective "${OBJECTIVE}" \\
    --scenarios ${PER_SIM} --regions ${REGIONS} --max_iter ${MAX_ITER} \\
    --sc-min ${SC_MIN} --sc-max ${SC_MAX} \\
    --mass-lo ${MASS_LO} --mass-hi ${MASS_HI} \\
    --pm-min ${PM_MIN} --pm-max ${PM_MAX} \\
    --timeout \$(( WALL_SEC - ${TRIAL_MARGIN} )) --python "${FB_PY}"
echo "[worker \${PBS_ARRAY_INDEX}] done \$(date)"
PBS_EOF
echo "Job script: $J"
# shellcheck disable=SC2086
JID="$(qsub ${DEPEND} "$J")"
echo "Array:    $JID"
echo "Logs:     ${LOG_DIR}/<jobid>.OU and .ER  (this cluster's naming)"
echo
echo "Track:    qstat -tan -u \$USER | grep -E 'fbretune|fbretsim'"
echo "Report:   ${PY} scripts/tuning/optuna_fb.py ${SIM_FLAGS} --study ${STUDY} \\"
echo "              --storage ${STORAGE} --objective ${OBJECTIVE} --report"
