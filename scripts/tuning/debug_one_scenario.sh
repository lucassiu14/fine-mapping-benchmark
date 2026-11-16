#!/bin/bash
# Run ONE scenario through fb_objective.R on a compute node and show EVERYTHING.
#
# Why this exists. optuna_fb.run_one_scenario invokes fb_objective.R with
# subprocess.run(capture_output=True), so the R script's stdout - which carries
# status=, message=, n_failed= - and the Python trainer's stderr are captured and
# then DISCARDED whenever status is 'ok'. Canary 956 completed 27 scenarios that
# way, returned the all-NaN sentinel ([0.0, 1.0, 1.0] for objective 'cells') and
# left no record anywhere of why. This calls the same script directly so the same
# bytes land in the log.
#
# It runs the scenario TWICE at an identical configuration, with posterior
# completion ON and then OFF, so the two outputs isolate whether completion is
# responsible. REGIONS defaults to 2 to keep it short: the metric VALUES are
# therefore not comparable with a 10-region canary, which is fine - what is being
# read here is status/message/n_failed and the presence of [completion] lines.
#
#   bash scripts/tuning/debug_one_scenario.sh
#   REGIONS=4 SCENARIO=5 STRATUM=sparse_inf bash scripts/tuning/debug_one_scenario.sh
set -euo pipefail
cd "$(dirname "$0")/../.." || exit 1

PROJECT_ROOT="$(pwd)"
TUNE_DIR="${TUNE_DIR:-${PROJECT_ROOT}/tuning}"
LOG_DIR="${FMB_LOG_DIR:-${TUNE_DIR}/logs}"
mkdir -p "$LOG_DIR"

STRATUM="${STRATUM:-sparse}"
case "$STRATUM" in
  sparse)     SIM="${SIM:-${TUNE_DIR}/sim_sparse_binary_e5.4_n5000.rds}" ;;
  sparse_inf) SIM="${SIM:-${TUNE_DIR}/sim_sparseinf_rho0.2_binary_e5.4_n5000.rds}" ;;
  *) echo "STRATUM must be sparse or sparse_inf (got '$STRATUM')" >&2; exit 1 ;;
esac
[[ -f "$SIM" ]] || { echo "sim not found: $SIM" >&2; ls -la "$TUNE_DIR"/*.rds 2>/dev/null; exit 1; }

SCENARIO="${SCENARIO:-1}"
REGIONS="${REGIONS:-2}"
MAX_ITER="${MAX_ITER:-2000}"
METHOD="${METHOD:-fb_xregion}"
R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
PY_MODULE="${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
FB_PY="${FB_PY:-$HOME/tools/py-venv-runner.sh}"
WALL="${WALL:-01:30:00}"
SELECT="${SELECT:-1:ncpus=1:mem=16gb}"
QUEUE="${QUEUE:-}"

# The canary's pinned configuration, prior_mass included (submit_fb_retune.sh:98
# plus the prior_mass=3 the canary banner reported).
PARAMS="${PARAMS:---sigma_sq 0.0899 --n_caus 10 --sparse_concrete 51 --prior_regularisation 0.1255 --lambda_l1 0.1223 --hierarchy_M 10.15 --gamma 0.1 --gamma_key 0.2 --prior_mass 3}"
TOL="${TOL:-0.001}"

J="$(mktemp -t fbdbg1_XXXXXX.sh)"
cat > "$J" <<PBS_EOF
#!/bin/bash
#PBS -N fbdbg1
$(if [[ -n "$QUEUE" ]]; then echo "#PBS -q $QUEUE"; else echo "# no queue named: the cx routing queue places this job"; fi)
#PBS -l select=${SELECT}
#PBS -l walltime=${WALL}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -uo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
module load ${PY_MODULE}
echo "[debug1 on \$(hostname)] start \$(date)"
echo "sim      : ${SIM}"
echo "stratum  : ${STRATUM}   scenario ${SCENARIO}   regions ${REGIONS}"
echo "params   : ${PARAMS}"
echo "git HEAD : \$(git rev-parse --short HEAD)"
echo
echo "checking the installed wrapper actually carries the flag:"
Rscript -e 'cat(any(grepl("posterior_completion", deparse(fmbenchmark:::.fb_joint_scenario_setup))), "\n")'
echo

for TOLV in ${TOL} 0; do
  echo "================================================================"
  echo "=== posterior_completion = \${TOLV} ==="
  echo "================================================================"
  Rscript scripts/tuning/fb_objective.R \\
      --sim "${SIM}" --stratum "${STRATUM}" --method "${METHOD}" \\
      --scenario ${SCENARIO} --regions ${REGIONS} \\
      --beatrice_dir "${PROJECT_ROOT}/BEATRICE_annot_sparse" \\
      --python "${FB_PY}" --max_iter ${MAX_ITER} \\
      ${PARAMS} \\
      --posterior_completion \${TOLV} \\
      --completion_seeds 8 --completion_width 32 --completion_rounds 4
  echo "--- Rscript exit status: \$? ---"
  echo
done
echo "[debug1] done \$(date)"
PBS_EOF

bash -n "$J" || { echo "generated job script has a syntax error: $J" >&2; exit 1; }
echo "Debug job script: $J"
JOBID="$(qsub "$J")"
echo "Debug job: ${JOBID}"
echo
# Print the REAL log path, fully substituted. This used to echo a literal
# ${JOBID%%.*}, which does not expand in the reader's shell - so the natural
# fallback is `ls -t *.OU | head -1`, and that silently returns the PREVIOUS
# job's log while this one is still queued. Two runs then look identical and a
# gate appears to pass when it was never checked.
LOGF="${LOG_DIR}/${JOBID}.OU"
echo "Reads EVERYTHING fb_objective.R prints, which optuna_fb normally captures"
echo "and throws away. This job's log, and no other (~4-60 min):"
echo "  grep -nE 'status=|n_failed=|^ap=|^fdr_at_95=|^total_mass_ratio=|\\[completion\\]' ${LOGF}"
echo "  tail -60 ${LOGF}"
echo
echo "If that path does not exist yet the job has not started - wait, do NOT"
echo "fall back to the newest .OU, which is a different job."
