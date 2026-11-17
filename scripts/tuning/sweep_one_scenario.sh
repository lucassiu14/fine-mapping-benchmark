#!/bin/bash
# Sweep ONE parameter on ONE scenario at the production region count.
#
# The point: a canary is 3-8h, but one scenario at 10 regions is ~300s, so a
# six-point sweep is half an hour. Diagnosis belongs here; canaries are for
# CONFIRMING a fix, not finding one.
#
# Default target is the worst cell in canary 4581's paired table: sparse
# scenario 7 (S=1, phi=0.6), where sum(PIP) ~ 11.8 with ONE true causal and AP
# is the lowest of the sparse scenarios (0.395). Under sparse there is no
# infinitesimal background, so the model is correctly specified and that
# overshoot is an inference failure, not misspecification - and it is the most
# direct candidate cause of the FDR floor.
#
#   bash scripts/tuning/sweep_one_scenario.sh                      # sigma_sq
#   SWEEP=sparse_concrete VALS=10,25,51,120,300 bash scripts/tuning/sweep_one_scenario.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

PROJECT_ROOT="$(pwd)"
TUNE_DIR="${TUNE_DIR:-${PROJECT_ROOT}/tuning}"
LOG_DIR="${FMB_LOG_DIR:-${TUNE_DIR}/logs}"
mkdir -p "$LOG_DIR"

STRATUM="${STRATUM:-sparse}"
case "$STRATUM" in
  sparse)     SIM="${SIM:-${TUNE_DIR}/sim_sparse_binary_e5.4_n5000.rds}" ;;
  sparse_inf) SIM="${SIM:-${TUNE_DIR}/sim_sparseinf_rho0.2_binary_e5.4_n5000.rds}" ;;
  *) echo "STRATUM must be sparse or sparse_inf" >&2; exit 1 ;;
esac
[[ -f "$SIM" ]] || { echo "sim not found: $SIM" >&2; exit 1; }

SCENARIO="${SCENARIO:-7}"
REGIONS="${REGIONS:-10}"
SWEEP="${SWEEP:-sigma_sq}"
VALS="${VALS:-0.005,0.02,0.0899,0.3,1.0}"
MAX_ITER="${MAX_ITER:-2000}"
METHOD="${METHOD:-fb_xregion}"
R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
PY_MODULE="${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
FB_PY="${FB_PY:-$HOME/tools/py-venv-runner.sh}"
WALL="${WALL:-04:00:00}"
SELECT="${SELECT:-1:ncpus=1:mem=16gb}"

# Completion OFF: this isolates the overshoot in the model as it has always
# run. Completion raised mass on 26/26 scenarios, so leaving it on would
# confound the thing being measured.
BASE_PARAMS="--sigma_sq 0.0899 --n_caus 10 --sparse_concrete 51 --prior_regularisation 0.1255 --lambda_l1 0.1223 --hierarchy_M 10.15 --gamma 0.1 --gamma_key 0.2 --prior_mass 3 --posterior_completion 0"

J="$(mktemp -t fbswp_XXXXXX.sh)"
cat > "$J" <<PBS_EOF
#!/bin/bash
#PBS -N fbsweep
#PBS -l select=${SELECT}
#PBS -l walltime=${WALL}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -uo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
module load ${PY_MODULE}
echo "[sweep on \$(hostname)] start \$(date)"
echo "commit  : \$(git rev-parse --short HEAD)"
echo "target  : ${STRATUM} scenario ${SCENARIO}, ${REGIONS} regions, completion OFF"
echo "sweeping: ${SWEEP} over ${VALS}"
echo
printf '%-16s %-6s %8s %8s %10s %8s %8s %7s\n' "$SWEEP" "S" "AP" "FDR@.95" "sum/causal" "ece_hi" "n@.95" "secs"
IFS=',' read -ra VV <<< "${VALS}"
for v in "\${VV[@]}"; do
  OUT="\$(Rscript scripts/tuning/fb_objective.R \\
      --sim "${SIM}" --stratum "${STRATUM}" --method "${METHOD}" \\
      --scenario ${SCENARIO} --regions ${REGIONS} \\
      --beatrice_dir "${PROJECT_ROOT}/BEATRICE_annot_sparse" \\
      --python "${FB_PY}" --max_iter ${MAX_ITER} \\
      ${BASE_PARAMS} --${SWEEP} "\$v" 2>&1)"
  g() { printf '%s\\n' "\$OUT" | sed -n "s/^\$1=//p" | tail -1; }
  printf '%-16s %-6s %8s %8s %10s %8s %8s %7s\n' \\
    "\$v" "\$(g scen_S)" "\$(g ap)" "\$(g fdr_at_95)" "\$(g total_mass_ratio)" \\
    "\$(g ece_hi)" "\$(g n_failed)" "\$(g seconds)"
  if [[ "\$(g status)" != "ok" ]]; then
    echo "   !! status=\$(g status) message=\$(g message)"
  fi
done
echo
echo "[sweep] done \$(date)"
PBS_EOF

bash -n "$J" || { echo "generated script has a syntax error: $J" >&2; exit 1; }
echo "Sweep script: $J"
JOBID="$(qsub "$J")" || { echo "qsub failed" >&2; exit 1; }
echo "Sweep job: ${JOBID}"
echo
echo "Reads as one table. ~${#VALS}-point sweep at ~300s each. Read it with:"
echo "  cat ${LOG_DIR}/${JOBID}.OU"
echo
echo "The column that matters is sum/causal: it is sum(PIP) over the number of"
echo "TRUE causals, so 11.8 at S=1 means the posterior expects ~12 causal"
echo "variants where there is one. If a setting brings that to ~1 while AP"
echo "holds, that is the FDR floor's cause."
