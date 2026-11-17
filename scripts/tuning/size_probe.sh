#!/bin/bash
# Measure the WEIGHTED CONFIGURATION-SIZE DISTRIBUTION the reported PIPs are
# formed from, for ONE scenario, with posterior completion OFF.
#
# Why this exists. sum_j PIP_j = E[|gamma| | gamma in M] exactly (calculate_pip
# normalises over the visited set M), so the number of variants that CAN hold a
# high PIP is set by how large the configurations in M are - and therefore the
# FDR@0.95 floor is too. Three candidates were tested against that:
#
#   sigma_sq   REFUTED locally. On a real 1000G locus (600 variants, 170 kb,
#              463 pairs at |r|>0.99) the greedy argmax of the real log-BF is a
#              SINGLETON at every kappa = n*sigma_sq from 10 to 450, when z is
#              drawn consistently with the LD. The near-null eigendirection leak
#              (0.5*kappa*c^2/(1+kappa*lambda), linear in kappa at lambda~0) is
#              real but needs z INCONSISTENT with R - i.e. a misspecified panel,
#              which is iteration 007's regime, not the in-sample sparse cells.
#   prior_mass REFUTED locally, and this reproduces the 411-trial observation:
#              sweeping it 0.5 -> 5.0 moved sum(PIP) from 1.01 to 1.08 at |z|=3
#              and 1.05 to 1.38 at |z|=8, with ZERO false positives at 0.95.
#              anchor_prior_mass fixes sum(p0), which is the PRIOR's E|gamma|;
#              it does not control E[|gamma| | M] because M is not the prior's
#              support.
#   gamma      NOT yet measured, and the only remaining candidate that acts on
#              configuration SIZE directly: joint_trainer.py:141-144 and :395
#              build each sampled configuration as
#                  |gamma| = #{j in top-`sparse_concrete` : c_j > eps}, eps=gamma
#              so gamma is the inclusion threshold on the Concrete draw. It has
#              been pinned at the Python default 0.1 since it was first exposed
#              (2026-11-10) and NEVER searched. Two independent records point the
#              same way: study v3 walked `sparse_concrete` to 2 and scored "5x
#              better FDR" (optuna_fb.py:205-210), and v4's lowest-mass trials
#              had gamma pinned at its 0.8 ceiling (optuna_fb.py:246-248).
#
# This probe does NOT change the model. It reads out the size distribution so
# the cause stops being a guess. Expect, per region, one line:
#   target=... configs=N sum_pip=X E_size=X mass_on_size1=F top_sizes=k:share;...
# E_size is computed independently of the PIP vector and MUST equal sum_pip; a
# mismatch is a bug, not a result. Read mass_on_size1 and top_sizes:
#   mass on size 1-3  -> the memo is sparse, and the FDR floor is NOT size-driven
#   mass on size ~12  -> scenario 7's sum(PIP)=11.8 is the sampler's |gamma|,
#                        and gamma/sparse_concrete are the levers
#
# Cost: one scenario, completion off. Job 4580 measured 877-2517 s per scenario
# at REGIONS=10 on SCENARIO=7, so the 2 h default has real margin.
#
#   bash scripts/tuning/size_probe.sh
#   STRATUM=sparse_inf SCENARIO=3 bash scripts/tuning/size_probe.sh
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

# SCENARIO=7 by default: scenario order is S-fastest, so 7 is S=1, phi=0.6 -
# ONE true causal with a strong signal, the cell whose sum(PIP) was 11.8 and
# whose sparse AP was the worst (0.395). It is the cleanest possible read.
SCENARIO="${SCENARIO:-7}"
REGIONS="${REGIONS:-10}"
MAX_ITER="${MAX_ITER:-2000}"
METHOD="${METHOD:-fb_xregion}"
R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
PY_MODULE="${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
FB_PY="${FB_PY:-$HOME/tools/py-venv-runner.sh}"
WALL="${WALL:-02:00:00}"
SELECT="${SELECT:-1:ncpus=1:mem=16gb}"
QUEUE="${QUEUE:-}"

PARAMS="${PARAMS:---sigma_sq 0.0899 --n_caus 10 --sparse_concrete 51 --prior_regularisation 0.1255 --lambda_l1 0.1223 --hierarchy_M 10.15 --gamma 0.1 --gamma_key 0.2 --prior_mass 3}"

J="$(mktemp -t fbsize_XXXXXX.sh)"
cat > "$J" <<PBS_EOF
#!/bin/bash
#PBS -N fbsize
$(if [[ -n "$QUEUE" ]]; then echo "#PBS -q $QUEUE"; else echo "# no queue named: the cx routing queue places this job"; fi)
#PBS -l select=${SELECT}
#PBS -l walltime=${WALL}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -uo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
module load ${PY_MODULE}
# Baked in, not exported from the submitting shell: PBS Pro here has no -V, so
# an export at submit time would leave the compute node on the defaults.
export FB_SIZE_LOG="${LOG_DIR}/size_diag.\${PBS_JOBID}.txt"
: > "\$FB_SIZE_LOG"
echo "[size on \$(hostname)] start \$(date)"
echo "sim      : ${SIM}"
echo "stratum  : ${STRATUM}   scenario ${SCENARIO}   regions ${REGIONS}"
echo "params   : ${PARAMS}"
echo "git HEAD : \$(git rev-parse --short HEAD)"
echo "size log : \$FB_SIZE_LOG"
echo
Rscript scripts/tuning/fb_objective.R \\
    --sim "${SIM}" --stratum "${STRATUM}" --method "${METHOD}" \\
    --scenario ${SCENARIO} --regions ${REGIONS} \\
    --beatrice_dir "${PROJECT_ROOT}/BEATRICE_annot_sparse" \\
    --python "${FB_PY}" --max_iter ${MAX_ITER} \\
    ${PARAMS} \\
    --posterior_completion 0
echo "--- Rscript exit status: \$? ---"
echo
echo "=== configuration-size distribution, one line per region ==="
cat "\$FB_SIZE_LOG" 2>/dev/null || echo "(no size diagnostic written)"
echo
echo "=== summary ==="
awk -F'[= ]' '/E_size=/ {
      for(i=1;i<=NF;i++){
        if(\$i=="sum_pip") s+=\$(i+1);
        if(\$i=="E_size")  e+=\$(i+1);
        if(\$i=="mass_on_size1") m+=\$(i+1);
      } n++
    } END {
      if(n) printf "regions=%d  mean sum_pip=%.4f  mean E_size=%.4f  mean mass_on_size1=%.4f\n", n, s/n, e/n, m/n;
      else print "no size lines parsed"
    }' "\$FB_SIZE_LOG" 2>/dev/null
echo "[size] done \$(date)"
PBS_EOF

bash -n "$J" || { echo "generated job script has a syntax error: $J" >&2; exit 1; }
echo "Size-probe job script: $J"
echo
echo "stratum ${STRATUM}  scenario ${SCENARIO}  regions ${REGIONS}  completion OFF"
echo "Reading: mass_on_size1 and top_sizes decide whether the FDR floor is"
echo "configuration-size driven. E_size must equal sum_pip exactly."
echo
JOBID="$(qsub "$J")"
echo "Submitted: $JOBID"
echo
echo "When it finishes:"
echo "  cat ${LOG_DIR}/fbsize.o\${JOBID%%.*} 2>/dev/null || cat ${LOG_DIR}/${JOBID%%.*}.OU"
