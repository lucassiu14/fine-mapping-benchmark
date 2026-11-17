#!/bin/bash
# One scenario at the PRODUCTION region count, both arms, then a receipt.
#
# Why this exists. Every hour wasted on this work was wasted the same way: code
# validated against something hand-built, then submitted to a multi-hour job.
#   - canary 956  (2h)   died on a 1-D vs 2-D z shape. The unit test built z
#                        2-D, so it could not fail.
#   - canary 4580 (7.6h) died on the 8h walltime at scenario 26/27, because
#                        width-512 cost was estimated from a 2-region probe on
#                        an already-saturated memo (0.2s) and was wrong by ~10x:
#                        real cost was 877-2517s per scenario against ~390s.
# debug_one_scenario.sh runs 2 regions, which catches crashes but CANNOT
# measure cost: the joint trainer holds every region's memo at once, so cost
# grows non-linearly in the region count. This runs ONE scenario at REGIONS
# (default 10, the production value), which both catches the crash and measures
# the real per-scenario seconds.
#
# It writes tuning/preflight/<commit>.receipt. submit_completion_canary.sh
# refuses to submit without a receipt for the CURRENT HEAD at the width it is
# about to use, and sizes the walltime from it.
#
#   bash scripts/tuning/preflight_canary.sh              # ~15-30 min
#   CWIDTH=512 bash scripts/tuning/preflight_canary.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

PROJECT_ROOT="$(pwd)"
TUNE_DIR="${TUNE_DIR:-${PROJECT_ROOT}/tuning}"
LOG_DIR="${FMB_LOG_DIR:-${TUNE_DIR}/logs}"
RECEIPT_DIR="${TUNE_DIR}/preflight"
mkdir -p "$LOG_DIR" "$RECEIPT_DIR"

COMMIT="$(git rev-parse HEAD 2>/dev/null)" || { echo "not a git checkout" >&2; exit 1; }
if [[ -n "$(git status --porcelain -- BEATRICE_annot_sparse scripts R 2>/dev/null)" ]]; then
  echo "WARNING: uncommitted changes under BEATRICE_annot_sparse/scripts/R." >&2
  echo "         The receipt is keyed by commit, so it would vouch for code" >&2
  echo "         that is not what runs. Commit first." >&2
  exit 1
fi

STRATUM="${STRATUM:-sparse_inf}"     # the harsher stratum: biggest memos
case "$STRATUM" in
  sparse)     SIM="${SIM:-${TUNE_DIR}/sim_sparse_binary_e5.4_n5000.rds}" ;;
  sparse_inf) SIM="${SIM:-${TUNE_DIR}/sim_sparseinf_rho0.2_binary_e5.4_n5000.rds}" ;;
  *) echo "STRATUM must be sparse or sparse_inf" >&2; exit 1 ;;
esac
[[ -f "$SIM" ]] || { echo "sim not found: $SIM" >&2; exit 1; }

# Scenario 7 by default: in canary 4580 the high-mass scenarios were the slow
# ones (mass 22.7 -> 2522s), so preflight must time a SLOW scenario, not an
# average one, or the walltime it certifies is optimistic.
SCENARIO="${SCENARIO:-7}"
REGIONS="${REGIONS:-10}"
CWIDTH="${CWIDTH:-128}"
TOL="${TOL:-0.001}"
MAX_ITER="${MAX_ITER:-2000}"
METHOD="${METHOD:-fb_xregion}"
R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
PY_MODULE="${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
FB_PY="${FB_PY:-$HOME/tools/py-venv-runner.sh}"
WALL="${WALL:-03:00:00}"
SELECT="${SELECT:-1:ncpus=1:mem=16gb}"
PARAMS="${PARAMS:---sigma_sq 0.0899 --n_caus 10 --sparse_concrete 51 --prior_regularisation 0.1255 --lambda_l1 0.1223 --hierarchy_M 10.15 --gamma 0.1 --gamma_key 0.2 --prior_mass 3}"

RECEIPT="${RECEIPT_DIR}/${COMMIT}.receipt"
J="$(mktemp -t fbpre_XXXXXX.sh)"
cat > "$J" <<PBS_EOF
#!/bin/bash
#PBS -N fbpreflt
#PBS -l select=${SELECT}
#PBS -l walltime=${WALL}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -uo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
module load ${PY_MODULE}
echo "[preflight on \$(hostname)] start \$(date)"
echo "commit   : ${COMMIT}"
echo "stratum  : ${STRATUM}  scenario ${SCENARIO}  regions ${REGIONS}  width ${CWIDTH}"
OK=0
MAXS=0
for TOLV in ${TOL} 0; do
  echo "=== posterior_completion=\${TOLV} ==="
  OUT="\$(Rscript scripts/tuning/fb_objective.R \\
      --sim "${SIM}" --stratum "${STRATUM}" --method "${METHOD}" \\
      --scenario ${SCENARIO} --regions ${REGIONS} \\
      --beatrice_dir "${PROJECT_ROOT}/BEATRICE_annot_sparse" \\
      --python "${FB_PY}" --max_iter ${MAX_ITER} \\
      ${PARAMS} \\
      --posterior_completion \${TOLV} \\
      --completion_seeds 8 --completion_width ${CWIDTH} \\
      --completion_rounds 4 2>&1)"
  echo "\$OUT"
  NF="\$(printf '%s\\n' "\$OUT" | sed -n 's/^n_failed=//p' | tail -1)"
  AP="\$(printf '%s\\n' "\$OUT" | sed -n 's/^ap=//p' | tail -1)"
  SEC="\$(printf '%s\\n' "\$OUT" | sed -n 's/^seconds=//p' | tail -1)"
  echo "--> parsed n_failed=\${NF:-?} ap=\${AP:-?} seconds=\${SEC:-?}"
  if [[ "\${NF:-1}" == "0" && -n "\${AP:-}" && "\${AP}" != "NA" && "\${AP}" != "NaN" ]]; then
    OK=\$((OK+1))
    MAXS="\$(python3 -c "print(max(float('\${MAXS}'), float('\${SEC:-0}')))")"
  else
    echo "ARM FAILED (completion=\${TOLV}): n_failed=\${NF:-?} ap=\${AP:-?}"
  fi
done
if [[ "\$OK" == "2" ]]; then
  printf 'commit=%s\nwidth=%s\nregions=%s\nstratum=%s\nscenario=%s\nmax_scenario_secs=%s\n' \\
    "${COMMIT}" "${CWIDTH}" "${REGIONS}" "${STRATUM}" "${SCENARIO}" "\$MAXS" > "${RECEIPT}"
  echo "PREFLIGHT PASS - receipt written: ${RECEIPT}"
  cat "${RECEIPT}"
else
  rm -f "${RECEIPT}"
  echo "PREFLIGHT FAIL (\$OK/2 arms ok) - no receipt. Do NOT submit a canary."
fi
echo "[preflight] done \$(date)"
PBS_EOF

bash -n "$J" || { echo "generated script has a syntax error: $J" >&2; exit 1; }
echo "Preflight script: $J"
JOBID="$(qsub "$J")" || { echo "qsub failed" >&2; exit 1; }
echo "Preflight job: ${JOBID}"
echo
echo "One scenario at ${REGIONS} regions (the production count) and width ${CWIDTH},"
echo "both arms, on the SLOW scenario ${SCENARIO}. ~15-45 min. Read it with:"
echo "  grep -E 'parsed|PREFLIGHT|max_scenario_secs' ${LOG_DIR}/${JOBID}.OU"
echo
echo "On PASS it writes ${RECEIPT}, which submit_completion_canary.sh requires"
echo "and uses to size the walltime. On FAIL no receipt is written and the"
echo "canary submitter will refuse."
