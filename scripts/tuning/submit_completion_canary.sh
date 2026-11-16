#!/bin/bash
# Posterior-completion canary: the deployed configuration, completion ON.
#
# Pairs against fb_cand01_canary, which is the SAME configuration with
# completion off (sparse AP 0.690 / FDR@0.95 0.336; sparse_inf 0.444 / 0.752),
# so the difference is attributable to completing the sum and nothing else.
#
# Safe to run from a fresh login session: loads R, reinstalls fmbenchmark (the
# objective prefers the INSTALLED package over the checkout, so a bare git pull
# leaves the old wrapper in place and the new flag never reaches the model),
# verifies the flag actually made it in, and only then submits.
#
#   bash scripts/tuning/submit_completion_canary.sh
#
# Env overrides: STUDY (default fb_complete01), TOL (default 0.001),
# R_MODULE, and anything submit_fb_retune.sh itself accepts.
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
STUDY="${STUDY:-fb_complete01}"
TOL="${TOL:-0.001}"

# The shipped configuration, unchanged, plus completion. Keep this identical to
# submit_fb_retune.sh's default FIX or the comparison stops being a comparison.
BASE="sigma_sq=0.0899 n_caus=10 sparse_concrete=51 prior_regularisation=0.1255"
BASE="${BASE} lambda_l1=0.1223 hierarchy_M=10.15 gamma=0.1 gamma_key=0.2"
BASE="${BASE} prior_mass=3"

echo "== 1/4  git pull =============================================="
git pull --ff-only || { echo "pull failed"; exit 1; }

echo
echo "== 2/4  module load + reinstall fmbenchmark ==================="
module load "${R_MODULE}" 2>/dev/null || true
command -v Rscript >/dev/null 2>&1 || {
  echo "ERROR: Rscript not on PATH. Try: module load ${R_MODULE}"; exit 1; }
Rscript -e 'install.packages(".", repos = NULL, type = "source")' \
  2>&1 | tail -3 || { echo "install failed"; exit 1; }

echo
echo "== 3/4  verify the flag reached the installed wrapper ========="
OK="$(Rscript -e 'cat(any(grepl("posterior_completion", deparse(fmbenchmark:::.fb_joint_scenario_setup))))' 2>/dev/null | tail -1)"
if [[ "$OK" != "TRUE" ]]; then
  echo "ERROR: the installed fmbenchmark does NOT carry --posterior_completion."
  echo "       The canary would silently run the old wrapper and measure"
  echo "       nothing. Fix the install before submitting."
  exit 1
fi
echo "  installed wrapper carries --posterior_completion: TRUE"

echo
echo "== 4/4  submit ================================================"
CANARY=1 STRATA="sparse sparse_inf" STUDY="${STUDY}" PM_MAX=10 \
  FIX="${BASE} posterior_completion=${TOL}" \
  bash scripts/tuning/submit_fb_retune.sh
rc=$?

echo
echo "Expect 'Pinned configuration:' to end with posterior_completion=${TOL}."
echo "It prints AFTER qsub, so if qsub fails you will not see it - that is the"
echo "submission failing, not the pin."
echo
echo "When it finishes (~3h), read it against the completion-off baseline:"
echo "  bash scripts/tuning/report_candidates.sh"
echo "  PY_BIN=\${PY_BIN:-\$HOME/tools/fmpy-venv/bin/python} \\"
echo "    \$PY_BIN scripts/tuning/read_candidate.py \\"
echo "      tuning/fb_cand01_canary.journal tuning/${STUDY}_canary.journal"
exit $rc
