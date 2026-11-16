#!/bin/bash
# Posterior-completion canary: a MATCHED PAIR at one configuration, completion
# ON in one arm and OFF in the other.
#
# It used to submit the ON arm alone and tell the reader to compare it with
# fb_cand01_canary. That comparison was invalid. fb_cand01_canary ran with the
# prior anchor DISABLED: optuna_fb.py:253 sets prior_mass to the literal 0.0
# whenever PM_MAX <= 0, which is why read_candidate shows it as neither a
# searched param nor a pinned one. This script pins prior_mass=3. Since
# sum(p0) "governs directly the posterior expected configuration size - hence
# sum_j PIP_j" (anchor_prior_mass's own docstring), the two runs differ in the
# prior scale as well as in completion, and every PIP, AP, FDR and mass moves
# with it. Two variables, one comparison.
#
# So both arms are now submitted here, from one BASE, differing ONLY in
# posterior_completion. Set PAIR=0 to submit the ON arm by itself.
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
OFF_STUDY="${OFF_STUDY:-fb_nocomplete01}"
TOL="${TOL:-0.001}"
PAIR="${PAIR:-1}"
# Held fixed across BOTH arms, so it cannot be confounded with completion. 3 is
# the provisionally deployed anchor value; 0 disables the anchor entirely.
PRIOR_MASS="${PRIOR_MASS:-3}"
# Completion's search window. Probe 4465 measured the trade directly: every
# pass stops by EXHAUSTING its neighbourhood (last_gain=0.0 at every width),
# so width is the lever and the round cap is not - complete_posterior already
# refreshes its seeds each round. Going 32 -> 512 raised the mass recovered
# from 0.165 to 0.252 in one region and 0.351 to 0.422 in another, at ~6s per
# region, about 27 min on a 175 min trial. Only the ON arm reads it.
CWIDTH="${CWIDTH:-512}"

# submit_fb_retune.sh's default FIX (:98), verbatim, plus the anchor. Anything
# added here must be added to BOTH arms - which is now automatic, since both
# arms are built from this one string.
BASE="sigma_sq=0.0899 n_caus=10 sparse_concrete=51 prior_regularisation=0.1255"
BASE="${BASE} lambda_l1=0.1223 hierarchy_M=10.15 gamma=0.1 gamma_key=0.2"
BASE="${BASE} prior_mass=${PRIOR_MASS}"

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
echo "Both arms are pinned at prior_mass=${PRIOR_MASS} and differ ONLY in"
echo "posterior_completion (the ON arm also sets completion_width=${CWIDTH},"
echo "which the OFF arm never reads because completion is disabled there),"
echo "so the difference between them is attributable to"
echo "completing the sum and to nothing else."
echo

echo "--- arm 1/2: completion ON (posterior_completion=${TOL}) ---"
CANARY=1 STRATA="sparse sparse_inf" STUDY="${STUDY}" PM_MAX=10 \
  FIX="${BASE} posterior_completion=${TOL} completion_width=${CWIDTH}" \
  bash scripts/tuning/submit_fb_retune.sh
rc=$?

if [[ "$PAIR" == "1" ]]; then
  echo
  echo "--- arm 2/2: completion OFF (posterior_completion=0) ---"
  # Explicitly 0 rather than omitted: 0 is the documented "keep the partial
  # sum" value, and passing it keeps the two argv vectors identical in length
  # and order, so nothing but the tolerance differs between the arms.
  CANARY=1 STRATA="sparse sparse_inf" STUDY="${OFF_STUDY}" PM_MAX=10 \
    FIX="${BASE} posterior_completion=0" \
    bash scripts/tuning/submit_fb_retune.sh
  rc2=$?
  [[ $rc -eq 0 ]] && rc=$rc2
fi

echo
echo "Expect each arm's 'Pinned configuration:' to end with its own"
echo "posterior_completion (${TOL} and 0). It prints AFTER qsub, so if qsub"
echo "fails you will not see it - that is the submission failing, not the pin."
echo
echo "The two arms run concurrently, so both land in ~3h. Read them together:"
echo "  bash scripts/tuning/report_candidates.sh"
echo
echo "A trial that reports values [0.0, 1.0, 1.0] produced NO PIPs at all and"
echo "is a failure, not a poor score - check the [WARN] lines for n_failed."
exit $rc
