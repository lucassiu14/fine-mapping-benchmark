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
# Completion's search window. Probe 4465: every pass stops by EXHAUSTING its
# neighbourhood (last_gain=0.0 at every width), so width is the lever and the
# round cap is not - complete_posterior already refreshes seeds each round.
#
# 128, not 512. The widths traded as mass_recovered 0.234 at 1.4s versus 0.252
# at 5.9s, i.e. 128 buys ~93% of the recovery for ~a quarter of the cost. The
# earlier default of 512 was set from that 5.9s figure, which was measured on
# 2-region probes whose memos held ~8,000 configs. In the real 10-region run
# the memo reaches 72,000+ and completion's cost scales with memo size, so
# canary 4580 ran 877-2517s per scenario against the OFF arm's ~390s - 2.7x
# overall, +280 min on the trial - and died on the 8h walltime at scenario
# 26/27. Hence also CANARY_WALLTIME below.
CWIDTH="${CWIDTH:-128}"
# submit_fb_retune.sh defaults this to 08:00:00, which 4580 exceeded. The OFF
# arm needs ~3h and is unaffected by the larger ceiling.
export CANARY_WALLTIME="${CANARY_WALLTIME:-16:00:00}" 

# submit_fb_retune.sh's default FIX (:98), verbatim, plus the anchor. Anything
# added here must be added to BOTH arms - which is now automatic, since both
# arms are built from this one string.
BASE="sigma_sq=0.0899 n_caus=10 sparse_concrete=51 prior_regularisation=0.1255"
BASE="${BASE} lambda_l1=0.1223 hierarchy_M=10.15 gamma=0.1 gamma_key=0.2"
BASE="${BASE} prior_mass=${PRIOR_MASS}"

echo "== 1/5  git pull =============================================="
git pull --ff-only || { echo "pull failed"; exit 1; }

echo
echo "== 2/5  module load + reinstall fmbenchmark ==================="
module load "${R_MODULE}" 2>/dev/null || true
command -v Rscript >/dev/null 2>&1 || {
  echo "ERROR: Rscript not on PATH. Try: module load ${R_MODULE}"; exit 1; }
Rscript -e 'install.packages(".", repos = NULL, type = "source")' \
  2>&1 | tail -3 || { echo "install failed"; exit 1; }

echo
echo "== 3/5  verify the flag reached the installed wrapper ========="
OK="$(Rscript -e 'cat(any(grepl("posterior_completion", deparse(fmbenchmark:::.fb_joint_scenario_setup))))' 2>/dev/null | tail -1)"
if [[ "$OK" != "TRUE" ]]; then
  echo "ERROR: the installed fmbenchmark does NOT carry --posterior_completion."
  echo "       The canary would silently run the old wrapper and measure"
  echo "       nothing. Fix the install before submitting."
  exit 1
fi
echo "  installed wrapper carries --posterior_completion: TRUE"

echo
echo "== 4/5  preflight receipt ====================================="
# A canary costs hours. It is not submitted unless ONE scenario has already
# run clean at the PRODUCTION region count, at this width, on this exact
# commit - and unless the walltime provably fits the cost that run measured.
# 956 (2h, a z-shape bug a 2-D test input could not catch) and 4580 (7.6h,
# killed at scenario 26/27 by a cost estimate that was 10x low) are both
# failures this check would have stopped.
COMMIT="$(git rev-parse HEAD 2>/dev/null)"
RECEIPT="${TUNE_DIR:-$PWD/tuning}/preflight/${COMMIT}.receipt"
if [[ "${SKIP_PREFLIGHT:-0}" == "1" ]]; then
  echo "  SKIP_PREFLIGHT=1 - submitting unverified, on your head."
elif [[ ! -f "$RECEIPT" ]]; then
  echo "  REFUSING: no preflight receipt for commit ${COMMIT:0:12}."
  echo
  echo "  Run this first (~15-45 min, one scenario at 10 regions, both arms):"
  echo "    CWIDTH=${CWIDTH} bash scripts/tuning/preflight_canary.sh"
  echo
  echo "  It writes the receipt only if both arms return n_failed=0 with a"
  echo "  real AP, and the receipt is what sizes the walltime below."
  echo "  SKIP_PREFLIGHT=1 overrides, which is how hours get wasted."
  exit 1
else
  R_WIDTH="$(sed -n 's/^width=//p' "$RECEIPT" | tail -1)"
  R_SECS="$(sed -n 's/^max_scenario_secs=//p' "$RECEIPT" | tail -1)"
  R_REG="$(sed -n 's/^regions=//p' "$RECEIPT" | tail -1)"
  echo "  receipt: commit ${COMMIT:0:12}  width=${R_WIDTH}  regions=${R_REG}"
  echo "           slowest scenario measured: ${R_SECS}s"
  if [[ "$R_WIDTH" != "$CWIDTH" ]]; then
    echo "  REFUSING: receipt is for width ${R_WIDTH}, this submits ${CWIDTH}."
    echo "            Cost scales with the window, so the receipt does not"
    echo "            vouch for this configuration. Re-run preflight with"
    echo "            CWIDTH=${CWIDTH}."
    exit 1
  fi
  # 27 scenarios, x1.4 margin. The receipt times a deliberately SLOW scenario,
  # so this is an over-estimate for the average one - which is the right
  # direction to err.
  NEED_H="$(python3 -c "
import math
print(max(4, math.ceil(27 * float('${R_SECS}') * 1.4 / 3600)))" 2>/dev/null)"
  HAVE_H="${CANARY_WALLTIME%%:*}"
  echo "           walltime needed ~${NEED_H}h, requested ${HAVE_H}h"
  if (( NEED_H > HAVE_H )); then
    export CANARY_WALLTIME="$(printf '%02d:00:00' "$NEED_H")"
    echo "  raised CANARY_WALLTIME to ${CANARY_WALLTIME} to fit the measured cost."
  fi
fi

echo
echo "== 5/5  submit ================================================"
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
