#!/bin/bash
# =============================================================================
# scripts/tuning/submit_fb_mass_line.sh
#
# prior_mass LINE SEARCH on the sparse front points, 2026-11-19.
#
# WHY. The 27-cell re-scoring (jobs 3030033-3030045) found ALL EIGHT sparse
# front points mass-infeasible on sparse_inf (3.12-4.38 against a 2.0 ceiling),
# while the only feasible candidate, i220, regresses sparse by 0.06 AP and 0.11
# on the worst cell. prior_mass is the parameter that EXISTS to control this:
# the configuration prior is independent-Bernoulli, so sum_j PIP_j =
# E[|gamma| | z] and is governed by sum_j p0_j = M. Nothing else in the model
# sets the posterior's total scale. So before accepting i220's sparse
# regression or abandoning the constraint, ask the one-parameter question:
# can a sparse-optimal config be made feasible on sparse_inf by lowering M?
#
# THE ARITHMETIC FIRST, because it decides which configs are worth probing.
# M scales BOTH strata together, so feasibility on both requires
#       ratio = mass_sparse_inf / mass_sparse  <  2.0 / 0.5 = 4
# and the measured ratios are:
#   s430 3.10   s270 3.72   s436 3.85   s339 3.93 | s460 4.22  s400 4.34
#   s389 4.41   s276 4.98                          (i220 2.18, for reference)
# The four on the left can be brought inside the window; the four on the right
# CANNOT, at any M, under exact proportionality. s400 - the best sparse config -
# is one of the impossible ones, so the probe cannot simply rescale the winner.
# s270 is excluded separately: its sparse_inf AP (0.385) has already regressed
# below the baseline's 0.418, which no change of M repairs.
#
# M* below is the M that puts sparse_inf mass exactly at the 2.0 ceiling, and
# mSp* the sparse mass that proportionality then predicts:
#   s430  M 2.420 -> M* 1.553, sparse mass 1.005 -> 0.645   (comfortable)
#   s339  M 2.667 -> M* 1.417, sparse mass 0.959 -> 0.509   (marginal)
#   s436  M 1.171 -> M* 0.704, sparse mass 0.864 -> 0.520   (marginal)
#   s400  M 1.910 -> M* 0.974, sparse mass 0.903 -> 0.461   (fails the floor)
#
# WHAT THIS CANNOT PREDICT, and is therefore the point of running it: M changes
# AP as well as mass. The anchor was adopted BECAUSE M=3 beat M=1 on sparse
# locally (mean AP 0.773 -> 0.817, run_benchmark_job.R:283-290), so lowering it
# may give back exactly the sparse AP the exercise is trying to protect. It
# should also pull sparse_inf FDR down, since i220 at M=0.648 sits at 0.503
# against these points' 0.671-0.745 at M=1.2-2.7. Neither effect is derivable.
#
# FALSIFIABLE PREDICTION. If the response is proportional: s430 at M=1.3-1.55
# comes back feasible and wins the selection rule, and s400 at M=0.95 comes back
# INFEASIBLE on the sparse floor. If s400m095 is feasible, the sparse_inf
# response is sub-linear in M and every config above has more room than the
# table says. If s430m155 is NOT feasible, it is super-linear and the whole
# line search fails - in which case the answer is i220 or dropping the
# constraint, not more probing.
#
#   bash scripts/tuning/submit_fb_mass_line.sh
#   DRY=1 bash scripts/tuning/submit_fb_mass_line.sh
#   ONLY="s430m155" bash scripts/tuning/submit_fb_mass_line.sh
# =============================================================================
set -euo pipefail
[[ -f scripts/tuning/optuna_fb.py ]] || { echo "ERROR: run from the project root" >&2; exit 1; }

DRY="${DRY:-0}"
ONLY="${ONLY:-}"

# tag source lambda_l1 prior_regularisation sigma_sq hierarchy_M n_caus sparse_concrete prior_mass gamma gamma_key
# Every column except prior_mass is the front point verbatim; prior_mass is the
# one thing being moved. Base points are s430 (most headroom and the lowest
# sparse_inf FDR of the eight), s339 and s436 (marginal), and s400 (the
# proportionality control, predicted to fail).
read -r -d '' CANDIDATES <<'TBL' || true
s430m130 s430_M2.420to1.30 0.7788834306113981 1.0050853288840873 0.09393351811229461 69.32383733527517 8 115 1.30 0.37581390099410567 0.10413246114459562
s430m155 s430_M2.420to1.55 0.7788834306113981 1.0050853288840873 0.09393351811229461 69.32383733527517 8 115 1.55 0.37581390099410567 0.10413246114459562
s339m120 s339_M2.667to1.20 0.5577273229404637 1.8856811642400804 0.1032473351358484 1.269432573354525 9 30 1.20 0.16239941433566002 0.09273161704121474
s436m070 s436_M1.171to0.70 0.21053958715166374 0.34049473982466116 0.24922956104339353 1.269432573354525 9 41 0.70 0.07208163585780197 0.4260666727023295
s400m095 s400_M1.910to0.95 0.28362249427122 0.8038180994579655 0.07612524744363866 31.7543967020142 9 20 0.95 0.06249756436559903 0.09185713974666732
TBL

n=0
while read -r tag src l1 pr sq hM nc sc pm g gk; do
  [[ -z "${tag:-}" ]] && continue
  if [[ -n "$ONLY" ]] && [[ " $ONLY " != *" $tag "* ]]; then continue; fi
  STUDY="fb_sel_${tag}"
  FIX="sigma_sq=${sq} n_caus=${nc} sparse_concrete=${sc} \
prior_regularisation=${pr} lambda_l1=${l1} hierarchy_M=${hM} \
gamma=${g} gamma_key=${gk} prior_mass=${pm}"
  echo "=========================================================="
  echo "candidate ${tag}  (${src})"
  echo "  journal: tuning/${STUDY}_canary.journal"
  if [[ "$DRY" == "1" ]]; then echo "  FIX: ${FIX}"; n=$((n+1)); continue; fi
  CANARY=1 STRATA="sparse sparse_inf" STUDY="$STUDY" \
    PM_MIN=0.5 PM_MAX=10 SIG_MIN=0.0001 SIG_MAX=2.0 \
    FIX="$FIX" \
    bash scripts/tuning/submit_fb_retune.sh < /dev/null
  # ^ stdin redirected: the loop is fed by `done <<< "$CANDIDATES"`, and any
  # command inside that reads stdin eats the remaining candidate lines.
  n=$((n+1))
done <<< "$CANDIDATES"

echo "=========================================================="
echo "${n} candidate job(s) handled."
cat <<'MSG'

These journals are named fb_sel_* so select_fb_config.py picks them up next to
the originals and applies the same rule across all of them:

  module load Python/3.12.3-GCCcore-13.3.0
  $HOME/tools/fmpy-venv/bin/python scripts/tuning/select_fb_config.py
MSG
