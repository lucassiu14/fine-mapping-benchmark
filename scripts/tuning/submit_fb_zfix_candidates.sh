#!/bin/bash
# =============================================================================
# scripts/tuning/submit_fb_zfix_candidates.sh
#
# JOINT RE-SCORING of the post-z-fix Pareto fronts. One deployable
# configuration has to serve every setup, and neither search produced one:
#   tuning/fb_zfix_sparse01.journal  483 trials, scored on the 9 SPARSE cells only
#   tuning/fb_zfix_sinf01.journal    232 trials, scored on the 18 SPARSE_INF cells only
# so every front point is missing half the grid. This re-scores each candidate
# on ALL 27 cells (9 sparse + 9 x 2 rho sparse_inf), same sims, same seed, every
# parameter pinned - a cost-free comparison in the sense that nothing is being
# searched.
#
# SELECTION RULE (user's decision, 2026-11-19): maximise the WORST SPARSE CELL's
# AP, subject to
#   (a) mass feasible on BOTH strata, sum(PIP)/#causal in [0.5, 2.0], and
#   (b) no regression against the deployed config on sparse_inf.
# The deployed config's own numbers on this exact grid are already measured -
# tuning/fb_zfix01_canary.journal, read with read_candidate.py:
#   sparse      AP 0.792  FDR@.95 0.023  ece_hi 0.104  mass 0.923   worst cell AP 0.499
#   sparse_inf  AP 0.418  FDR@.95 0.744  ece_hi 0.565  mass 4.267   worst cell FDR 0.860
# so it is NOT re-run here. Note it FAILS (a) on sparse_inf at mass 4.267, which
# is why feasibility there is a real constraint and not a formality.
#
# WHY NOT just take the sparse front's best point: the two fronts disagree on
# sigma_sq and only on sigma_sq. Sparse settles in 0.026-0.249 (the deployed
# 0.0899 is inside it); sparse_inf settles in 0.32-1.75. Under sparse_inf
# sigma_sq stops being the prior variance and acts as a misspecification
# absorber - a larger prior variance lets one sparse effect stand in for diffuse
# mass the model has no component for - which is why it moves OPPOSITE to the
# oracle phi/(S(1-phi)), that being smaller under sparse_inf since only rho*phi
# of the genetic variance sits in the sparse component. How much a single
# sigma_sq costs is therefore a measurement, not a derivation.
#
#   bash scripts/tuning/submit_fb_zfix_candidates.sh          # submit all 10
#   DRY=1 bash scripts/tuning/submit_fb_zfix_candidates.sh    # print, submit nothing
#   ONLY="s400 i220" bash scripts/tuning/submit_fb_zfix_candidates.sh
# =============================================================================
set -euo pipefail
[[ -f scripts/tuning/optuna_fb.py ]] || { echo "ERROR: run from the project root" >&2; exit 1; }

DRY="${DRY:-0}"
ONLY="${ONLY:-}"

# tag  source                       the 9 pinned parameters
# Columns are lambda_l1 prior_regularisation sigma_sq hierarchy_M n_caus
# sparse_concrete prior_mass gamma gamma_key, read verbatim off each front point.
# s* = tuning/fb_zfix_sparse01 front (8 points, ALL at worst_fdr_xs 0.0000).
# i* = tuning/fb_zfix_sinf01 front, kept as a fallback for the case where no
#      sparse point clears constraint (b): i220 is its best min-cell AP with
#      feasible mass, i146 its lowest worst-cell FDR that has not collapsed AP.
read -r -d '' CANDIDATES <<'TBL' || true
s400 sparse#400_minAP0.6713 0.28362249427122 0.8038180994579655 0.07612524744363866 31.7543967020142 9 20 1.9099459067836553 0.06249756436559903 0.09185713974666732
s270 sparse#270_minAP0.6597 0.7788834306113981 0.605602232973392 0.025602259202236972 69.32383733527517 8 102 2.419689920785676 0.37581390099410567 0.09778742750581951
s339 sparse#339_minAP0.6470 0.5577273229404637 1.8856811642400804 0.1032473351358484 1.269432573354525 9 30 2.6672716896851467 0.16239941433566002 0.09273161704121474
s276 sparse#276_minAP0.6467 0.21053958715166374 1.8856811642400804 0.04318688161586183 1.269432573354525 9 30 1.1714569161071122 0.07208163585780197 0.12350347740807749
s430 sparse#430_minAP0.6327 0.7788834306113981 1.0050853288840873 0.09393351811229461 69.32383733527517 8 115 2.419689920785676 0.37581390099410567 0.10413246114459562
s460 sparse#460_minAP0.6283 0.21053958715166374 0.2865600794242421 0.24922956104339353 3.1184642070816286 6 41 2.627213681115266 0.07212086517631781 0.3986031479784162
s436 sparse#436_minAP0.6181 0.21053958715166374 0.34049473982466116 0.24922956104339353 1.269432573354525 9 41 1.1714569161071122 0.07208163585780197 0.4260666727023295
s389 sparse#389_minAP0.6056 0.30631910208919255 3.8899086232684494 0.11609985261154512 2.873121117789073 6 19 0.9253955812609035 0.03387261029286738 0.11664232466100681
i220 sinf#220_minAP0.2035 0.008163369572113525 44.923255939442164 0.3227016375881618 21.026350169323493 5 12 0.6475698183954328 0.2932567447107542 0.2905907897290682
i146 sinf#146_fdrxs0.5015 0.0024406410150025776 6.444340168909987 0.7381707222764375 13.36758583916022 3 196 2.2210229395325958 0.49819782110645583 0.03546382297529168
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
  if [[ "$DRY" == "1" ]]; then
    echo "  FIX: ${FIX}"
    n=$((n+1)); continue
  fi
  # PM_MAX must stay above 0 or the search space forces prior_mass=0.0 BEFORE
  # --fix is consulted, silently discarding the pinned anchor.
  CANARY=1 STRATA="sparse sparse_inf" STUDY="$STUDY" \
    PM_MIN=0.5 PM_MAX=10 SIG_MIN=0.0001 SIG_MAX=2.0 \
    FIX="$FIX" \
    bash scripts/tuning/submit_fb_retune.sh < /dev/null
  # ^ stdin redirected: this loop is fed by `done <<< "$CANDIDATES"`, and any
  # command inside it that reads stdin (qsub has been known to) consumes the
  # remaining candidate lines, so only the first job would ever be submitted.
  n=$((n+1))
done <<< "$CANDIDATES"

echo "=========================================================="
echo "${n} candidate job(s) handled."
cat <<'MSG'

Read them when they finish - NOT with optuna_fb.py --report, whose headline
block comes from study.best_trials and is EMPTY whenever the mass constraint
fails, which is the expected outcome for several of these:

  module load Python/3.12.3-GCCcore-13.3.0
  $HOME/tools/fmpy-venv/bin/python scripts/tuning/read_candidate.py \
      tuning/fb_sel_*_canary.journal
MSG
