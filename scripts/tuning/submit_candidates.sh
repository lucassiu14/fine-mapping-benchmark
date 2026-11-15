#!/bin/bash
# Re-score the shortlisted configurations on ALL 27 cells, one canary each.
#
#   bash scripts/tuning/submit_candidates.sh          # submit them all
#   bash scripts/tuning/submit_candidates.sh c071     # just one
#   bash scripts/tuning/report_candidates.sh          # read them afterwards
#
# WHY this step exists. The v5 studies searched ONE stratum each, so a front
# point's numbers describe that stratum only: c013/c071/c119 have never been
# scored on `sparse` at all. One configuration has to serve every setup
# (feedback 2026-11-14), so nothing here is deployable until it has been
# measured on the whole grid. A canary pins every parameter, so each of these
# is a measurement, not a search.
#
# The baseline to beat is the DEPLOYED configuration, already measured in
# fb_cand01_canary.journal:
#   sparse      AP 0.690  FDR@.95 0.336  ece 0.231  mass 2.28
#   sparse_inf  AP 0.444  FDR@.95 0.752  ece 0.557  mass 4.43
#   min-AP cell 0.246   worst-FDR cell 0.870 (sparse_inf|S=1)
#
# These five span the accuracy<->calibration frontier the sparse_inf front
# exposed: roughly 0.06 of worst-cell AP buys 0.22 of worst-cell FDR. Which end
# of that is wanted is a scientific call, which is why the whole span is being
# measured rather than one pick.
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

# Re-score each shortlisted configuration on ALL 27 cells.
# Generated - do not hand-edit the numbers.
declare -a CANDS=(
  # c177: trial #177 - best sparse compromise: minAP 0.630, worst-cell FDR 0.608, ece 0.217, mass 1.80
  "c177|sigma_sq=0.05837500572 n_caus=8 sparse_concrete=23 prior_regularisation=0.66209811 lambda_l1=0.3475319643 hierarchy_M=29.21091121 gamma=0.5899192553 gamma_key=0.05311019606 prior_mass=1.089373315"
  # c047: trial #47 - lowest ece on the sparse front: minAP 0.605, FDR 0.636, ece 0.203
  "c047|sigma_sq=0.04727297319 n_caus=7 sparse_concrete=52 prior_regularisation=43.79562756 lambda_l1=0.0001080049698 hierarchy_M=5.356595133 gamma=0.2243637687 gamma_key=0.0712233959 prior_mass=1.194177719"
  # c013: trial #13 - modest sinf calibration gain: minAP 0.226, worst-cell FDR 0.800
  "c013|sigma_sq=0.04913522991 n_caus=9 sparse_concrete=63 prior_regularisation=0.3545711875 lambda_l1=0.0001159125522 hierarchy_M=96.15779983 gamma=0.5641093937 gamma_key=0.05202014161 prior_mass=2.506077188"
  # c119: trial #119 - further along the trade: minAP 0.195, worst-cell FDR 0.733
  "c119|sigma_sq=0.2851370545 n_caus=9 sparse_concrete=55 prior_regularisation=0.3545711875 lambda_l1=0.274733521 hierarchy_M=96.15779983 gamma=0.5641093937 gamma_key=0.05202014161 prior_mass=2.506077188"
  # c071: trial #71 - best calibration found anywhere: minAP 0.183, worst-cell FDR 0.650, ece 0.382
  "c071|sigma_sq=0.05597826499 n_caus=8 sparse_concrete=26 prior_regularisation=1.541194772 lambda_l1=0.2942645457 hierarchy_M=29.61128793 gamma=0.72773119 gamma_key=0.5742523269 prior_mass=1.160827725"
)

WANT="${1:-}"
for entry in "${CANDS[@]}"; do
  tag="${entry%%|*}"; fix="${entry#*|}"
  if [[ -n "$WANT" && "$WANT" != "$tag" ]]; then continue; fi
  echo
  echo "=================== $tag ==================="
  # PM_MAX must stay above 0: the search space forces prior_mass to 0.0 when
  # PM_MAX <= 0, BEFORE --fix is consulted, so a pinned prior_mass is ignored.
  CANARY=1 STRATA="sparse sparse_inf" STUDY="fb_$tag" PM_MAX=10 \
    FIX="$fix" bash scripts/tuning/submit_fb_retune.sh
done
