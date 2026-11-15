#!/bin/bash
# Read everything the v5 round produced, in one command, on a login node.
#
#   bash scripts/tuning/report_v5.sh
#
# Three reads, because they answer different questions:
#   fb_cand01            the DEPLOYED config (shipped + prior_mass=3) scored on
#                        all 27 cells. This is the number a new config has to
#                        beat. Read with read_candidate.py, because the trial
#                        violates the default mass window and so never appears
#                        in --report's Pareto front.
#   fb_retune_sparse_v5  the sparse front, searched under MASS_HI=3.0
#   fb_retune_sinf_v5    the sparse_inf front, searched with the mass
#                        constraint disabled
# The two fronts are EXPLORATION ONLY. One configuration has to serve every
# setup, so nothing here gets deployed until its candidate is re-scored on all
# 27 cells (CANARY=1 STRATA="sparse sparse_inf" FIX="...").
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

TD="$PWD/tuning"
PY_BIN="${PY_BIN:-$HOME/tools/fmpy-venv/bin/python}"
# The venv python is dynamically linked against the module's libpython, so the
# module has to be loaded before it will start (this has bitten before).
module load "${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}" 2>/dev/null || true

S_SPARSE="sparse=$TD/sim_sparse_binary_e5.4_n5000.rds"
S_INF02="sparse_inf=$TD/sim_sparseinf_rho0.2_binary_e5.4_n5000.rds"
S_INF06="sparse_inf=$TD/sim_sparseinf_rho0.6_binary_e5.4_n5000.rds"

hdr () { printf '\n\n########## %s ##########\n' "$1"; }

hdr "fb_cand01 - deployed config (prior_mass=3) on all 27 cells"
# A CANARY writes "${STUDY}_canary.journal", not "${STUDY}.journal" - the
# submitter's "Storage" banner line shows the search study's path even in
# canary mode, which sent an earlier read of this to the wrong file.
CAND_J=""
for c in "$TD/fb_cand01_canary.journal" "$TD/fb_cand01.journal"; do
  [[ -f "$c" ]] && CAND_J="$c" && break
done
if [[ -n "$CAND_J" ]]; then
  echo "(reading $CAND_J)"
  "$PY_BIN" scripts/tuning/read_candidate.py "$CAND_J"
else
  echo "NO canary journal. Looked for:"
  echo "  $TD/fb_cand01_canary.journal"
  echo "  $TD/fb_cand01.journal"
  echo "Journals that DO exist:"; ls -la "$TD"/*.journal 2>/dev/null
  # Pick the canary's own log, not merely the newest one (a search worker's log
  # is newer and says nothing about the canary).
  echo "Most recent canary log:"
  for ext in ER OU; do
    f="$(grep -l canary tuning/logs/*."$ext" 2>/dev/null | xargs -r ls -t 2>/dev/null | head -1)"
    [[ -n "$f" ]] && { echo "--- $f ---"; tail -30 "$f"; }
  done
fi

hdr "fb_retune_sparse_v5 - sparse front (searched under MASS_HI=3.0)"
"$PY_BIN" scripts/tuning/optuna_fb.py --sim "$S_SPARSE" \
  --study fb_retune_sparse_v5 --storage "$TD/fb_retune_sparse_v5.journal" \
  --objective cells --mass-hi 3.0 --report

hdr "fb_retune_sinf_v5 - sparse_inf front (mass constraint disabled)"
"$PY_BIN" scripts/tuning/optuna_fb.py --sim "$S_INF02" --sim "$S_INF06" \
  --study fb_retune_sinf_v5 --storage "$TD/fb_retune_sinf_v5.journal" \
  --objective cells --mass-lo -1000 --mass-hi 1000 --report
