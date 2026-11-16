#!/bin/bash
# Read every candidate canary on one screen, worst-cell first.
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
TD="$PWD/tuning"
PY_BIN="${PY_BIN:-$HOME/tools/fmpy-venv/bin/python}"
module load "${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}" 2>/dev/null || true
shopt -s nullglob
# Match EVERY canary journal. The old pattern was
#   fb_cand01_canary.journal fb_c[0-9][0-9][0-9]_canary.journal
# which silently omitted fb_complete01_canary.journal - the run actually
# under test - and printed only the stale candidates. The table then looks
# complete while missing the one row that matters, which is worse than an
# empty table. Any new STUDY name is covered by the wildcard.
J=("$TD"/fb_*_canary.journal)
if [[ ${#J[@]} -eq 0 ]]; then
  echo "No candidate journals yet. Expected ${TD}/fb_c*_canary.journal"
  ls -la "$TD"/*.journal 2>/dev/null
  exit 0
fi
"$PY_BIN" scripts/tuning/read_candidate.py "${J[@]}"
