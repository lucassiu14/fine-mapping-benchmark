#!/bin/bash
# Read every candidate canary on one screen, worst-cell first.
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
TD="$PWD/tuning"
PY_BIN="${PY_BIN:-$HOME/tools/fmpy-venv/bin/python}"
module load "${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}" 2>/dev/null || true
shopt -s nullglob
J=("$TD"/fb_cand01_canary.journal "$TD"/fb_c[0-9][0-9][0-9]_canary.journal)
if [[ ${#J[@]} -eq 0 ]]; then
  echo "No candidate journals yet. Expected ${TD}/fb_c*_canary.journal"
  ls -la "$TD"/*.journal 2>/dev/null
  exit 0
fi
"$PY_BIN" scripts/tuning/read_candidate.py "${J[@]}"
