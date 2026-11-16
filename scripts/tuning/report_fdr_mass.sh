#!/bin/bash
# FDR@0.95 against posterior mass, over EVERY completed trial in every journal.
# Same module-load + venv dance as report_candidates.sh: the venv python is
# dynamically linked against the module's libpython, so the load must precede it.
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
TD="$PWD/tuning"
PY_BIN="${PY_BIN:-$HOME/tools/fmpy-venv/bin/python}"
module load "${PY_MODULE:-Python/3.12.3-GCCcore-13.3.0}" 2>/dev/null || true
shopt -s nullglob
J=("$TD"/*.journal)
if [[ ${#J[@]} -eq 0 ]]; then
  echo "No journals under ${TD}"
  exit 0
fi
echo "Journals: ${#J[@]}"
printf '  %s\n' "${J[@]##*/}"
"$PY_BIN" scripts/tuning/fdr_vs_mass.py --min-nsel "${MIN_NSEL:-20}" "${J[@]}"
