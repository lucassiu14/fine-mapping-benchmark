#!/bin/bash
# Triage a run: which tasks failed, where, and why.
#
#   bash triage.sh            # summary across all logs
#   bash triage.sh 02_prep    # only that stage
#
# Every abnormal line in a stage log is prefixed FATAL or WARN by lib.sh, so
# this is just a structured grep - but it saves reading 265 files by hand.
set -uo pipefail
cd "${APP_ROOT:?source env.sh first}"
FILTER="${1:-}"
shopt -s nullglob
LOGS=(logs/*${FILTER}*)
[ ${#LOGS[@]} -gt 0 ] || { echo "no logs matching '${FILTER}' in $PWD/logs"; exit 0; }

echo "== ${#LOGS[@]} log files =="
ok=0; failed=(); warned=()
for f in "${LOGS[@]}"; do
  if grep -q "FATAL" "$f"; then failed+=("$f")
  elif grep -q "===== done" "$f"; then ok=$((ok+1))
  else failed+=("$f")      # no FATAL and no completion: killed, or still running
  fi
  grep -q "WARN" "$f" && warned+=("$f")
done
echo "  completed : $ok"
echo "  failed    : ${#failed[@]}"
echo "  with warnings: ${#warned[@]}"

if [ ${#failed[@]} -gt 0 ]; then
  echo
  echo "== failures (first 15) =="
  for f in "${failed[@]:0:15}"; do
    echo "--- $f"
    if grep -q FATAL "$f"; then
      grep -m3 "FATAL" "$f" | sed 's/^/    /'
    else
      echo "    no FATAL and no completion marker - killed by the scheduler?"
      echo "    (walltime or memory; check the PBS epilogue at the end of the file)"
      tail -3 "$f" | sed 's/^/    | /'
    fi
  done
  echo
  echo "== the distinct failure modes =="
  grep -h "FATAL" "${failed[@]}" 2>/dev/null \
    | sed 's/^\[[0-9:]*\] [^ ]* *FATAL //' \
    | sed 's/[0-9]\{3,\}/N/g' \
    | sort | uniq -c | sort -rn | head -10 | sed 's/^/  /'
fi

if [ ${#warned[@]} -gt 0 ]; then
  echo
  echo "== warnings =="
  grep -h "WARN" "${warned[@]}" 2>/dev/null \
    | sed 's/^\[[0-9:]*\] [^ ]* *WARN *//' | sort | uniq -c | sort -rn | head -5 | sed 's/^/  /'
fi
