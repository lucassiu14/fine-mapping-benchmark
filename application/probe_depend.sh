#!/bin/bash
# Discover which array-job dependency syntax THIS scheduler accepts.
#
# Reasoning about PBS dependency syntax from documentation proved unreliable:
# qsub returns an array id as "1234[].server", and both
# "afterokarray:1234[].server" and "afterokarray:1234[].server[]" are
# rejected with the uninformative "qsub: illegal -W value". So measure it.
#
# Submits a HELD 2-task array, then tries each candidate dependency on a HELD
# one-line job, reporting accepted/rejected per form. Everything is deleted
# again; nothing runs. Held jobs consume a queue slot for seconds, not nodes.
#
#   bash probe_depend.sh
set -uo pipefail
cd "${APP_ROOT:?source env.sh first}"

TMP=$(mktemp -d "${TMPDIR:-/tmp}/depprobe.XXXXXX")
cleanup () { [ -n "${AID:-}" ] && qdel "$AID" 2>/dev/null; rm -rf "$TMP"; }
trap cleanup EXIT

cat > "$TMP/arr.pbs" <<'EOF'
#!/bin/bash
#PBS -N depprobe_arr
#PBS -l select=1:ncpus=1:mem=100mb
#PBS -l walltime=00:01:00
#PBS -J 1-2
#PBS -j oe
true
EOF
cat > "$TMP/dep.pbs" <<'EOF'
#!/bin/bash
#PBS -N depprobe_dep
#PBS -l select=1:ncpus=1:mem=100mb
#PBS -l walltime=00:01:00
#PBS -j oe
true
EOF

echo "qsub version: $(qsub --version 2>&1 | head -1)"
echo
AID=$(cd "$TMP" && qsub -h arr.pbs) || { echo "FATAL: could not submit the probe array"; exit 1; }
echo "probe array (held): $AID"

# Decompose the id. For "1234[].pbs-7": SEQ=1234  SRV=pbs-7
SEQ=${AID%%[*}
SRV=${AID##*].}
[ "$SRV" = "$AID" ] && SRV=""      # no "]." in the id
echo "  sequence number : $SEQ"
echo "  server suffix   : ${SRV:-<none>}"
echo

try () {
  local label="$1" dep="$2" out rc
  out=$(cd "$TMP" && qsub -h -W "depend=$dep" dep.pbs 2>&1); rc=$?
  if [ "$rc" -eq 0 ]; then
    printf '  ACCEPTED   %-28s  %s\n' "$label" "$dep"
    qdel "$out" 2>/dev/null
  else
    printf '  rejected   %-28s  %s\n' "$label" "$dep"
    printf '             %-28s  -> %s\n' "" "$(echo "$out" | head -1)"
  fi
}

echo "array dependency forms:"
try "as qsub returned it"      "afterokarray:${AID}"
try "brackets moved to end"    "afterokarray:${SEQ}${SRV:+.$SRV}[]"
try "no server suffix"         "afterokarray:${SEQ}[]"
try "no brackets"              "afterokarray:${SEQ}${SRV:+.$SRV}"
try "bare sequence number"     "afterokarray:${SEQ}"
try "doubled brackets"         "afterokarray:${AID}[]"
echo
echo "non-array type, for contrast:"
try "afterok, as returned"     "afterok:${AID}"
try "afterok, brackets at end" "afterok:${SEQ}${SRV:+.$SRV}[]"
echo
echo "Report the ACCEPTED line(s); submit_all.sh will be built to match."
