#!/bin/bash
# =============================================================================
# scripts/hpc/check_toolchain.sh
#
# Verify EVERY external path run_benchmark_job.R expects, before spending an
# array on it.
#
# WHY THIS EXISTS. Iteration 004 lost its whole sparse_inf arm because ~/tools
# was a symlink onto $EPHEMERAL (created Jul 11 10:38) and Imperial's rolling
# 30-day purge ate it file by file during the run. Sourcing the venv's activate
# script kept succeeding while the packages inside it stopped importing, so
# seven Python-backed methods recorded 40-50% failed fits for hours instead of
# stopping. By the time the run finished, finemap and PAINTOR's binaries had
# gone too - they had been fine DURING the run, so a naive re-run would have
# failed in a fresh set of methods.
#
# The lesson is that "the tool directory exists" proves nothing. Check the leaf
# artefacts, and check that the interpreter can actually import.
#
#   bash scripts/hpc/check_toolchain.sh
#   TOOLS_ROOT=/rds/general/project/.../tools bash scripts/hpc/check_toolchain.sh
#
# Exit 0 = safe to submit. Exit 1 = do not submit.
# =============================================================================
set -uo pipefail

TOOLS_ROOT="${TOOLS_ROOT:-$(readlink -f "$HOME/tools" 2>/dev/null || echo "$HOME/tools")}"
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
FAIL=0

say()  { printf "  %-10s %s\n" "$1" "$2"; }
okay() { say "[ok]"   "$1"; }
bad()  { say "[MISS]" "$1"; FAIL=1; }

echo "Toolchain check"
echo "  TOOLS_ROOT = $TOOLS_ROOT"

# A tools root on ephemeral is the actual root cause, so say so loudly even if
# every file happens to be present today.
case "$TOOLS_ROOT" in
  *ephemeral*)
    echo
    echo "  !! TOOLS_ROOT resolves onto EPHEMERAL storage."
    echo "  !! Imperial purges ephemeral on a rolling ~30-day window. A toolchain"
    echo "  !! here will decay MID-RUN and the failures surface as per-method"
    echo "  !! ImportErrors, not as a job crash. Move it to project storage."
    FAIL=1
    ;;
esac
echo

echo "Binaries and tool directories"
for p in \
  "py-venv-runner.sh" \
  "fmpy-venv/bin/activate" \
  "fmpy-venv/pyvenv.cfg" \
  "SparsePro" \
  "fine-mapping-inf" \
  "finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64" \
  "PAINTOR_V3.0/PAINTOR"
do
  [[ -e "$TOOLS_ROOT/$p" ]] && okay "$p" || bad "$p"
done

# BEATRICE lives in the PROJECT, not in tools - it is a working copy, so it is
# checked separately and against the project root.
echo
echo "Project-local"
[[ -d "$PROJECT_ROOT/BEATRICE_annot_sparse" ]] \
  && okay "BEATRICE_annot_sparse" || bad "BEATRICE_annot_sparse"

# The decisive test. Directory presence proved nothing in Iteration 004; what
# matters is whether the interpreter the methods are handed can import.
echo
echo "Python imports (via the runner the methods are actually given)"
RUNNER="$TOOLS_ROOT/py-venv-runner.sh"
if [[ -x "$RUNNER" ]]; then
  for m in numpy scipy.special pandas torch absl.flags; do
    if "$RUNNER" -c "import $m" 2>/dev/null; then okay "import $m"; else bad "import $m"; fi
  done
else
  bad "py-venv-runner.sh not executable - cannot test imports"
fi

echo
if (( FAIL )); then
  echo "RESULT: FAIL - do not submit. Rebuild with scripts/hpc/rebuild_toolchain.sh"
  exit 1
fi
echo "RESULT: PASS - toolchain is complete and importable."
