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
  # DERIVE the required imports from the tools' own source rather than listing
  # them here. A hand-written list is how the second array was lost: it read
  # numpy/scipy/pandas/torch/absl-py, every one of which imported fine, while
  # BEATRICE also needs imageio, seaborn and tqdm and FINEMAP-inf needs bgzip.
  # The check passed, the array ran, and all 2,520 fits died on
  # ModuleNotFoundError. The source is on disk; there is no reason to guess.
  DEPS="$TOOLS_ROOT/py_deps_check.py"
  [[ -f "$DEPS" ]] || DEPS="$PROJECT_ROOT/scripts/hpc/py_deps_check.py"
  if [[ -f "$DEPS" ]]; then
    dep_out="$("$RUNNER" "$DEPS" \
        "$PROJECT_ROOT/BEATRICE_annot_sparse" \
        "$TOOLS_ROOT/SparsePro" "$TOOLS_ROOT/fine-mapping-inf" \
        "$TOOLS_ROOT/Funmap" 2>&1)"; dep_rc=$?
    while IFS= read -r line; do
      [[ -n "$line" ]] && printf "             %s\n" "$line"
    done < <(printf '%s\n' "$dep_out")
    if (( dep_rc == 0 )); then
      okay "every import the tool sources actually use resolves"
    else
      bad "missing packages - see the pip line above"
    fi
  else
    bad "py_deps_check.py not found - cannot verify the real dependency set"
  fi
  # ASSERT the interpreter is the venv's, do not merely report it. Imports
  # succeeding proves nothing on its own: if the venv's interpreter is broken,
  # a bare `python` resolves to whatever else is on PATH, and should that
  # interpreter happen to carry numpy/scipy/torch the whole check passes while
  # every method silently runs against an environment nobody configured. This
  # assertion is the difference between a preflight and a placebo.
  echo
  got_prefix="$("$RUNNER" -c 'import sys; print(sys.prefix)' 2>&1 | tail -1)"
  want_prefix="$TOOLS_ROOT/fmpy-venv"
  printf "  %-10s %s\n" "[info]" "prefix: $got_prefix"
  if [[ "$got_prefix" == "$want_prefix" ]]; then
    okay "interpreter is the benchmark venv"
  else
    bad "WRONG INTERPRETER - expected $want_prefix"
    say "" "the runner fell through to another python; the venv is broken"
  fi
else
  bad "py-venv-runner.sh not executable - cannot test imports"
fi

echo
if (( FAIL )); then
  echo "RESULT: FAIL - do not submit. Rebuild with scripts/hpc/rebuild_toolchain.sh"
  exit 1
fi
echo "RESULT: PASS - toolchain is complete and importable."
