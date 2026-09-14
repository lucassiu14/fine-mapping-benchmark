#!/bin/bash
# =============================================================================
# scripts/analysis/submit_iter006_analysis.sh
#
# >>> ITERATION 006 ONLY - TEMPORARY. See docs/autoresearch/iteration-006-REVERT.md
#
# submit_iter005_analysis.sh pointed at Iteration 006. It drives the same R
# scripts - iter005_collect.R, then iter005_report.R and its three senses - which
# read the design from the grid rather than hardcoding it: strata come from
# the relationships and annotation types present. The differences are these:
#   - defaults: $EPHEMERAL/fmbench_iter006 and results/iter006
#   - the grid is regenerated from generate_params_grid_iter006.R
#   - the completeness count is bounded to job_*/scenario_*/results.rds, and
#     SKIP_COUNT=1 skips it (the Iteration 005 script walks the whole tree)
#   - job names fmb6collect / fmb6analyse
#
# Two chained PBS jobs. Stage B is held on afterok, so a partial collect is never
# analysed as if it were complete.
#
#   SKIP_COUNT=1 bash scripts/analysis/submit_iter006_analysis.sh
# =============================================================================
set -euo pipefail

PROJECT_ROOT="$(pwd)"
SCRATCH="${FMB_SCRATCH:-${EPHEMERAL:-/rds/general/user/$USER/ephemeral}/fmbench_iter006}"
BENCH_ROOT="${BENCH_ROOT:-${SCRATCH}/results/benchmark}"

# Analysis output goes to PROJECT SPACE, not ephemeral - ephemeral auto-purges,
# and losing the L1 table is losing the analysis.
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/results/iter006}"
L1_DIR="${L1_DIR:-${OUT_DIR}/L1}"
LOG_DIR="${LOG_DIR:-${OUT_DIR}/logs}"

R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
RSCRIPT="${RSCRIPT:-Rscript}"
# v1_small72a, not v1_small24. On 2026-08-24 v1_small24 rejected even the
# 12-element Stage A array with "qsub: Array job exceeds server or queue size
# limit", while v1_small72a had accepted a 240-element array the same day. The
# precise limit was not diagnosed - it may be array size or a resource cap - so
# this simply defaults to the queue that is known to work. Both stages ask for
# 6h walltime, comfortably inside a 72h queue; override if a shorter queue
# starts sooner for you.
QUEUE_A="${QUEUE_A:-v1_small72a}"
QUEUE_B="${QUEUE_B:-v1_small72a}"
SELECT_A="${SELECT_A:-1:ncpus=1:mem=24gb}"
WALL_A="${WALL_A:-06:00:00}"
SELECT_B="${SELECT_B:-1:ncpus=1:mem=16gb}"
WALL_B="${WALL_B:-06:00:00}"

if [[ ! -d "$BENCH_ROOT" ]]; then
  echo "ERROR: benchmark output not found at $BENCH_ROOT" >&2; exit 1
fi
mkdir -p "$L1_DIR" "$LOG_DIR" "$OUT_DIR"

# Regenerate the grid from the Iteration 006 generator into OUT_DIR. Do NOT
# default to scripts/hpc/params_grid.csv: the submitter overwrites that file on
# every run, so it holds whichever design was submitted last. The report joins
# L1 to it by job_dir and stops on any unmatched row, but a wrong grid that
# happened to match would be worse than a loud failure.
GRID_CSV="${GRID_CSV:-${OUT_DIR}/params_grid_iter006.csv}"
GEN="${PROJECT_ROOT}/scripts/hpc/generate_params_grid_iter006.R"
if [[ ! -f "$GRID_CSV" ]]; then
  [[ -f "$GEN" ]] || { echo "ERROR: generator not found: $GEN" >&2; exit 1; }
  echo "Regenerating the Iteration 006 grid -> $GRID_CSV"
  "$RSCRIPT" "$GEN" "$GRID_CSV" >/dev/null
fi
"$RSCRIPT" "${PROJECT_ROOT}/scripts/hpc/check_grid_columns.R" "$GRID_CSV" \
    "${PROJECT_ROOT}/scripts/hpc/generate_params_grid.R" \
  || { echo "ABORTING: grid is missing columns the pipeline reads." >&2; exit 1; }

N_ROWS=$(find "$BENCH_ROOT" -maxdepth 1 -type d -name 'job_*' | wc -l | tr -d ' ')
(( N_ROWS >= 1 )) || { echo "ERROR: no job_* directories under $BENCH_ROOT" >&2; exit 1; }

# Scenarios per row = |S| * |phi| * n_iter, read from the grid itself.
# fixed = TRUE, NOT the regex "\\|" used elsewhere in the .R files: Rscript -e
# collapses the double backslash before R parses it, and R then rejects "\|" as
# an unrecognised escape. Splitting on a literal pipe sidesteps the quoting
# entirely.
SCEN_PER_ROW=$("$RSCRIPT" --vanilla -e '
  g <- read.csv(commandArgs(TRUE)[1], stringsAsFactors = FALSE)
  cat(length(strsplit(g$S_values[1],   "|", fixed = TRUE)[[1]]) *
      length(strsplit(g$phi_values[1], "|", fixed = TRUE)[[1]]) *
      as.integer(g$n_iter[1]))
' "$GRID_CSV" 2>/dev/null | tail -1)
if ! [[ "$SCEN_PER_ROW" =~ ^[0-9]+$ ]] || (( SCEN_PER_ROW < 1 )); then
  echo "ERROR: could not read scenarios-per-row from $GRID_CSV" >&2; exit 1
fi
EXPECTED=$(( N_ROWS * SCEN_PER_ROW ))
if [[ "${SKIP_COUNT:-0}" == "1" ]]; then
  # The count reads one directory per scenario - not for a login node. Skip it
  # once scripts/hpc/check_iter006.sh shows every fitting task finished; Stage A
  # logs report skipped scenarios per row, and the checker sums them.
  N_SCEN=$EXPECTED
else
  N_SCEN=$(find "$BENCH_ROOT" -mindepth 3 -maxdepth 3 -path '*/job_*/scenario_*/results.rds' | wc -l | tr -d ' ')
fi

echo "Benchmark root : $BENCH_ROOT"
echo "Grid           : $GRID_CSV"
echo "Rows           : $N_ROWS   (scenarios/row: $SCEN_PER_ROW)"
if [[ "${SKIP_COUNT:-0}" == "1" ]]; then echo "Scenarios      : count skipped (SKIP_COUNT=1)"; else echo "Scenarios      : $N_SCEN of $EXPECTED expected"; fi

# A partial grid is not a smaller version of the answer: the design is balanced
# by construction and the noise correction assumes every cell is full.
if (( N_SCEN < EXPECTED )); then
  echo
  echo "WARNING: the run is INCOMPLETE ($(( 100 * N_SCEN / EXPECTED ))%)."
  if [[ "${ALLOW_PARTIAL:-0}" != "1" ]]; then
    echo "         Re-run with ALLOW_PARTIAL=1 to proceed anyway." >&2; exit 1
  fi
  echo "         ALLOW_PARTIAL=1 set - proceeding."
fi

# --- Stage A: collect, one task per grid row ---------------------------------
A_SCRIPT="$(mktemp -t fmb6a_XXXXXX.sh)"
cat > "$A_SCRIPT" <<EOF
#!/bin/bash
#PBS -N fmb6collect
#PBS -q ${QUEUE_A}
#PBS -l select=${SELECT_A}
#PBS -l walltime=${WALL_A}
#PBS -J 1-${N_ROWS}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
echo "[collect \${PBS_ARRAY_INDEX} on \$(hostname)] start \$(date)"
Rscript scripts/analysis/iter005_collect.R \\
    "\${PBS_ARRAY_INDEX}" "${BENCH_ROOT}" "${L1_DIR}"
echo "[collect \${PBS_ARRAY_INDEX}] done \$(date)"
EOF
A_ID="$(qsub "$A_SCRIPT")"
echo
echo "Stage A (collect, ${N_ROWS} tasks): ${A_ID}"

# --- Stage B: bind, validate, drive the three senses -------------------------
B_SCRIPT="$(mktemp -t fmb6b_XXXXXX.sh)"
cat > "$B_SCRIPT" <<EOF
#!/bin/bash
#PBS -N fmb6analyse
#PBS -q ${QUEUE_B}
#PBS -l select=${SELECT_B}
#PBS -l walltime=${WALL_B}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/
set -euo pipefail
cd "${PROJECT_ROOT}"
module load ${R_MODULE}
echo "[analyse on \$(hostname)] start \$(date)"
Rscript scripts/analysis/iter005_report.R "${L1_DIR}" "${GRID_CSV}" "${OUT_DIR}"
echo "[analyse] done \$(date)"
EOF
B_ID="$(qsub -W depend=afterok:"${A_ID}" "$B_SCRIPT")"
echo "Stage B (analyse, held on afterok): ${B_ID}"

cat <<EOF

Output    : ${OUT_DIR}
  combined_fit_metrics.rds        L1 - every statistic is re-derivable from it
  combined_replicate_metrics.rds  L2 - paired differences for the Sense D gate
  combined_scenario_metrics.rds   L3 - cell means, the ANOVA input
  validity_checks.txt             read this FIRST
  iter005_sense_v/sense_v.txt     which variables move the metric
  iter005_sense_d/sense_d_ap.txt  which variables change the winner
  iter005_sense_g/sense_g_ap.txt  which variables hurt THIS method specifically

Strata    : relationship x annotation_type = 8 cells (continuous only). NOT
            model x annotation_type - every row here is model="sparse", so those
            keys would pool the seven non-null relationships into one number.

Track     : qstat -tan -u \$USER | grep -E 'fmb6collect|fmb6analyse'
EOF
