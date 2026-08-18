#!/bin/bash
# scripts/hpc/submit_benchmark_pbs.sh
# PBS Pro driver for the fine-mapping benchmark array (Imperial HPC).
# SLURM users: see scripts/hpc/submit_benchmark.sh.

set -euo pipefail

PROJECT_ROOT="$(pwd)"
RSCRIPT="${RSCRIPT:-Rscript}"
PARAMS_CSV="${PROJECT_ROOT}/scripts/hpc/params_grid.csv"

# Output + logs go on scratch, NOT home. Each task writes ~100 MB+ of RDS
# files; a full array overflows a 1 TB home quota and saveRDS() dies with
# "error writing to connection". Default to personal ephemeral (multi-TB).
# Override with FMB_SCRATCH, or FMB_OUTPUT_ROOT / FMB_LOG_DIR individually.
FMB_SCRATCH="${FMB_SCRATCH:-${EPHEMERAL:-/rds/general/user/$USER/ephemeral}/fmbench}"
OUTPUT_ROOT="${FMB_OUTPUT_ROOT:-${FMB_SCRATCH}/results/benchmark}"
LOG_DIR="${FMB_LOG_DIR:-${FMB_SCRATCH}/logs/benchmark}"
export FMB_OUTPUT_ROOT="${OUTPUT_ROOT}"

# Config (override with env vars: PBS_QUEUE=... bash submit_benchmark_pbs.sh)
PBS_QUEUE="${PBS_QUEUE:-v1_small72a}"
PBS_WALLTIME="${PBS_WALLTIME:-72:00:00}"       # queue max; BEATRICE-family
                                                # is compute-heavy so start
                                                # generous, cut later if fine.
# 32gb: a task holds the row's full region panel (all regions must stay
# together for cross-region pooling in scenario_setup) + ONE scenario's
# results at a time (the chunk loop processes scenarios one-by-one, not all
# at once) + a torch subprocess for BEATRICE/FB. Peak memory is set by the
# sim + a single scenario, NOT the chunk size, so chunking does not raise it;
# Iteration 002's 10-region sim is lighter than Iteration 001's 20-region one.
# The earlier 8gb OOM was the OLD full-row task accumulating all scenarios'
# results in memory. 32gb keeps a safe margin and still packs ~4 tasks per
# 128gb node.
PBS_SELECT="${PBS_SELECT:-1:ncpus=1:mem=32gb}"
ARRAY_RANGE="${ARRAY_RANGE:-}"                  # e.g. "1-2" canary, "" full

R_MODULE="${R_MODULE:-R/4.5.2-gfbf-2025b}"
GSL_MODULE="${GSL_MODULE:-GSL/2.8-GCC-14.3.0}"
PYTHON_MODULE="${PYTHON_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
PY_VENV_ACTIVATE="${PY_VENV_ACTIVATE:-$HOME/tools/fmpy-venv/bin/activate}"

# ALWAYS regenerate the grid.
#
# This used to be "generate if missing". params_grid.csv is untracked but
# persists on disk between iterations, so after the Iteration 004 grid revision
# the stale Iteration 003 file (135 rows, n=1000, regions 100-1000, three LD
# levels) was silently reused and a whole array was submitted against the wrong
# simulation design. The generator is deterministic and takes under a second.
echo "Regenerating params_grid.csv from generate_params_grid.R ..."
"$RSCRIPT" "${PROJECT_ROOT}/scripts/hpc/generate_params_grid.R" "$PARAMS_CSV" >/dev/null
# The array is subdivided by scenario, with all regions kept together in each
# task and (by default) several scenarios batched per task. SCENARIOS_PER_ROW =
# |S| * |phi| * n_iter; the chunk size collapses that to TASKS_PER_ROW tasks.
# Columns are located by header name so column order can change.
N_ROWS=$(( $(wc -l < "$PARAMS_CSV") - 1 ))
if (( N_ROWS < 1 )); then
  echo "ERROR: params_grid.csv has no rows" >&2
  exit 1
fi
_col() { head -1 "$PARAMS_CSV" | tr ',' '\n' | grep -n "\"$1\"" | cut -d: -f1; }
S_COL=$(_col S_values); PHI_COL=$(_col phi_values); NITER_COL=$(_col n_iter)
if [[ -z "$S_COL" || -z "$PHI_COL" || -z "$NITER_COL" ]]; then
  echo "ERROR: params_grid.csv missing S_values/phi_values/n_iter columns" >&2
  exit 1
fi
SCENARIOS_PER_ROW=$(awk -F, -v s="$S_COL" -v ph="$PHI_COL" -v ni="$NITER_COL" '
  NR==2 {
    gsub(/"/,"",$s); gsub(/"/,"",$ph); gsub(/"/,"",$ni);
    ns = split($s, a, "[|]"); np = split($ph, b, "[|]");
    print ns * np * ($ni + 0);
  }' "$PARAMS_CSV")
# Scenario chunking: each task processes SCENARIOS_PER_TASK scenarios from one
# row (see run_benchmark_job.R). This collapses N_ROWS*SCENARIOS_PER_ROW down to
# N_ROWS*ceil(SCENARIOS_PER_ROW/chunk) tasks so a large grid fits under the PBS
# array cap while each (fatter) task still finishes inside walltime. The worker
# reads the SAME value from FMB_SCENARIOS_PER_TASK (exported into the job below),
# so the two always agree.
SCENARIOS_PER_TASK="${FMB_SCENARIOS_PER_TASK:-25}"
if (( SCENARIOS_PER_TASK < 1 )); then SCENARIOS_PER_TASK=1; fi
if (( SCENARIOS_PER_TASK > SCENARIOS_PER_ROW )); then SCENARIOS_PER_TASK=$SCENARIOS_PER_ROW; fi
export FMB_SCENARIOS_PER_TASK="${SCENARIOS_PER_TASK}"
TASKS_PER_ROW=$(( (SCENARIOS_PER_ROW + SCENARIOS_PER_TASK - 1) / SCENARIOS_PER_TASK ))
N_JOBS=$(( N_ROWS * TASKS_PER_ROW ))
# Echo the simulation design that will actually run. A grid mismatch is far
# cheaper to catch here than after 72h of wall-clock. Done in R rather than
# awk because the CSV is quoted and column lookup by name is trivial there.
"$RSCRIPT" -e '
  g <- read.csv(commandArgs(TRUE)[1], stringsAsFactors = FALSE)
  u <- function(x) paste(sort(unique(x)), collapse = ", ")
  cat(sprintf("DESIGN: n=%s   regions={%s}\n", u(g$n),
              gsub("[|]", ",", g$p_values[1])))
  cat(sprintf("        S={%s}   phi={%s}\n",
              gsub("[|]", ",", g$S_values[1]), gsub("[|]", ",", g$phi_values[1])))
  cat(sprintf("        LD: %s\n",
              u(ifelse(is.na(g$n_ref), "in-sample", paste0("ref", g$n_ref)))))
  cat(sprintf("        p_causal: %s   annotations: %s\n",
              u(g$p_causal[!is.na(g$p_causal)]), u(g$annotation_type)))
' "$PARAMS_CSV"
echo "Grid: ${N_ROWS} rows x ${SCENARIOS_PER_ROW} scenarios; chunk=${SCENARIOS_PER_TASK} scenario(s)/task"
echo "      -> ${TASKS_PER_ROW} tasks/row x ${N_ROWS} rows = ${N_JOBS} array tasks"

# PBS caps the number of array elements per submission. If the chunked count
# still exceeds it, tell the user to raise the chunk rather than silently
# submitting an array PBS will reject.
MAX_ARRAY="${FMB_MAX_ARRAY:-10000}"
if (( N_JOBS > MAX_ARRAY )); then
  NEED=$(( (SCENARIOS_PER_ROW * N_ROWS + MAX_ARRAY - 1) / MAX_ARRAY ))
  echo "ERROR: ${N_JOBS} tasks exceeds the array cap FMB_MAX_ARRAY=${MAX_ARRAY}." >&2
  echo "       Raise the chunk, e.g. FMB_SCENARIOS_PER_TASK=${NEED} bash $0" >&2
  exit 1
fi
if [[ -z "$ARRAY_RANGE" ]]; then ARRAY_RANGE="1-${N_JOBS}"; fi

mkdir -p "$LOG_DIR" "$OUTPUT_ROOT"
echo "Output root: $OUTPUT_ROOT"
echo "Log dir:     $LOG_DIR"

JOB_SCRIPT="$(mktemp -t fmbench_pbs_XXXXXX.sh)"
cat > "$JOB_SCRIPT" <<PBS_EOF
#!/bin/bash
#PBS -N fmbench
#PBS -q ${PBS_QUEUE}
#PBS -l select=${PBS_SELECT}
#PBS -l walltime=${PBS_WALLTIME}
#PBS -J ${ARRAY_RANGE}
#PBS -o ${LOG_DIR}/
#PBS -e ${LOG_DIR}/

set -euo pipefail
cd "${PROJECT_ROOT}"
# CRITICAL: export the scratch output root INTO the compute-node job.
# The login-side export does not propagate to the node, so without this
# the worker falls back to results/ under home and overflows the quota.
# \${OUTPUT_ROOT} is expanded here at submit time to the literal path.
export FMB_OUTPUT_ROOT="${OUTPUT_ROOT}"
# Supplemental re-run selector, forwarded the same way (a login-side export
# does NOT reach the node). Empty string = normal full-method run.
export FMB_METHODS="${FMB_METHODS:-}"
# Scenario chunk size. MUST match the value the submit script used to size the
# array, or the worker's (row, scenario-block) decode won't line up. Expanded
# here at submit time to the literal integer.
export FMB_SCENARIOS_PER_TASK="${SCENARIOS_PER_TASK}"
module load ${R_MODULE}
module load ${GSL_MODULE}
module load ${PYTHON_MODULE}
source ${PY_VENV_ACTIVATE}
# PREFLIGHT: sourcing activate proves the venv PATH exists, NOT that its packages
# import. Iteration 004 lost the entire sparse_inf arm to this gap - activation
# succeeded, then every Python-backed method (beatrice, functional_beatrice,
# fb_pooled, fb_xregion, sparsepro, finemap_inf, funmap) died on ImportError
# inside its subprocess. The R wrapper caught each one and recorded a failed fit,
# so the job ran for hours writing 40-50% NA instead of stopping, and the damage
# was only visible after the analysis. Fail the task at second zero instead.
# Preflight via the SAME derived scan check_toolchain.sh uses. An earlier
# version listed five modules by hand - numpy/scipy/torch/absl/pandas - all of
# which imported fine while BEATRICE's imageio, seaborn and tqdm and
# FINEMAP-inf's bgzip were absent. The preflight passed and every fit in the
# array died on ModuleNotFoundError. Derive the list from the tool sources so
# the check cannot drift from what the methods actually import.
_TOOLS="\$(readlink -f "\$HOME/tools" 2>/dev/null || echo "\$HOME/tools")"
if [[ -f "${PROJECT_ROOT}/scripts/hpc/py_deps_check.py" ]]; then
  python "${PROJECT_ROOT}/scripts/hpc/py_deps_check.py" \\
      "${PROJECT_ROOT}/BEATRICE_annot_sparse" \\
      "\$_TOOLS/SparsePro" "\$_TOOLS/fine-mapping-inf" "\$_TOOLS/Funmap" \\
    || { echo "PREFLIGHT FAILED - see the missing/broken packages above" >&2; exit 1; }
  echo "preflight ok - dependency scan clean" >&2
else
  echo "PREFLIGHT FAILED - py_deps_check.py missing at ${PROJECT_ROOT}/scripts/hpc/" >&2
  exit 1
fi
echo "[node:\$(hostname)] task \${PBS_ARRAY_INDEX} of ${ARRAY_RANGE} starting at \$(date)"
echo "[node:\$(hostname)] FMB_OUTPUT_ROOT=\${FMB_OUTPUT_ROOT}"
${RSCRIPT} scripts/hpc/run_benchmark_job.R "\${PBS_ARRAY_INDEX}"
echo "[node:\$(hostname)] task \${PBS_ARRAY_INDEX} finished at \$(date)"
PBS_EOF

echo "Job script: $JOB_SCRIPT"
echo "Submitting array ${ARRAY_RANGE} to queue ${PBS_QUEUE} ..."
qsub "$JOB_SCRIPT"

echo
echo "Track with:"
echo "  qstat -tan \$USER          # per-array-element status"
echo "  ls -lh ${LOG_DIR}          # log files as they land"
