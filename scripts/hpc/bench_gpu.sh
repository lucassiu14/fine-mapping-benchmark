#!/bin/bash
# =============================================================================
# bench_gpu.sh - one-scenario runtime benchmark for fb_xregion's joint trainer.
#
# Answers three questions on REAL cluster hardware, on identical data:
#   1. how much the arithmetic changes (contracted Woodbury + eigvalsh + batched
#      prior head) actually buy on a CX3 core, vs the pre-change baseline;
#   2. whether the GPU path runs at all;
#   3. whether the GPU is faster or slower than the core for this workload.
#
# It does NOT touch the benchmark grid, results/ or anything under $EPHEMERAL's
# iteration directories. Everything happens in $TMPDIR on the compute node.
#
# Usage, from the repo root on a login node:
#     qsub scripts/hpc/bench_gpu.sh
#     # CPU-only node (skips the GPU leg):
#     qsub -q v1_small72a -l select=1:ncpus=1:mem=32gb scripts/hpc/bench_gpu.sh
#
# The results land in the job's .o file as a single table.
# =============================================================================
#PBS -N fbgpu_bench
#PBS -q v1_gpu72
#PBS -l select=1:ncpus=4:mem=64gb:ngpus=1
#PBS -l walltime=03:00:00
#PBS -j oe

set -uo pipefail

PROJECT_ROOT="${PBS_O_WORKDIR:-$PWD}"
BASE_COMMIT="${BASE_COMMIT:-d8e8e73}"        # pre-change reference
MAX_ITER="${MAX_ITER:-1500}"                 # production value for fb_xregion
WORK="${TMPDIR:-/tmp}/fbbench.$$"
mkdir -p "$WORK"

# Overridable so the script can be exercised off-cluster before it costs queue
# time; the defaults match scripts/hpc/submit_benchmark_pbs.sh.
PYTHON_MODULE="${PYTHON_MODULE:-Python/3.12.3-GCCcore-13.3.0}"
PY_VENV_ACTIVATE="${PY_VENV_ACTIVATE:-$HOME/tools/fmpy-venv/bin/activate}"
if [[ "${SKIP_MODULES:-0}" != "1" ]]; then
  module load "$PYTHON_MODULE" 2>/dev/null
  source "$PY_VENV_ACTIVATE" || { echo "venv activate failed: $PY_VENV_ACTIVATE"; exit 1; }
fi

cd "$PROJECT_ROOT" || exit 1
echo "commit under test : $(git rev-parse --short HEAD) ($(git rev-parse --abbrev-ref HEAD))"
echo "baseline commit   : $BASE_COMMIT"
echo "node              : $(hostname)"
python - <<'PY'
import torch
print(f"torch             : {torch.__version__}  cuda_build={torch.version.cuda}  "
      f"available={torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"gpu               : {torch.cuda.get_device_name(0)}")
PY
echo

# --- baseline package at $BASE_COMMIT, extracted read-only -------------------
mkdir -p "$WORK/base"
git archive "$BASE_COMMIT" BEATRICE_annot_sparse | tar -x -C "$WORK/base" || {
  echo "could not extract $BASE_COMMIT"; exit 1; }

# --- one production-shaped scenario: 10 regions, 2 at each of the 5 sizes ----
python - "$WORK" <<'PY'
import numpy as np, os, sys
W = sys.argv[1]; d = os.path.join(W, "data"); os.makedirs(d, exist_ok=True)
rng = np.random.default_rng(21)
rows, sizes = [], [500, 500, 750, 750, 1000, 1000, 1500, 1500, 2000, 2000]
for r, m in enumerate(sizes, start=1):
    names = np.array([f"rs{r}_{j}" for j in range(m)])
    X = rng.normal(size=(2500, m))
    for j in range(1, m):
        X[:, j] = 0.6 * X[:, j - 1] + 0.8 * X[:, j]
    X = (X - X.mean(0)) / X.std(0)
    LD = np.corrcoef(X, rowvar=False)
    z = rng.normal(size=m); z[[10, m // 2]] += 6.0
    A = rng.integers(0, 2, size=(m, 10)).astype(float); A[[10, m // 2], :5] = 1.0
    np.savetxt(f"{d}/r{r}.z", np.column_stack([names, z]), fmt="%s", delimiter=" ")
    np.savetxt(f"{d}/r{r}.ld", LD, fmt="%.6f", delimiter=" ")
    np.savetxt(f"{d}/r{r}.annot", np.column_stack([names, A]), fmt="%s", delimiter=" ")
    rows.append(f"{d}/r{r}.z\t{d}/r{r}.ld\t{d}/r{r}.annot\t{{OUT}}/out{r}\t5000")
open(f"{d}/manifest.template", "w").write("\n".join(rows) + "\n")
print(f"scenario: 10 regions {sizes}, 10 annotations")
PY
echo

# --- run one configuration and report wall time ------------------------------
run_cfg () {   # $1=tag  $2=label  $3=package root  $4=FMBENCH_DEVICE  $5=threads
  local TAG="$1" LABEL="$2" PKG="$3" DEV="$4" THREADS="$5"
  local OUT="$WORK/out_$TAG"; rm -rf "$OUT"; mkdir -p "$OUT"
  sed "s|{OUT}|$OUT|g" "$WORK/data/manifest.template" > "$OUT/manifest.tsv"
  local t0 t1 rc rate
  # Cap each leg. If the accelerator turns out to be launch-bound and slower
  # than the core - which results.tex sec 5.6.1 predicts for the un-batched
  # region loop - we want that reported as a timeout, not discovered when PBS
  # kills the job at walltime with nothing printed.
  local TO=""; command -v timeout >/dev/null && TO="timeout ${LEG_TIMEOUT_S:-5400}"
  t0=$(date +%s)
  ( cd "$PKG" && OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" \
      FMBENCH_DEVICE="$DEV" $TO python beatrice_joint.py \
        --manifest "$OUT/manifest.tsv" --prior_head lassonet --max_iter "$MAX_ITER" \
        --n_caus 5 --sigma_sq 0.05 --sparse_concrete 50 ) > "$OUT/log.txt" 2>&1
  rc=$?; t1=$(date +%s)
  rate=$(grep -o '[0-9.]*it/s' "$OUT/log.txt" | tail -1)
  [ "$rc" -eq 124 ] && rate="TIMED OUT >${LEG_TIMEOUT_S:-5400}s"
  printf '  %-34s %5ss   %-12s exit=%s\n' "$LABEL" "$((t1-t0))" "${rate:--}" "$rc"
  [ "$rc" -ne 0 ] && { echo "    --- last 12 lines ---"; tail -12 "$OUT/log.txt" | sed 's/^/    /'; }
  return 0
}

echo "=== one scenario, ${MAX_ITER} iterations ==="
printf '  %-34s %6s   %-12s\n' "configuration" "wall" "train rate"
run_cfg base "baseline @ $BASE_COMMIT (cpu, 1 thr)" "$WORK/base/BEATRICE_annot_sparse" cpu 1
run_cfg newcpu "this branch (cpu, 1 thr)"       "$PROJECT_ROOT/BEATRICE_annot_sparse" cpu 1
if python -c "import torch,sys; sys.exit(0 if torch.cuda.is_available() else 1)"; then
  run_cfg newgpu "this branch (cuda)"           "$PROJECT_ROOT/BEATRICE_annot_sparse" cuda 4
else
  echo "  (no GPU on this node - cuda leg skipped)"
fi

echo
echo "=== outputs written by the branch run ==="
ls "$WORK/out_newcpu/out1" 2>/dev/null | sed 's/^/  /'
echo
echo "work dir (node-local, purged with the job): $WORK"
