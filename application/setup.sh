#!/bin/bash
# One-shot bootstrap: data, upstream checkouts, and the three Python
# environments. Safe to re-run - every step skips work already done.
#
#   export APP_ROOT=/rds/general/user/<you>/home/chol_app
#   bash setup.sh
#   source "$APP_ROOT/env.sh"      # then use submit_all.sh
set -euo pipefail
: "${APP_ROOT:?export APP_ROOT to the application directory first}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_BOOT="${PY_BOOT:-python3}"

mkdir -p "$APP_ROOT"/{ld,annot,regions,polyfun_priors,logs,results}
cd "$APP_ROOT"

say() { printf '\n=== %s ===\n' "$1"; }
bail() { printf '\nFATAL  %s\n' "$1" >&2; exit 1; }

say "pipeline scripts"
# lib.sh in particular: every stage does `source ./lib.sh` relative to
# APP_ROOT, so without it nothing runs at all.
if [ "$(cd "$HERE" && pwd)" != "$(cd "$APP_ROOT" && pwd)" ]; then
  cp -f "$HERE"/*.py "$HERE"/*.pbs "$HERE"/*.sh "$APP_ROOT"/
else
  echo "  (running from APP_ROOT; scripts already in place)"
fi
cp -f "$HERE"/blocks_scored.csv "$HERE"/blocks.csv "$APP_ROOT"/ 2>/dev/null || true

# ld_windows.txt: the distinct LD windows, one per line, which 01's array
# indexes by line number. It is .gitignore'd (it is derived, not source) and
# only 00_make_blocks.py wrote it -- so it never reached APP_ROOT, and every
# task of 01 exited with "no window at line N" having downloaded nothing.
# Derive it from blocks.csv here instead of committing it, so the two cannot
# drift. The column is looked up by NAME, not position, and duplicates are
# dropped in first-appearance order, matching 00_make_blocks.py's
# drop_duplicates().
if [ -s "$APP_ROOT/blocks.csv" ]; then
  awk -F, 'NR==1 { for (i=1; i<=NF; i++) if ($i == "ld_region") c = i;
                   if (!c) { print "blocks.csv has no ld_region column" > "/dev/stderr"; exit 1 }
                   next }
           !seen[$c]++ { print $c }' "$APP_ROOT/blocks.csv" > "$APP_ROOT/ld_windows.txt.tmp" \
    && mv "$APP_ROOT/ld_windows.txt.tmp" "$APP_ROOT/ld_windows.txt"
  echo "  ld_windows.txt: $(wc -l < "$APP_ROOT/ld_windows.txt") distinct LD windows"
fi
cp -f "$HERE"/../inst/extdata/gwfm_regions_ldetect_EUR.csv "$APP_ROOT"/ 2>/dev/null \
  || echo "  NOTE: copy gwfm_regions_ldetect_EUR.csv here if you want to re-derive blocks_scored.csv"

say "GWAS summary statistics (526 MB)"
if [ -s chol_sumstats.parquet ]; then echo "  already present"; else
  curl -fL --retry 5 --retry-delay 10 -o chol_sumstats.parquet.part \
    "https://storage.googleapis.com/seqnn-share/sniff/_20260810_tmp_files/fixed_n_337k/sample_munged_gwas/biochemistry_Cholesterol.sumstats_munged.parquet"
  mv chol_sumstats.parquet.part chol_sumstats.parquet
fi

say "upstream comparators"
[ -d "$APP_ROOT/Funmap"  ] || git clone -q --depth 1 https://github.com/LeeHITsz/Funmap.git  "$APP_ROOT/Funmap"
[ -d "$APP_ROOT/polyfun" ] || git clone -q --depth 1 https://github.com/omerwe/polyfun.git "$APP_ROOT/polyfun"
echo "  Funmap and polyfun checked out (upstream, unmodified)"

# Every pip call below runs with PYTHONPATH and PYTHONHOME CLEARED, and this
# is the whole point rather than tidiness.
#
# EasyBuild's Python modules put their own site-packages on PYTHONPATH, and a
# venv does NOT shadow PYTHONPATH. So with a module loaded, pip saw
# `packaging` already present at
#   /sw-eb/.../Python-bundle-PyPI/2023.10-GCCcore-13.2.0/.../packaging/
# declared it satisfied, and installed nothing into appenv. pip check agreed,
# and `import packaging` worked on the login node. The compute node loads only
# APP_MODULES -- not the awscli bundle that path came from -- so the import
# failed there, twice, after a queue wait each time. absl-py had the same
# shape.
#
# Clearing PYTHONPATH makes pip see only the venv, so every listed package
# genuinely lands in it and the venv behaves the same wherever it runs.
mkenv () {  # name, then pip args
  local n="$1"; shift
  local py="$APP_ROOT/$n/bin/python" log="$APP_ROOT/logs/pip_$n.log"
  if [ -x "$py" ]; then
    # Already built - but the package list may have GROWN since, so install
    # the current list rather than skipping, otherwise a dependency added here
    # never reaches a cluster that ran setup.sh before it was added.
    #
    # The failure is NOT swallowed. This previously read
    #   pip install -q "$@" && echo ok || echo "WARNING: could not reconcile"
    # which, with -q hiding the output, let a broken reconcile pass as a
    # warning and setup continue. appenv then carried matplotlib without
    # `packaging`, and 03 died on the import an hour later.
    if env -u PYTHONPATH -u PYTHONHOME "$py" -m pip install -q "$@" >"$log" 2>&1; then
      echo "  $n already built (package list reconciled)"
    else
      sed 's/^/    /' "$log" >&2
      bail "$n: could not reconcile packages - see $log"
    fi
  else
    "$PY_BOOT" -m venv "$APP_ROOT/$n"
    env -u PYTHONPATH -u PYTHONHOME "$py" -m pip install -q --upgrade pip >"$log" 2>&1 || true
    if ! env -u PYTHONPATH -u PYTHONHOME "$py" -m pip install -q "$@" >>"$log" 2>&1; then
      sed 's/^/    /' "$log" >&2
      bail "$n: pip install failed - see $log"
    fi
    echo "  $n built"
  fi
  # VERIFY, do not assume. `pip check` reports any installed package whose own
  # declared dependencies are missing - exactly how appenv ended up holding
  # matplotlib but not `packaging`, passing an `import matplotlib` test at one
  # moment and failing it later. A package LIST is not a working environment,
  # and scanning an application's direct imports does not cover its
  # dependencies' dependencies.
  env -u PYTHONPATH -u PYTHONHOME "$py" -m pip check >"$log.check" 2>&1 || true
  # Two very different findings come out of pip check, and only one has ever
  # broken anything here:
  #   "X requires Y, which is not installed."   -> a MISSING package. This is
  #        what absl, packaging and typing-extensions were, and each one killed
  #        a job. Fatal.
  #   "X has requirement Y>=2, but you have Y 1." -> a version floor not met.
  #        polyfunenv carries four of these (xarray/packaging,
  #        scikit-learn/joblib, scikit-learn/threadpoolctl, pandas-plink/cffi)
  #        and nevertheless produced both PolyFun priors end to end, so
  #        failing on them would block work that demonstrably runs. Warn.
  if grep -q 'which is not installed' "$log.check" 2>/dev/null; then
    grep 'which is not installed' "$log.check" | sed 's/^/    /' >&2
    bail "$n: a required package is not installed - see $log.check"
  fi
  if [ -s "$log.check" ] && ! grep -q '^No broken requirements found' "$log.check"; then
    echo "    pip check: version floors not met (not fatal; this env has been"
    echo "               exercised as-is). See $log.check:"
    sed 's/^/      /' "$log.check"
  else
    echo "    pip check: dependencies satisfied"
  fi
}

say "python environments"
# absl-py, tqdm and matplotlib are imported by the FB entry points, not just
# by its plotting: beatrice_joint.py does `from absl import app/flags` at
# module scope, so without absl-py the run dies immediately with
# ModuleNotFoundError before any work starts. Found by scanning every import
# in the FB tree rather than one error at a time. LassoNet needs no package --
# LassoNetPrior is FB's own class.
mkenv appenv    numpy pandas pyarrow scipy torch absl-py tqdm matplotlib packaging

say "verifying appenv can import what fb_xregion needs"
# beatrice_joint.py and its imports, checked as a real import rather than as a
# package list. Two separate runs have died here: absl (not installed at all)
# and packaging (a matplotlib dependency that went missing). Failing in setup
# costs seconds; failing in 03 costs a queue wait.
# PYTHONPATH cleared so this tests what a COMPUTE NODE sees: only the venv.
# Importing successfully is not enough -- each module must resolve INSIDE the
# venv. `packaging` imported fine on the login node from a module's
# site-packages and was absent in the job.
env -u PYTHONPATH -u PYTHONHOME "$APP_ROOT/appenv/bin/python" - "$APP_ROOT/appenv" <<'PYCHK' \
  || bail "appenv cannot import FB's dependencies from inside the venv"
import importlib, sys, os
venv = os.path.abspath(sys.argv[1])
mods = ["absl.app", "absl.flags", "matplotlib", "matplotlib.pyplot",
        "tqdm.auto", "torch", "torch.nn", "torch.optim", "numpy", "pandas",
        "pyarrow.parquet", "scipy.sparse", "packaging.version"]
failed, outside = [], []
for m in mods:
    try:
        mod = importlib.import_module(m)
    except Exception as e:
        failed.append(f"{m}: {type(e).__name__}: {e}")
        continue
    f = getattr(mod, "__file__", None)
    if f and not os.path.abspath(f).startswith(venv):
        outside.append(f"{m} -> {f}")
if failed:
    print("  FAILED to import:")
    for b in failed: print("   ", b)
if outside:
    print("  resolved OUTSIDE the venv (will not exist on a compute node):")
    for b in outside: print("   ", b)
if failed or outside:
    sys.exit(1)
import torch, numpy
print(f"  all {len(mods)} imports OK and inside the venv"
      f"  (torch {torch.__version__}, numpy {numpy.__version__})")
PYCHK

# Funmap pins these exactly; do not relax them.
# typing-extensions: funmapenv picked up rich-click transitively, and
# rich-click 1.7.0 declares typing-extensions, which was not installed. The
# new `pip check` in mkenv found this on its first run -- a latent break that
# would have surfaced in 04_funmap rather than here. Satisfying the declared
# dependency is safer than removing a package something pulled in.
mkenv funmapenv "numpy==1.26.2" "scipy==1.11.4" "pandas==2.1.3" "matplotlib==3.8.2" \
                typing-extensions
# PolyFun pins, each found the hard way:
#   pyarrow<=17   newer pyarrow requires numpy>=2, which conflicts with the rest
#   rpy2 3.5.x    3.6 turned numpy2ri.activate() into a raised DeprecationWarning
#                 and finemapper.py calls it
#   bitarray<3    3.x changed bitarray.decode() to return an iterator; PolyFun's
#                 plink reader passes the result straight to np.array(dtype=float64)
#                 and dies with "float() argument must be ... not decodeiterator"
mkenv polyfunenv "numpy<2" scipy pandas "pyarrow==17.0.0" scikit-learn pandas-plink \
                 "bitarray==2.9.2" networkx packaging tqdm "rpy2==3.5.17"

say "susieR 0.11.92 for PolyFun"
RLIB="$APP_ROOT/Rlib_polyfun"; mkdir -p "$RLIB"
if [ -d "$RLIB/susieR" ] && [ -d "$RLIB/Ckmeans.1d.dp" ]; then
  echo "  already installed"
else
  # Install from a CRAN SNAPSHOT contemporary with the R module, not from
  # today's CRAN. Two problems otherwise, both seen in practice:
  #   - the R/4.3.2 module ships no recommended packages, so Matrix must be
  #     installed, but current Matrix (>= 1.7) requires R >= 4.4 and is
  #     therefore "not available" on 4.3.2. irlba, then mixsqp, then susieR
  #     all fail behind it.
  #   - today's CRAN resolves a 2022 package's dependencies to 2026 versions
  #     (ggplot2 4.x, S7, RcppArmadillo 15.x), which is a lot of compiling for
  #     combinations nobody has tested together.
  # A snapshot dated just after R 4.3.2's release gives a dependency set that
  # was current when both R 4.3.2 and susieR 0.11.92 were in use.
  RSNAP="${R_CRAN_SNAPSHOT:-https://packagemanager.posit.co/cran/2023-11-01}"
  echo "  using CRAN snapshot: $RSNAP"

  # Safety net for the snapshot's C++11-era packages. R picks the compiler
  # standard from each package's own CXX_STD, and a package pinning CXX11
  # cannot build against the Armadillo that current RcppArmadillo bundles
  # (Armadillo >= 12 hard-errors: "C++14 compiler required"). Raising what
  # CXX11 resolves to keeps such a package buildable; gnu++14 is backward
  # compatible with gnu++11 for anything these packages do.
  MKV="$APP_ROOT/Rmakevars_polyfun"
  printf 'CXX11STD = -std=gnu++14\n' > "$MKV"
  export R_MAKEVARS_USER="$MKV"
  echo "  R_MAKEVARS_USER=$MKV (CXX11STD -> gnu++14)"

  Rscript --vanilla -e "
    .libPaths(c('$RLIB', .libPaths()))
    options(repos = c(CRAN = '$RSNAP'))
    # Matrix explicitly and first: it is the root of the failure chain.
    for (p in c('Matrix','irlba','matrixStats','reshape'))
      if (!requireNamespace(p, quietly=TRUE))
        install.packages(p, lib='$RLIB')
    stopifnot(requireNamespace('Matrix', quietly=TRUE))
    cat('  Matrix', as.character(packageVersion('Matrix')), 'available\n')

    # mixsqp is the ONE package taken from current CRAN rather than the
    # snapshot, and it is deliberate. The snapshot's 0.3-48 declares
    # CXX_STD = CXX11 in src/Makevars, so it compiles with -std=gnu++11 and
    # the bundled Armadillo 12 headers reject it outright. Upstream 0.3-54
    # drops the CXX_STD pin for exactly this reason and is otherwise
    # identical: same Depends, same Imports (utils, stats, irlba,
    # Rcpp >= 0.12.15), same LinkingTo. susieR 0.11.92 imports mixsqp with
    # no version constraint, so either satisfies it.
    if (!requireNamespace('mixsqp', quietly=TRUE))
      install.packages('mixsqp', lib='$RLIB',
                       repos = c(CRAN = 'https://cloud.r-project.org'))
    stopifnot(requireNamespace('mixsqp', quietly=TRUE))
    cat('  mixsqp', as.character(packageVersion('mixsqp')), 'available\n')
    # Ckmeans.1d.dp: PolyFun's DEFAULT SNP-binning step calls it through rpy2
    # (polyfun.py partition_snps_Ckmedian -> importr('Ckmeans.1d.dp')), and
    # --compute-h2-L2 dies there, AFTER the ridge regression has already run
    # and written snpvar_ridge, with
    #   'Could not load the R package Ckmeans.1d.dp. Either install it or
    #    rerun PolyFun with --skip-Ckmedian'
    # The --skip-Ckmedian fallback substitutes sklearn KMeans for the exact 1-D
    # dynamic-programming k-median, i.e. a different algorithm from the one the
    # paper describes. Install the package rather than take the fallback.
    if (!requireNamespace('Ckmeans.1d.dp', quietly=TRUE))
      install.packages('Ckmeans.1d.dp', lib='$RLIB')
    stopifnot(requireNamespace('Ckmeans.1d.dp', quietly=TRUE))
    cat('  Ckmeans.1d.dp', as.character(packageVersion('Ckmeans.1d.dp')), 'available\n')
    if (!requireNamespace('remotes', quietly=TRUE))
      install.packages('remotes', lib='$RLIB')
    remotes::install_github('stephenslab/susieR@v0.11.92', lib='$RLIB', upgrade='never')
    v <- as.character(packageVersion('susieR', lib.loc='$RLIB'))
    stopifnot(identical(v, '0.11.92'))
    cat('susieR', v, 'installed\n')"
fi

say "LDSC reference data for computing the PolyFun prior (0.30 GB)"
Z=https://zenodo.org/records/10515792/files
if [ -d "$APP_ROOT/1000G_EUR_Phase3_plink" ]; then echo "  plink panel already present"; else
  curl -fL --retry 5 -o "$APP_ROOT/plink.tgz" "$Z/1000G_Phase3_plinkfiles.tgz?download=1"
  tar xzf "$APP_ROOT/plink.tgz" -C "$APP_ROOT" && rm -f "$APP_ROOT/plink.tgz"
fi
if [ -d "$APP_ROOT/1000G_Phase3_weights_hm3_no_MHC" ]; then echo "  weights already present"; else
  curl -fL --retry 5 -o "$APP_ROOT/w.tgz" "$Z/1000G_Phase3_weights_hm3_no_MHC.tgz?download=1"
  tar xzf "$APP_ROOT/w.tgz" -C "$APP_ROOT" && rm -f "$APP_ROOT/w.tgz"
fi
echo "  1000G EUR panel and HapMap3 weights ready"

say "environment file"
cat > "$APP_ROOT/env.sh" <<ENVEOF
# source this before submitting
export APP_ROOT="$APP_ROOT"
# Modules the COMPUTE nodes need. The venvs below were built against the
# python active when setup.sh ran, so that python must be loadable there too.
export APP_MODULES="${APP_MODULES:-}"
# awscli separately: it is typically built against a different toolchain than
# Python and the two may refuse to load together. Only stage 01 needs it.
export AWS_MODULE="${AWS_MODULE:-}"
export N_LIKELIHOOD=373220
export APP_VENV="$APP_ROOT/appenv"
export FUNMAP_VENV="$APP_ROOT/funmapenv"
export FUNMAP_DIR="$APP_ROOT/Funmap"
export POLYFUN_ENV="$APP_ROOT/polyfunenv"
export POLYFUN_DIR="$APP_ROOT/polyfun"
export POLYFUN_RLIB="$RLIB"
  # The R library rpy2's EMBEDDED R must see. Its .libPaths() is otherwise
  # just the module's own, so a package installed into Rlib_polyfun is
  # invisible to it even though Rscript finds it. Each stage also sets this
  # via lib.sh's rpy2_env, so this is for interactive use.
  # Prepend only if absent: env.sh is sourced repeatedly in a shell, and an
  # unconditional prepend grows R_LIBS by one copy each time (observed at 5x).
  # R dedupes it so .libPaths() stays correct, but an unbounded variable is a
  # false lead waiting to be chased. lib.sh's rpy2_env guards the same way.
  case ":\${R_LIBS:-}:" in
    *":$RLIB:"*) ;;
    *) export R_LIBS="$RLIB\${R_LIBS:+:\$R_LIBS}" ;;
  esac
export R_HOME="\${R_HOME:-\$(R RHOME 2>/dev/null)}"
export RPY2_CFFI_MODE=ABI
# functional-beatrice: main at d0fb022 or later. 03_fb_xregion.pbs verifies this.
export FB_DIR="\${FB_DIR:-$APP_ROOT/functional-beatrice}"
ENVEOF
echo "  wrote $APP_ROOT/env.sh"

say "next"
cat <<'NEXT'
  1. Put functional-beatrice at $FB_DIR (main at d0fb022 or later):
       git clone https://github.com/lucassiu14/functional-beatrice.git "$FB_DIR"
  2. source "$APP_ROOT/env.sh"
  3. bash submit_all.sh
NEXT
