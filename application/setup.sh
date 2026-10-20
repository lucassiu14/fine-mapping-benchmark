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

say "pipeline scripts"
# lib.sh in particular: every stage does `source ./lib.sh` relative to
# APP_ROOT, so without it nothing runs at all.
if [ "$(cd "$HERE" && pwd)" != "$(cd "$APP_ROOT" && pwd)" ]; then
  cp -f "$HERE"/*.py "$HERE"/*.pbs "$HERE"/*.sh "$APP_ROOT"/
else
  echo "  (running from APP_ROOT; scripts already in place)"
fi
cp -f "$HERE"/blocks_scored.csv "$HERE"/blocks.csv "$APP_ROOT"/ 2>/dev/null || true
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

mkenv () {  # name, then pip args
  local n="$1"; shift
  if [ -x "$APP_ROOT/$n/bin/python" ]; then echo "  $n already built"; return; fi
  "$PY_BOOT" -m venv "$APP_ROOT/$n"
  "$APP_ROOT/$n/bin/python" -m pip install -q --upgrade pip
  "$APP_ROOT/$n/bin/python" -m pip install -q "$@"
  echo "  $n built"
}

say "python environments"
mkenv appenv    numpy pandas pyarrow scipy torch
# Funmap pins these exactly; do not relax them.
mkenv funmapenv "numpy==1.26.2" "scipy==1.11.4" "pandas==2.1.3" "matplotlib==3.8.2"
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
if [ -d "$RLIB/susieR" ]; then echo "  already installed"; else
  Rscript --vanilla -e "
    .libPaths('$RLIB')
    if (!requireNamespace('remotes', quietly=TRUE))
      install.packages('remotes', lib='$RLIB', repos='https://cloud.r-project.org')
    remotes::install_github('stephenslab/susieR@v0.11.92', lib='$RLIB', upgrade='never')
    cat('susieR', as.character(packageVersion('susieR', lib.loc='$RLIB')), 'installed\n')"
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
export N_LIKELIHOOD=373220
export APP_VENV="$APP_ROOT/appenv"
export FUNMAP_VENV="$APP_ROOT/funmapenv"
export FUNMAP_DIR="$APP_ROOT/Funmap"
export POLYFUN_ENV="$APP_ROOT/polyfunenv"
export POLYFUN_DIR="$APP_ROOT/polyfun"
export POLYFUN_RLIB="$RLIB"
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
