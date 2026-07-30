# Adding methods: CARMA, CAVIAR, FiniMOM, FINEMAP-inf

Four summary-statistic methods added after the Iteration 003 deep dive. The
selection is not arbitrary: the deep dive found that **LD mis-specification, not
the annotation prior, is the binding constraint**, so the priority is methods
built for that, plus one that resolves an open question about the infinitesimal
model.

| method | why it was added |
|---|---|
| **CARMA** | Its defining feature is **LD-discrepancy outlier detection** — exactly the failure mode we identified. Already wrapped since Iteration 001 but **never installed**, so it has contributed 100% NA to every run so far. |
| **CAVIAR** | The original enumeration-based summary-stat method; a distinct algorithmic family from everything in the benchmark, and the ancestor of PAINTOR/FINEMAP. |
| **FiniMOM** | Non-local (product inverse-moment) prior + beta-binomial model-size prior, designed to separate **multiple causal variants in high LD**. Uniquely, it has an explicit `insampleLD` switch that maps directly onto our LD arm. |
| **FINEMAP-inf** | susie_inf trails plain susie in *every* stratum — but it inherits SuSiE's LD fragility (LD explains 56% of its reliability variance vs **<1%** of finemap's AUPRC variance). FINEMAP-inf puts the same infinitesimal extension on the **LD-robust finemap backbone**, separating "infinitesimal modelling doesn't help" from "the SuSiE backbone is the problem". |

All four are summary-statistic based (z-scores/betas + LD), so they are in scope.

---

## Installation on the HPC — VERIFIED RECIPE

Every step below was executed on Imperial's HPC on 2026-07-29. **Five separate
environment failures** were hit and fixed; none of them is documented upstream,
which is why this section is prescriptive rather than a link.

> **Do not run `module purge`.** This cluster has `LMOD_SYSTEM_DEFAULT_MODULES`
> empty, so `module purge` (and `module reset`) leaves an empty hierarchy in
> which `R/4.5.2-gfbf-2025b` cannot even be found. If that happens, start a
> fresh login shell.

```bash
module load R/4.5.2-gfbf-2025b GSL/2.8-GCC-14.3.0
cd ~/fine-mapping-benchmark          # ESSENTIAL for the R packages: this project
                                     # uses renv, and only here is the library
                                     # writable. Installing from ~/tools fails
                                     # with "lib = ... is not writable".
```

`R/4.5.2-gfbf-2025b` is on the **GCCcore/14.3.0** lineage, so `GSL/2.8-GCC-14.3.0`
is the matching GSL. Lmod prints a long "dependent module(s) not currently
loaded" warning when these load together — it is narrating module swaps and is
**benign**.

### 1. CARMA — needs GSL at build time

`LinkingTo: Rcpp, RcppArmadillo, RcppGSL`, so the GSL module must be loaded
*before* installing. This is why the Iteration 001 install failed and CARMA sat
at 100% NA for three iterations.

```bash
cd ~/fine-mapping-benchmark
Rscript -e 'remotes::install_github("ZikunY/CARMA")'
Rscript -e 'library(CARMA); cat("CARMA", as.character(packageVersion("CARMA")), "OK\n")'
```

`RcppGSL` does not need to remain installed afterwards — `LinkingTo` packages are
compile-time only.

### 2. FiniMOM — two traps

**(a) C++ standard.** The package pins `CXX_STD = CXX11` in `src/Makevars`, but
current RcppArmadillo needs C++14+. Setting `CXX_STD` in `~/.R/Makevars` does
**not** help — the package's own Makevars wins. **(b)** it needs `roptim`, which
is not pulled in automatically by a local source install.

```bash
cd ~/tools
git clone --depth 1 https://github.com/vkarhune/finimom.git
sed -i 's/^CXX_STD = CXX11/CXX_STD = CXX17/' finimom/src/Makevars
cd ~/fine-mapping-benchmark
Rscript -e 'install.packages("roptim")'
Rscript -e 'install.packages("~/tools/finimom", repos = NULL, type = "source")'
Rscript -e 'library(finimom); cat("finimom OK\n")'
```

### 3. CAVIAR — pre-C++17 code, and no generic BLAS

Two problems: its `struct data` / `size()` collide with `std::data` / `std::size`
once compiled as C++17 (GCC's default), and it links `-llapack -lblas`, which
this cluster does not provide — the gfbf toolchain ships **FlexiBLAS**. Editing
the Makefile is not enough (it sets `CC=g++` and compiles via `$(CC)`); link
explicitly:

```bash
cd ~/tools
git clone https://github.com/fhormoz/caviar.git
cd caviar/CAVIAR-C++
g++ -std=gnu++11 caviar.cpp PostCal.cpp Util.cpp TopKSNP.cpp \
    -I armadillo/include/ -DARMA_DONT_USE_WRAPPER \
    -L$EBROOTFLEXIBLAS/lib -lflexiblas \
    -L$EBROOTGSL/lib -lgsl -lgslcblas \
    -o CAVIAR
./CAVIAR 2>&1 | head -3      # usage text = success
```

### 4. FINEMAP-inf — three traps

There is **no `requirements.txt`**; it is two sub-packages installed separately.
The venv python also needs the Python module loaded (it is linked against the
module's libpython), the `finemapinf` build imports numpy so pip's **build
isolation must be disabled**, and `run_fine_mapping.py` imports `pkg_resources`,
which setuptools >= 81 removed.

```bash
module load Python/3.12.3-GCCcore-13.3.0     # else: libpython3.12.so.1.0 not found
PY=~/tools/fmpy-venv/bin/python
cd ~/tools && git clone https://github.com/FinucaneLab/fine-mapping-inf.git
$PY -m pip install setuptools wheel numpy pandas scipy bgzip
$PY -m pip install ~/tools/fine-mapping-inf/susieinf
$PY -m pip install --no-build-isolation ~/tools/fine-mapping-inf/finemapinf
```

`pkg_resources` is used only to log a version string (lines 106/113), so shim it
rather than downgrading setuptools in a venv shared with torch and optuna:

```bash
$PY - <<'EOF'
p = "/rds/general/user/<you>/home/tools/fine-mapping-inf/run_fine_mapping.py"
s = open(p).read()
if "_PkgShim" not in s:
    shim = """try:
    import pkg_resources
except ImportError:
    class _PkgShim:
        @staticmethod
        def require(name):
            class _V: version = "unknown"
            return [_V()]
    pkg_resources = _PkgShim()"""
    open(p, "w").write(s.replace("import pkg_resources", shim, 1))
EOF
$PY ~/tools/fine-mapping-inf/run_fine_mapping.py -h | head -5
```

> **Author caveat worth carrying into the analysis.** The FINEMAP-inf README
> states that it and SuSiE-inf "are developed for use in single-cohort
> fine-mapping with **in-sample LD**; for fine-mapping of meta-analyzed GWAS with
> reference panel LD, please see methods like SLALOM or **CARMA**". So
> FINEMAP-inf's n_ref arms are outside its intended use and must be reported with
> that flag — and the authors independently point at CARMA for the
> reference-panel regime, which is where our own deep dive landed.

## Verifying before you commit compute

**Do not launch an array until each method has been shown to recover a planted
causal.** A method that silently returns NAs looks exactly like a method that
performs badly, and that is how CARMA went unnoticed for three iterations.

```bash
Rscript scripts/analysis/verify_new_methods.R
```

It plants known causal variants in a synthetic locus and checks, per method,
that the output is non-NA, correctly sized, and puts high PIP on the truth.
Methods whose dependency is absent report **SKIP**, never PASS.

---

## Running them

Supplemental run against the cached Iteration 002 sims, as usual:

```bash
export FMB_SCRATCH=$EPHEMERAL/fmbench_iter004     # NEW root: the iter002 root
                                                  # already has evaluation_supp.rds
                                                  # for every scenario, so resume
                                                  # would skip everything
export FMB_METHODS="carma,caviar,finimom,finemap_inf,paintor"
export FMB_SCENARIOS_PER_TASK=25
bash scripts/hpc/submit_benchmark_pbs.sh
```

`paintor` is included because its cross-loci pooling fix (see
`iteration-003.md`) invalidates its Iterations 001–003 numbers — those measured
PAINTOR with enrichment estimated from a single locus, which is not how the
method works.

### Cost warning

CAVIAR and CARMA are both **combinatorial in the number of causal variants**.
CAVIAR is capped at `max_causal = 2` and CARMA at `num.causal = 5` for that
reason. On the 1000-SNP regions these will be the slowest methods in the
benchmark — fire an `ARRAY_RANGE=1-2` canary and check the per-scenario timing
in the log before releasing the full array.
