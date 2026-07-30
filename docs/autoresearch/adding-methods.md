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

## Installation on the HPC

Run from a login node with the toolchain modules loaded. `~/tools` is the
established location (`TOOLS_ROOT` in `scripts/hpc/run_benchmark_job.R`).

```bash
module load R/4.5.2-gfbf-2025b
module load GSL/2.8-GCC-14.3.0
module load Python/3.12.3-GCCcore-13.3.0
mkdir -p ~/tools && cd ~/tools
```

### 1. CARMA — why it failed before

CARMA is an R package (`ZikunY/CARMA`) built on Rcpp/RcppArmadillo/RcppGSL. It
needs **GSL headers** at compile time, which is almost certainly why the
Iteration 001 install failed silently and left the method at 100% NA. Load the
GSL module *before* installing:

```bash
module load GSL/2.8-GCC-14.3.0
Rscript -e 'remotes::install_github("ZikunY/CARMA")'
Rscript -e 'library(CARMA); cat("CARMA", as.character(packageVersion("CARMA")), "OK\n")'
```

If it still fails, capture the real error (the default install is quiet about
compile failures):

```bash
Rscript -e 'remotes::install_github("ZikunY/CARMA", quiet = FALSE)' 2>&1 | tail -40
```

Common causes: missing `gsl/gsl_*.h` (GSL module not loaded), or an
RcppArmadillo C++ standard mismatch — the same class of problem as FiniMOM
below, fixed the same way.

### 2. FiniMOM — the C++ standard trap (verified 2026-07-29)

A plain `remotes::install_github("vkarhune/finimom")` **fails**:

```
error: "*** C++14 compiler required; enable C++14 mode in your compiler"
```

The package's own `src/Makevars` pins `CXX_STD = CXX11`, while current
RcppArmadillo requires C++14 or later. **Setting `CXX_STD` in `~/.R/Makevars`
does not help — the package's own Makevars wins.** Patch the source:

```bash
cd ~/tools
git clone --depth 1 https://github.com/vkarhune/finimom.git
sed -i 's/^CXX_STD = CXX11/CXX_STD = CXX17/' finimom/src/Makevars
Rscript -e 'install.packages("~/tools/finimom", repos = NULL, type = "source")'
Rscript -e 'library(finimom); cat("finimom", as.character(packageVersion("finimom")), "OK\n")'
```

Verified working locally with this fix (finimom 0.2.0).

### 3. CAVIAR

```bash
cd ~/tools
git clone https://github.com/fhormoz/caviar.git
cd caviar/CAVIAR-C++ && make
./CAVIAR 2>&1 | head -5      # prints usage if the build succeeded
```

The worker expects the binary at `~/tools/caviar/CAVIAR-C++/CAVIAR`; override
with `METHOD_ARGS$caviar$caviar_path` if you put it elsewhere.

### 4. FINEMAP-inf

```bash
cd ~/tools
git clone https://github.com/FinucaneLab/fine-mapping-inf.git
cd fine-mapping-inf
~/tools/fmpy-venv/bin/python -m pip install -r requirements.txt
~/tools/fmpy-venv/bin/python run_fine_mapping.py -h | head -20
```

The worker expects `~/tools/fine-mapping-inf` and calls it through the existing
`py-venv-runner.sh`.

---

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
