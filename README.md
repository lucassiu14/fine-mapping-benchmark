# Fine-Mapping Benchmarking Framework

Benchmarking framework for evaluating statistical fine-mapping methods using
simulated genetic data. Supports 8 methods, flexible simulation parameters,
automatic evaluation (AUPRC, credible set metrics, PIP calibration), and
multi-page PDF plots with stratified results.

## Requirements

- R >= 4.1.0
- [renv](https://rstudio.github.io/renv/) (installed automatically with the project)
- [htslib](https://github.com/samtools/htslib) for `tabix` and `bgzip` (required by `scripts/prepare_vcfs.R`)

## Installation

```bash
git clone https://github.com/lucassiu14/fine-mapping-benchmark.git
cd fine-mapping-benchmark
```

Open R in the project directory and restore all R package dependencies:

```r
renv::restore()
```

> **Note:** Two packages (CARMA and susieR) install from GitHub. If you hit
> rate limits, set a GitHub personal access token first:
> `Sys.setenv(GITHUB_PAT = "your_token_here")` — or create one at
> github.com/settings/tokens (no scopes needed for public repos).

This installs everything needed to run SuSiE, SuSiE-inf, ABF, and CARMA
immediately. The other methods require additional steps described below.

## One-time reference data setup

The genotype simulator draws from real 1000 Genomes Phase 3 haplotypes.
Download 50 diverse 300 kb regions (one per genomic window in `data/regions.csv`,
covering all 22 autosomes):

```bash
Rscript scripts/prepare_vcfs.R
```

This streams each window from the 1000 Genomes EBI FTP via tabix — no
whole-chromosome downloads. Total download is ~150 MB. Files are saved to
`data/vcf/` and `data/genetic_maps/` (both gitignored).

> **Requires:** `tabix` and `bgzip` from htslib.
> Install with `brew install htslib` (macOS) or `conda install -c bioconda htslib`.

## Quick start

```r
source("R/simulate_genotypes.R")
source("R/simulate_phenotypes.R")
source("R/run_simulation.R")
source("R/run_methods.R")
source("R/evaluate.R")
source("R/plot_results.R")
source("R/wrappers/susie.R")
source("R/wrappers/abf.R")

# 1. Simulate genotypes + phenotypes
#    vcf_dir randomly samples n_regions from the 50 downloaded regions.
#    Set seed for reproducibility.
sim <- run_simulation(
  n_regions     = 20,
  n             = 500,
  p             = 400,
  n_iter        = 20,
  S             = c(1, 2, 3, 5),
  phi           = c(0.1, 0.2, 0.4),
  model         = "sparse",
  annotations   = "binary",
  n_annotations = 3,
  vcf_dir       = "data/vcf",
  seed          = 42
)

# 2. Run fine-mapping methods
results <- run_methods(
  simulation  = sim,
  methods     = c("susie", "abf"),
  method_args = list(
    susie = list(L = 10, coverage = 0.95),
    abf   = list(prior_variance = 0.04)
  )
)

# 3. Evaluate
eval_out <- evaluate_methods(sim, results, save = TRUE, output_dir = "results/run1")

# 4. Plot
plot_results(eval_out, output_file = "results/run1/results.pdf")
```

> **Without `prepare_vcfs.R`:** omit `vcf_dir` and the simulator falls back to
> the small bundled VCF (one chr4 region, ~500 SNPs). Useful for quick tests.

## Supported methods

| Method | Type | Dependencies |
|---|---|---|
| **SuSiE** | R package | None (installed via renv) |
| **SuSiE-inf** | R package | None (installed via renv) |
| **ABF** | R (built-in) | None |
| **CARMA** | R package | None (installed via renv) |
| **FINEMAP** | C++ binary | Auto-downloaded by `setup_finemap()` |
| **PAINTOR** | C++ binary | `conda install -c bioconda paintor` |
| **Funmap** | Python package | `conda env create -f environment.yml` |
| **BEATRICE** | Python script | `conda env create -f environment.yml` + BEATRICE repo |

Methods that fail (binary not found, Python error, etc.) are skipped gracefully
and reported in the results summary. They do not crash the pipeline.

## Method-specific setup

### FINEMAP

The binary is downloaded automatically the first time you call `setup_finemap()`:

```r
source("R/wrappers/finemap.R")
fp <- setup_finemap()   # downloads to R user cache dir; returns path
```

To disable auto-download and install manually, download the binary for your OS
from [http://www.christianbenner.com](http://www.christianbenner.com) and put it on your PATH.
Note: there is no official Windows binary — use WSL on Windows.

### PAINTOR

```bash
conda install -c bioconda paintor
```

Then in R:

```r
source("R/wrappers/paintor.R")
pp <- setup_paintor()   # finds PAINTOR on PATH
```

### Funmap and BEATRICE (shared Python environment)

Both Funmap and BEATRICE require Python. A single conda environment covers both.
An `environment.yml` is included in the repo:

```bash
conda env create -f environment.yml
conda activate finemapping-python
```

> **Apple Silicon:** remove the `cpuonly` line from `environment.yml` before
> creating the environment — PyTorch has native arm64 support.

> **GPU:** also remove `cpuonly` if you have a CUDA-capable GPU.

BEATRICE additionally requires its own repository (a Python script, not a package):

```bash
git clone https://github.com/sayangsep/Beatrice-Finemapping ~/Beatrice-Finemapping
```

Then get the Python path for use in R:

```bash
conda run -n finemapping-python which python   # copy this path
```

Pass it via `method_args` in `run_methods()`:

```r
PYTHON <- "/path/to/envs/finemapping-python/bin/python"   # from above

results <- run_methods(
  simulation  = sim,
  methods     = c("funmap", "beatrice"),
  method_args = list(
    funmap   = list(python = PYTHON, L = 10),
    beatrice = list(python = PYTHON, beatrice_dir = "~/Beatrice-Finemapping")
  )
)
```

## Project structure

```
fine-mapping-benchmark/
├── R/
│   ├── simulate_genotypes.R    # Genotype simulation (sim1000G + 1000G haplotypes)
│   ├── simulate_phenotypes.R   # Phenotype simulation (sparse / infinitesimal)
│   ├── run_simulation.R        # Orchestration: simulates over a parameter grid
│   ├── run_methods.R           # Runs fine-mapping methods on a simulation object
│   ├── evaluate.R              # Computes AUPRC, CS metrics, PIP calibration
│   ├── plot_results.R          # Generates multi-page PDF plots
│   └── wrappers/               # One file per method
│       ├── susie.R
│       ├── susie_inf.R
│       ├── abf.R
│       ├── carma.R
│       ├── finemap.R
│       ├── paintor.R
│       ├── funmap.R
│       └── beatrice.R
├── data/
│   ├── regions.csv             # 50 genomic regions used for genotype simulation
│   ├── vcf/                    # Downloaded 1000G VCF files (prepare_vcfs.R)
│   └── genetic_maps/           # Cached HapMap GRCh37 genetic maps (auto-downloaded)
├── scripts/
│   ├── prepare_vcfs.R          # One-time download of 1000G reference VCFs
│   ├── test_pipeline.R         # End-to-end pipeline test (all methods)
│   └── test_evaluate.R         # Unit tests for evaluation module (125 tests)
├── docs/
│   ├── methods.md              # Method descriptions and wrapper API
│   ├── evaluation.md           # Evaluation metrics: formulas and implementation
│   ├── Benchmarking.pdf        # Benchmarking design document
│   └── simulation_documentation.pdf  # Simulation methodology documentation
├── environment.yml             # conda environment for Funmap + BEATRICE
├── renv.lock                   # R package lockfile (use renv::restore())
└── README.md
```

Results are written to `results/` (gitignored).

## Simulation parameters

`run_simulation()` accepts:

| Parameter | Description |
|---|---|
| `n_regions` | Number of independent genomic regions |
| `n` | Sample size per region |
| `p` | Number of variants per region |
| `n_iter` | Number of simulation replicates per scenario |
| `S` | Vector of causal variant counts (one scenario per value) |
| `phi` | Vector of per-causal heritability values (crossed with S) |
| `model` | `"sparse"` or `"sparse_inf"` |
| `annotations` | `"binary"`, `"continuous"`, or `"none"` |
| `n_annotations` | Number of annotation columns |
| `vcf_dir` | Directory of prepared 1000G VCF files (see `scripts/prepare_vcfs.R`) |
| `genetic_map_dir` | Cache directory for HapMap genetic maps (default: `data/genetic_maps`) |
| `seed` | Random seed |
| `save` | Write simulation RDS to `output_dir` if `TRUE` |

## Evaluation output

`evaluate_methods()` returns metrics for each method, globally and stratified
by S, phi, and p_causal (sparse_inf model). Key metrics:

- **AUPRC** — area under the precision-recall curve (PIPs vs causal truth)
- **CS coverage** — proportion of credible sets containing at least one causal variant
- **CS power** — proportion of causal variants captured by any credible set
- **Median CS size** — median number of variants per credible set
- **PIP calibration** — binned mean PIP vs observed fraction causal

Standard errors are computed across replicates (`n_iter`). With `save = TRUE`,
results are written as `evaluation.rds` and `evaluation_summary.csv`.

See [`docs/evaluation.md`](docs/evaluation.md) for full details.

## Adding a new method

1. Create `R/wrappers/mymethod.R` with a `run_mymethod_region(region, params)` function.
2. Return a named list with at minimum: `pip` (numeric vector, length p), `credible_sets` (list of integer vectors), `method`, `runtime_seconds`.
3. Register the method in `R/run_methods.R` in `.FM_REGISTRY`.
4. Add a `setup_mymethod()` function if external dependencies are required.

See [`docs/methods.md`](docs/methods.md) for the full wrapper API specification.
