# Fine-Mapping Benchmarking Framework

Benchmarking framework for evaluating statistical fine-mapping methods using
simulated genetic data. Supports 8 methods, flexible simulation parameters,
automatic evaluation (AUPRC, credible set metrics, PIP calibration), and
multi-page PDF plots with stratified results.

## Requirements

- R >= 4.1.0
- [renv](https://rstudio.github.io/renv/) (installed automatically with the project)

For methods that rely on external tools, see [Method-specific setup](#method-specific-setup) below.

## Installation

```bash
git clone https://github.com/lucassiu14/fine-mapping-benchmark.git
cd fine-mapping-benchmark
```

Open R in the project directory and restore all R package dependencies:

```r
renv::restore()
```

This installs everything needed to run SuSiE, SuSiE-inf, ABF, and CARMA immediately. The other methods require additional steps described below.

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
sim <- run_simulation(
  n_regions     = 5,
  n             = 500,
  p             = 200,
  n_iter        = 10,
  S             = c(1, 2, 3),      # number of causal variants per region
  phi           = c(0.2, 0.5),     # per-causal heritability
  model         = "sparse",
  annotations   = "binary",
  n_annotations = 3,
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

Or run the end-to-end test pipeline (all 6 methods, 3 regions, save + plot):

```bash
Rscript scripts/test_pipeline.R
```

Output is written to `results/test_run/`.

## Supported methods

| Method | Type | Dependencies |
|---|---|---|
| **SuSiE** | R package | None (installed via renv) |
| **SuSiE-inf** | R package | None (installed via renv) |
| **ABF** | R (built-in) | None |
| **CARMA** | R package | None (installed via renv) |
| **FINEMAP** | C++ binary | Auto-downloaded by `setup_finemap()` |
| **PAINTOR** | C++ binary | `conda install -c bioconda paintor` |
| **Funmap** | Python package | Python + Funmap (see below) |
| **BEATRICE** | Python script | Python + BEATRICE repo (see below) |

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

### Funmap

Funmap is a Python package called from R via reticulate.

```bash
git clone https://github.com/LeeHITsz/Funmap.git
cd Funmap
pip install -r requirements.txt
pip install .
```

Then in R:

```r
source("R/wrappers/funmap.R")
setup_funmap(python = "/path/to/python3")   # pass your Python path
```

### BEATRICE

BEATRICE is a Python script. Clone the repository and create its conda environment:

```bash
git clone https://github.com/sayangsep/Beatrice-Finemapping
cd Beatrice-Finemapping
conda env create -f conda_environment.yml
conda activate beatrice
```

Then in R:

```r
source("R/wrappers/beatrice.R")
setup_beatrice(
  beatrice_dir = "~/Beatrice-Finemapping",
  python       = "~/anaconda3/envs/beatrice/bin/python"
)
```

Pass `beatrice_dir` and `python` via `method_args` in `run_methods()`:

```r
results <- run_methods(
  simulation  = sim,
  methods     = "beatrice",
  method_args = list(
    beatrice = list(
      beatrice_dir = "~/Beatrice-Finemapping",
      python       = "~/anaconda3/envs/beatrice/bin/python"
    )
  )
)
```

## Project structure

```
fine-mapping-benchmark/
├── R/
│   ├── simulate_genotypes.R    # Genotype simulation (sim1000G)
│   ├── simulate_phenotypes.R   # Phenotype simulation (sparse / infinitesimal)
│   ├── run_simulation.R        # Orchestration: runs simulation over scenarios
│   ├── run_methods.R           # Runs all methods on a simulation object
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
├── scripts/
│   ├── test_pipeline.R         # End-to-end pipeline test (all methods)
│   └── test_evaluate.R         # Unit tests for evaluation module (127 tests)
├── docs/
│   ├── methods.md              # Method descriptions and wrapper API
│   └── evaluation.md           # Evaluation metrics: formulas and implementation details
├── renv.lock                   # R package lockfile (use renv::restore())
└── README.md
```

Results and figures are written to `results/` and `figures/` (both gitignored).

## Simulation parameters

`run_simulation()` accepts:

| Parameter | Description |
|---|---|
| `n_regions` | Number of independent genomic regions |
| `n` | Sample size per region |
| `p` | Number of variants per region (capped at 500 for bundled VCF data) |
| `n_iter` | Number of simulation replicates per scenario |
| `S` | Vector of causal variant counts (one scenario per value) |
| `phi` | Vector of per-causal heritability values (crossed with S) |
| `model` | `"sparse"` or `"infinitesimal"` |
| `annotations` | `"binary"` or `"none"` |
| `n_annotations` | Number of binary annotation columns |
| `seed` | Random seed |
| `save` | Write simulation RDS to `output_dir` if `TRUE` |

## Evaluation output

`evaluate_methods()` returns metrics for each method, globally and stratified by S, phi, and p_causal (sparse/infinitesimal). Key metrics:

- **AUPRC** — area under the precision-recall curve (PIPs vs causal truth)
- **CS coverage** — proportion of credible sets containing at least one causal variant
- **CS power** — proportion of causal variants captured by any credible set
- **Median CS size** — median number of variants per credible set
- **PIP calibration** — binned mean PIP vs observed fraction causal

Standard errors are computed across replicates (`n_iter`) and stored alongside each metric. With `save = TRUE`, results are written as `evaluation.rds` and `evaluation_summary.csv`.

See [`docs/evaluation.md`](docs/evaluation.md) for full details.

## Adding a new method

1. Create `R/wrappers/wrap_mymethod.R` with a `run_mymethod_region(region, params)` function.
2. Return a named list with at minimum: `pip` (numeric vector, length p), `credible_sets` (list of integer vectors), `method`, `runtime_seconds`.
3. Register the method in `R/run_methods.R` in the `METHOD_REGISTRY` list.
4. Add a `setup_mymethod()` function if external dependencies are required.

See [`docs/methods.md`](docs/methods.md) for the full wrapper API specification.
