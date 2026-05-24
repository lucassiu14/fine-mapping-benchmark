# Fine-Mapping Benchmarking Framework

Benchmarking framework for evaluating statistical fine-mapping methods using
simulated genetic data. Supports 12 methods, flexible simulation parameters,
automatic evaluation (AUPRC, credible set metrics, PIP calibration), and
multi-page PDF plots with stratified results.

## Requirements

- R >= 4.1.0
- [renv](https://rstudio.github.io/renv/) (installed automatically with the project)

Some methods and the reference-data setup have additional system requirements
(C++ binaries, Python, `htslib`) — see the relevant tier or section below.

## Installation

Installation is tiered. **You only need to install what you actually want to
use.** Most of the experiment — including 7 of the 9 methods, all evaluation,
and all plotting — works at Tier 1 with just R packages.

### Clone the repo

```bash
git clone https://github.com/lucassiu14/fine-mapping-benchmark.git
cd fine-mapping-benchmark
```

### Tier 1 — Pure R (~15 min)

Open R in the project directory and restore all R package dependencies:

```r
renv::restore()
```

> **Note:** CARMA and susieR install from GitHub. A GitHub personal access
> token is required (even for public repos, unauthenticated installs hit rate
> limits quickly). Set one before running `renv::restore()`:
>
> ```r
> Sys.setenv(GITHUB_PAT = "your_token_here")
> renv::restore()
> ```
>
> Create a token at github.com/settings/tokens — no scopes needed.
>
> `sim1000G` and `hapsim` are archived from CRAN and are fetched directly
> from the CRAN archive; this is handled automatically by the lockfile.

**Tier 1 unlocks:**

- Methods: **SuSiE**, **SuSiE-inf**, **ABF**, **CARMA**, **marginal_z**\*,
  **polyfun_oracle**\*, **polyfun_est**\*
- All of `evaluate_methods()` (AUPRC, CS metrics, PIP calibration)
- All of `plot_results()` (multi-page PDF output)
- All simulation (`simulate_genotypes()`, `simulate_phenotypes()`,
  `run_simulation()`, `simulate_gwfm_data()`)

\*Added in Phase 1 of the development roadmap.

### Tier 2 — Adds C++ binary methods (~30 min more)

#### FINEMAP

The binary is downloaded automatically the first time you call
`setup_finemap()`:

```r
source("R/wrappers/finemap.R")
fp <- setup_finemap()   # downloads to R user cache dir; returns path
```

To disable auto-download and install manually, download the binary for your OS
from [http://www.christianbenner.com](http://www.christianbenner.com) and put
it on your PATH. Note: there is no official Windows binary — use WSL on
Windows.

#### PAINTOR

```bash
conda install -c bioconda paintor
```

Then in R:

```r
source("R/wrappers/paintor.R")
pp <- setup_paintor()   # finds PAINTOR on PATH
```

### Tier 3 — Adds Python methods (~1 hour more)

Funmap, BEATRICE, and Functional BEATRICE all require Python. A single conda
environment covers all three; an `environment.yml` is included in the repo:

```bash
conda env create -f environment.yml
conda activate finemapping-python
```

> **Apple Silicon:** remove the `cpuonly` line from `environment.yml` before
> creating the environment — PyTorch has native arm64 support.

> **GPU:** also remove `cpuonly` if you have a CUDA-capable GPU.

BEATRICE additionally requires its own repository (a Python script, not a
package):

```bash
git clone https://github.com/sayangsep/Beatrice-Finemapping ~/Beatrice-Finemapping
```

**Functional BEATRICE** is bundled in this repo (`BEATRICE_annot_sparse/`) and
needs no separate clone — it uses the same conda env.

Then get the Python path for use in R:

```bash
conda run -n finemapping-python which python   # copy this path
```

Pass it via `method_args` in `run_methods()`:

```r
PYTHON <- "/path/to/envs/finemapping-python/bin/python"   # from above

results <- run_methods(
  simulation  = sim,
  methods     = c("funmap", "beatrice", "functional_beatrice"),
  method_args = list(
    funmap              = list(python = PYTHON, L = 10),
    beatrice            = list(python = PYTHON, beatrice_dir = "~/Beatrice-Finemapping"),
    functional_beatrice = list(python = PYTHON, beatrice_dir = "BEATRICE_annot_sparse")
  )
)
```

## Reference data (optional)

The default genotype simulator can draw from the bundled `sim1000G` example
VCF (one chr4 region, ~500 SNPs) with no extra setup — useful for quick
tests. For realistic LD across multiple regions, download the per-locus or
genome-wide reference VCFs.

**Requires:** `tabix` and `bgzip` from
[htslib](https://github.com/samtools/htslib).

Install with `brew install htslib` (macOS) or
`conda install -c bioconda htslib`.

Then either:

```bash
# 50 diverse 300 kb regions (per-locus benchmark; ~150 MB)
Rscript scripts/prepare_vcfs.R

# OR 128 genome-wide regions (genome-wide benchmark; ~400 MB)
Rscript scripts/prepare_gwfm_vcfs.R
```

Each script streams the requested windows from the 1000 Genomes EBI FTP via
tabix — no whole-chromosome downloads. Files are saved to `data/vcf/` (or
`data/gwfm_vcf/`) and `data/genetic_maps/` (all gitignored).

## Quick start

```r
source("R/utils.R")
source("R/simulate_genotypes.R")
source("R/simulate_phenotypes.R")
source("R/run_simulation.R")
source("R/run_methods.R")
source("R/evaluate.R")
source("R/plot_results.R")
source("R/wrappers/susie.R")
source("R/wrappers/abf.R")

# 1. Simulate genotypes + phenotypes across a parameter grid
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
  seed          = 42,
  save          = TRUE,
  output_dir    = "results/run1"
)

# 2. Run fine-mapping methods
results <- run_methods(
  simulation  = sim,
  methods     = c("susie", "abf"),
  method_args = list(
    susie = list(L = 10, coverage = 0.95),
    abf   = list(prior_variance = 0.04)
  ),
  save       = TRUE,
  output_dir = "results/run1"
)

# 3. Evaluate
eval_out <- evaluate_methods(
  sim, results,
  save       = TRUE,
  output_dir = "results/run1"
)

# 4. Plot
plot_results(eval_out, output_dir = "results/run1")
```

> **Without `prepare_vcfs.R`:** omit `vcf_dir` and the simulator falls back to
> the small bundled VCF (one chr4 region, ~500 SNPs). Useful for quick tests.

## Supported methods

| Method | Tier | Type | Dependencies |
|---|---|---|---|
| **SuSiE** | 1 | R package | None (installed via renv) |
| **SuSiE-inf** | 1 | R package | None (installed via renv) |
| **ABF** | 1 | R (built-in) | None |
| **CARMA** | 1 | R package | None (installed via renv) |
| **FINEMAP** | 2 | C++ binary | Auto-downloaded by `setup_finemap()` |
| **PAINTOR** | 2 | C++ binary | `conda install -c bioconda paintor` |
| **Funmap** | 3 | Python package | `conda env create -f environment.yml` |
| **BEATRICE** | 3 | Python script | `conda env create -f environment.yml` + BEATRICE repo |
| **Functional BEATRICE** | 3 | Python script (in this repo) | Same conda env as BEATRICE; code is bundled in `BEATRICE_annot_sparse/` |
| **marginal_z** | 1 | R (built-in) | None |
| **polyfun_oracle** | 1 | R package | susieR (already in renv) |
| **polyfun_est** | 1 | R package | susieR (already in renv) |

See [Installation](#installation) for what each tier requires.

Methods that fail (binary not found, Python error, etc.) are skipped gracefully
and reported in the results summary. They do not crash the pipeline.

## Project structure

```
fine-mapping-benchmark/
├── R/
│   ├── utils.R                 # Shared helpers (sourced first)
│   ├── simulate_genotypes.R    # Per-locus genotype simulation (sim1000G + 1000G haplotypes)
│   ├── simulate_phenotypes.R   # Per-locus phenotype simulation (sparse / sparse+inf)
│   ├── simulate_gwfm_data.R    # Genome-wide simulation (shared y across regions)
│   ├── run_simulation.R        # Orchestrates per-locus simulations over a parameter grid
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
│       ├── beatrice.R
│       └── functional_beatrice.R
├── BEATRICE_annot_sparse/      # Functional BEATRICE source (Python + training scripts)
├── data/
│   ├── regions.csv                  # 50 per-locus regions (prepare_vcfs.R)
│   ├── gwfm_regions.csv             # 128 genome-wide regions (prepare_gwfm_vcfs.R)
│   ├── gwfm_regions_ldetect_EUR.csv # ~1,703 LDetect EUR blocks (download_ldetect_regions.R)
│   ├── vcf/                         # Downloaded 1000G VCF slices for per-locus regions
│   └── genetic_maps/                # Cached HapMap GRCh37 maps (auto-downloaded)
├── scripts/
│   ├── prepare_vcfs.R               # Download VCFs for the 50 per-locus regions
│   ├── prepare_gwfm_vcfs.R          # Download VCFs for the 128 genome-wide regions
│   ├── download_ldetect_regions.R   # (Optional) fetch LDetect block partition + VCFs
│   ├── test_pipeline.R              # End-to-end pipeline test (all methods)
│   ├── test_evaluate.R              # Unit tests for evaluation module
│   ├── test_comprehensive.R         # Argument-level tests for all functions
│   └── hpc/                         # SLURM job array for the benchmark
│       ├── generate_params_grid.R   # Build params_grid.csv (40 jobs)
│       ├── params_grid.csv          # One row per HPC array task
│       ├── run_benchmark_job.R      # Worker script: runs one (model, p, annot) combination
│       ├── collect_results.R        # Combines per-job outputs into one data frame
│       ├── smoke_test.R             # Per-method smoke test before submitting
│       └── submit_benchmark.sh      # Submits the SLURM array
├── docs/
│   ├── methods.md                       # Method descriptions and wrapper API
│   ├── evaluation.md                    # Evaluation metrics: formulas and implementation
│   ├── gw_simulation_documentation.md   # Technical spec of genome-wide simulation
│   └── testing_report.md                # Auto-generated argument-level test report
├── environment.yml             # conda environment for Funmap + BEATRICE
├── renv.lock                   # R package lockfile (use renv::restore())
└── README.md
```

Results are written to `results/` (gitignored).

---

## API Reference

All six main functions share a consistent interface: every function accepts
`save`, `output_dir`, and `verbose` arguments. Setting `save = TRUE` writes
output to `output_dir` (which is created automatically if it does not exist).
Setting `verbose = FALSE` suppresses all progress messages.

---

### `simulate_genotypes()`

Simulates genotype matrices from 1000 Genomes haplotypes using sim1000G.
Returns one matrix per region with realistic LD structure.

```r
genotypes <- simulate_genotypes(
  n_regions       = 3,
  n               = 500,
  p               = 200,
  vcf_files       = NULL,
  min_maf         = 0.01,
  max_maf         = NA,
  standardise     = TRUE,
  genetic_map_dir = "data/genetic_maps",
  seed            = NULL,
  save            = FALSE,
  output_dir      = "results",
  verbose         = TRUE
)
```

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `n_regions` | integer ≥ 1 | `3` | Number of independent genomic regions to simulate |
| `n` | integer ≥ 1 | `500` | Number of unrelated individuals |
| `p` | integer or integer vector | `200` | Target number of SNPs per region. Scalar = same for all regions; vector of length `n_regions` = per-region targets. Capped at 500 when using the bundled VCF |
| `vcf_files` | character vector or `NULL` | `NULL` | VCF files to simulate from (one per region, or one path reused for all). `NULL` uses the bundled sim1000G example VCF |
| `min_maf` | numeric ∈ [0, 0.5] | `0.01` | Minimum minor allele frequency filter |
| `max_maf` | numeric or `NA` | `NA` | Maximum MAF filter; `NA` = no upper limit |
| `standardise` | logical | `TRUE` | If `TRUE`, standardise each column to mean 0, variance 1. If `FALSE`, return raw 0/1/2 dosages |
| `genetic_map_dir` | character or `NULL` | `NULL` | Directory for caching HapMap GRCh37 genetic maps. `NULL` uses `tempdir()` (re-downloaded each session) |
| `seed` | integer or `NULL` | `NULL` | Random seed for reproducibility |
| `save` | logical | `FALSE` | If `TRUE`, write the returned list as `genotypes_{n_regions}regions_n{n}_p{p}_{seed}.rds` inside `output_dir` |
| `output_dir` | character | `"results"` | Directory for saved output; created automatically if absent |
| `verbose` | logical | `TRUE` | Print progress messages |

**Returns:** A list of length `n_regions`. Each element contains:
`X` (standardised genotype matrix, n × p), `X_raw` (raw 0/1/2 matrix),
`n`, `p`, `maf`, `variant_ids`, `region_id`, `vcf_source`.

---

### `simulate_phenotypes()`

Takes the output of `simulate_genotypes()` and adds phenotypes, summary
statistics, LD matrices, and ground truth to each region.

```r
sim <- simulate_phenotypes(
  genotypes              = genotypes,
  S                      = 1,
  phi                    = 0.1,
  model                  = "sparse",
  p_causal               = 0.5,
  inf_model              = "beatrice",
  effect_distribution    = "normal",
  effect_variance        = 0.36,
  annotations            = "none",
  n_annotations          = 3,
  annotation_proportions = NULL,
  enrichment             = NULL,
  seed                   = NULL,
  save                   = FALSE,
  output_dir             = "results",
  verbose                = TRUE
)
```

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `genotypes` | list | — | Output from `simulate_genotypes()` |
| `S` | integer or integer vector | `1` | Number of causal variants per region. Scalar = same for all; vector of length `n_regions` = per-region |
| `phi` | numeric or numeric vector | `0.1` | Proportion of variance explained (PVE). Must be in (0, 1). Scalar or per-region vector |
| `model` | `"sparse"` or `"sparse_inf"` | `"sparse"` | Genetic architecture. `"sparse"` = causal effects only; `"sparse_inf"` = sparse + infinitesimal background |
| `p_causal` | numeric ∈ (0, 1] | `0.5` | Fraction of total genetic variance from the sparse component. Only used when `model = "sparse_inf"` |
| `inf_model` | `"beatrice"` or `"susie_inf"` | `"beatrice"` | Infinitesimal formulation: `"beatrice"` draws background effects from non-causal variants only; `"susie_inf"` uses all variants. Only used when `model = "sparse_inf"` |
| `effect_distribution` | `"normal"` or `"equal"` | `"normal"` | Distribution for causal effect sizes. `"normal"` draws from N(0, `effect_variance`); `"equal"` partitions variance equally |
| `effect_variance` | numeric > 0 | `0.36` | Variance of the normal effect size distribution (SD ≈ 0.6). Only used when `effect_distribution = "normal"` |
| `annotations` | `"none"`, `"binary"`, `"continuous"`, or matrix | `"none"` | Functional annotation mode. A user-supplied p × m matrix is also accepted |
| `n_annotations` | integer ≥ 1 | `3` | Number of annotation columns (for `"binary"` or `"continuous"`) |
| `annotation_proportions` | numeric, vector, or `NULL` | `NULL` | Proportion of SNPs with value 1 per binary annotation. `NULL` = random from Uniform(0.01, 0.30). Scalar or vector of length `n_annotations` |
| `enrichment` | numeric, vector, or `NULL` | `NULL` | Fold-enrichment of each annotation for causal variant selection. `NULL` = random from Uniform(2, 10). Scalar or vector of length `n_annotations` |
| `seed` | integer or `NULL` | `NULL` | Random seed |
| `save` | logical | `FALSE` | If `TRUE`, write the returned list as `phenotypes_{model}_S{S}_phi{phi}_{seed}.rds` inside `output_dir` |
| `output_dir` | character | `"results"` | Directory for saved output; created automatically if absent |
| `verbose` | logical | `TRUE` | Print progress messages |

**Returns:** The input `genotypes` list with additional fields per region:
`y`, `z`, `beta_hat`, `se`, `LD`, `annotations_matrix`, `truth`.
The `truth` sub-list records `causal_indices`, `causal_effects`, `beta_true`, `pve`, `S`, `phi`, `model`, and annotation settings.

---

### `run_simulation()`

Orchestrates a full benchmarking simulation. Simulates genotypes once, then
sweeps over all combinations of S, phi (and optionally p_causal), generating
`n_iter` independent replicates per combination.

```r
sim <- run_simulation(
  n_regions              = 3,
  n                      = 500,
  p                      = 200,
  n_iter                 = 5,
  S                      = c(1, 2, 3, 5),
  phi                    = c(0.1, 0.2, 0.4, 0.6),
  model                  = "sparse",
  p_causal               = c(0.1, 0.2, 0.4),
  inf_model              = "beatrice",
  effect_distribution    = "normal",
  effect_variance        = 0.36,
  annotations            = "none",
  n_annotations          = 3,
  annotation_proportions = NULL,
  enrichment             = NULL,
  vcf_dir                = NULL,
  vcf_files              = NULL,
  genetic_map_dir        = "data/genetic_maps",
  min_maf                = 0.01,
  max_maf                = NA,
  standardise            = TRUE,
  seed                   = NULL,
  save                   = FALSE,
  output_dir             = "results",
  verbose                = TRUE
)
```

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `n_regions` | integer ≥ 1 | `3` | Number of independent genomic regions |
| `n` | integer ≥ 1 | `500` | Number of individuals |
| `p` | integer or vector | `200` | Target SNPs per region |
| `n_iter` | integer ≥ 1 | `5` | Independent replicates per parameter combination |
| `S` | integer vector | `c(1,2,3,5)` | Causal-variant counts to sweep over |
| `phi` | numeric vector ∈ (0, 1) | `c(0.1,0.2,0.4,0.6)` | PVE values to sweep over |
| `model` | `"sparse"` or `"sparse_inf"` | `"sparse"` | Genetic architecture model |
| `p_causal` | numeric vector ∈ (0, 1] | `c(0.1,0.2,0.4)` | Sparse-component fractions to sweep. Only used when `model = "sparse_inf"` |
| `inf_model` | `"beatrice"` or `"susie_inf"` | `"beatrice"` | Infinitesimal formulation. Only used when `model = "sparse_inf"` |
| `effect_distribution` | `"normal"` or `"equal"` | `"normal"` | Effect size distribution |
| `effect_variance` | numeric > 0 | `0.36` | Variance for normal effect sizes |
| `annotations` | `"none"`, `"binary"`, `"continuous"`, or matrix | `"none"` | Functional annotation mode |
| `n_annotations` | integer ≥ 1 | `3` | Number of annotation columns |
| `annotation_proportions` | numeric, vector, or `NULL` | `NULL` | Proportion of 1s per binary annotation (scalar or per-annotation vector) |
| `enrichment` | numeric, vector, or `NULL` | `NULL` | Fold-enrichment per annotation for causal selection |
| `vcf_dir` | character or `NULL` | `NULL` | Directory of VCF files from `scripts/prepare_vcfs.R`. `n_regions` files are sampled at random (reproducibly if `seed` is set) |
| `vcf_files` | character vector or `NULL` | `NULL` | Explicit VCF paths; overrides `vcf_dir` when supplied |
| `genetic_map_dir` | character or `NULL` | `"data/genetic_maps"` | Cache directory for HapMap genetic maps |
| `min_maf` | numeric | `0.01` | Minimum MAF filter |
| `max_maf` | numeric or `NA` | `NA` | Maximum MAF filter |
| `standardise` | logical | `TRUE` | Standardise genotype columns |
| `seed` | integer or `NULL` | `NULL` | Master random seed |
| `save` | logical | `FALSE` | If `TRUE`, write the full result as an `.rds` file inside `output_dir` |
| `output_dir` | character | `"results"` | Output directory; created automatically if absent |
| `verbose` | logical | `TRUE` | Print progress messages |

**Returns:** A list with:
- `genotypes` — output of `simulate_genotypes()`, with pre-computed LD matrices
- `scenarios` — list of simulation scenarios, each containing `S`, `phi`, `p_causal`, `iter`, `model`, and `regions` (per-region phenotypes + truth)
- `params` — all simulation parameters for reproducibility

The total number of scenarios is `length(S) × length(phi) × n_iter` (sparse) or
`length(S) × length(phi) × length(p_causal) × n_iter` (sparse_inf).

---

### `run_methods()`

Applies one or more fine-mapping methods to every scenario × region combination
in a simulation object.

```r
results <- run_methods(
  simulation  = sim,
  methods     = "susie",
  method_args = list(),
  save        = FALSE,
  output_dir  = "results",
  verbose     = TRUE
)
```

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `simulation` | list | — | Output of `run_simulation()` |
| `methods` | character vector | `"susie"` | Methods to run (case-insensitive). Any combination of the supported methods below |
| `method_args` | named list | `list()` | Per-method argument overrides. Each entry is a named list passed to that method's region runner. Arguments not listed use the method's own defaults |
| `save` | logical | `FALSE` | If `TRUE`, write each method's results as `{method}.rds` and a `run_metadata.rds` file inside a sub-directory of `output_dir` |
| `output_dir` | character | `"results"` | Root output directory; sub-directory is named after simulation parameters |
| `verbose` | logical | `TRUE` | Print per-fit progress |

**Supported methods and their tuneable arguments (via `method_args`):**

| Method key | Key arguments |
|-----------|--------------|
| `"susie"` | `L` (int, default 10), `coverage` (0–1, default 0.95), `min_abs_corr` (default 0.5), `max_iter` (default 100), `estimate_residual_variance` (logical), `estimate_prior_variance` (logical), `prior_variance` (numeric) |
| `"susie_inf"` | `L` (int, default 10), `coverage` (0–1, default 0.95), `max_iter` (default 100) |
| `"abf"` | `prior_variance` (default 0.04), `coverage` (0–1, default 0.95) |
| `"carma"` | `rho.index` (default 0.95), `num.causal` (default 10) |
| `"finemap"` | `finemap_path` (path to binary), `n_causal` (default 5), `prior_std` (default 0.05) |
| `"paintor"` | `paintor_path` (path to binary), `max_causal` (default 2) |
| `"beatrice"` | `beatrice_dir` (path to repo), `python` (Python executable), `max_iter` (default 2000) |
| `"funmap"` | `python` (Python executable), `L` (default 10), `max_iter` (default 100) |

**Returns:** A named list, one entry per method plus `methods_run`, `simulation_params`, `run_timestamp`.
Each method entry contains `results` (flat list of per-fit outputs), `n_total`, `n_failed`, `total_runtime_seconds`.
Each per-fit result has `pip`, `credible_sets`, `method`, `input_type`, `params`, `runtime_seconds`, `additional`, plus scenario metadata (`scenario_id`, `region_id`, `S`, `phi`, `iter`).

---

### `evaluate_methods()`

Computes evaluation metrics for each method against the simulation ground truth.
Results are stratified globally and by S, phi, and p_causal (sparse_inf only).

```r
eval_out <- evaluate_methods(
  simulation     = sim,
  results        = results,
  pip_thresholds = seq(0, 1, by = 0.005),
  n_pip_cal_bins = 10L,
  save           = FALSE,
  output_dir     = "results",
  verbose        = TRUE
)
```

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `simulation` | list | — | Output of `run_simulation()` |
| `results` | list | — | Output of `run_methods()` |
| `pip_thresholds` | numeric vector | `seq(0, 1, by = 0.005)` | PIP thresholds at which to compute power and FDR for the precision-recall curve |
| `n_pip_cal_bins` | integer ≥ 1 | `10` | Number of equal-width bins for the PIP calibration plot |
| `save` | logical | `FALSE` | If `TRUE`, write `evaluation.rds` (full object) and `evaluation_summary.csv` (per-method global metrics) to `output_dir` |
| `output_dir` | character | `"results"` | Output directory; created automatically if absent |
| `verbose` | logical | `TRUE` | Print progress messages |

**Returns:** A named list (one entry per method) containing `global`, `by_S`, `by_phi`, and `by_p_causal` strata.
Each stratum contains:

| Field | Description |
|-------|-------------|
| `fdr_power_curve` | Data frame with columns `threshold`, `tp`, `fp`, `fn`, `fdr`, `power`, `precision`, `recall`, `power_se`, `precision_se` |
| `auprc` | Area under the precision-recall curve (±`auprc_se`) |
| `pip_calibration` | Data frame: binned mean PIP vs observed fraction causal (±`frac_causal_se`) |
| `cs_coverage` | Proportion of reported credible sets containing ≥1 true causal variant (±`cs_coverage_se`) |
| `cs_power` | Proportion of true causal variants captured by any credible set (±`cs_power_se`) |
| `cs_size_median` | Median number of variants per credible set (±SE) |
| `cs_size_mean` | Mean number of variants per credible set (±SE) |
| `n_cs_reported` | Total credible sets reported across all fits |
| `runtime_mean` | Mean runtime in seconds across successful fits (±`runtime_mean_se`) |
| `n_fits` | Total fits attempted |
| `n_failed` | Fits that errored |

Standard errors are computed across replicates (`n_iter`). SE fields are `NA` when `n_iter < 2`.

---

### `plot_results()`

Generates a multi-page PDF with precision-recall curves, PIP calibration plots,
and summary metric panels — globally and stratified by S, phi, and p_causal.

```r
plot_results(
  eval_out    = eval_out,
  output_file = NULL,
  output_dir  = "results",
  save        = TRUE,
  methods     = NULL,
  verbose     = TRUE
)
```

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `eval_out` | list | — | Output of `evaluate_methods()` |
| `output_file` | character or `NULL` | `NULL` | Full path for the PDF. When specified, takes precedence over `output_dir`. When `NULL`, the file is written as `evaluation.pdf` inside `output_dir` |
| `output_dir` | character | `"results"` | Directory in which to write `evaluation.pdf` when `output_file` is `NULL`. Created automatically if absent |
| `save` | logical | `TRUE` | If `FALSE`, skip writing the PDF entirely (useful for dry runs) |
| `methods` | character vector or `NULL` | `NULL` | Methods to include in the plots. `NULL` = all evaluated methods |
| `verbose` | logical | `TRUE` | Print section progress |

**Returns:** Invisibly returns the path of the PDF that was (or would have been) written.

**PDF sections:**
1. **Global** — PR curve (all methods), PIP calibration (faceted by method), summary table with ±SE
2. **By S** — PR grid, calibration grid, metric line plots vs S with error bars
3. **By phi** — same structure as By S
4. **By p_causal** — same structure (sparse_inf model only)

---

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

1. Create `R/wrappers/mymethod.R` with a `run_mymethod_region(region_geno, region_pheno, ...)` function.
2. Return a named list with at minimum: `pip` (numeric vector, length p), `credible_sets` (list of integer vectors), `method`, `input_type`, `params`, `runtime_seconds`, `additional`.
3. Register the method in `R/run_methods.R` in `.FM_REGISTRY`.
4. Add a `setup_mymethod()` function if external dependencies are required.

See [`docs/methods.md`](docs/methods.md) for the full wrapper API specification.

## Testing

Three test scripts are provided:

| Script | Coverage | Tests |
|--------|----------|-------|
| `scripts/test_comprehensive.R` | Every argument of every public function | 201 PASS / 2 SKIP |
| `scripts/test_evaluate.R` | Evaluation module unit tests | 125 PASS |
| `scripts/test_pipeline.R` | End-to-end pipeline (all methods) | Informational |

Run any test from the project root:

```bash
Rscript scripts/test_comprehensive.R
```

The 2 SKIPs in `test_comprehensive.R` are the FINEMAP and PAINTOR binary tests,
which require external binaries unavailable on Apple Silicon without additional setup
(see method-specific setup above).

A human-readable argument-level test report is generated automatically at
[`docs/testing_report.md`](docs/testing_report.md).

## Genome-wide simulation

`simulate_gwfm_data()` replaces the locus-based `run_simulation()` with a
genome-wide model: a single shared phenotype is generated across all regions
(as in a real GWAS), causal variants are assigned genome-wide via Bernoulli(π),
and per-region summary statistics are computed from the shared phenotype.

```r
source("R/utils.R")
source("R/simulate_genotypes.R")
source("R/simulate_phenotypes.R")
source("R/simulate_gwfm_data.R")

sim <- simulate_gwfm_data(
  n        = 2000,
  n_iter   = 3,
  pi       = c(1e-4, 1e-3),
  h2       = c(0.1, 0.3),
  regions  = "representative",   # 128 bundled regions; or "1" for one chromosome
  vcf_dir  = "data/gwfm_vcf",
  seed     = 1
)
```

The resulting object is compatible with `run_methods()` and `evaluate_methods()`.
See [`docs/gw_simulation_documentation.md`](docs/gw_simulation_documentation.md)
for the full statistical specification.

**Memory:** the simulator holds one dense p×p LD matrix per region in memory.
A warning is emitted when combined LD storage exceeds ~4 GB. For the full LDetect
partition (~1,700 blocks), run on a high-memory node or reduce `coverage`.

## Running on an HPC cluster

The `scripts/hpc/` directory contains a self-contained SLURM job array that
sweeps the benchmark across `(model, p, annotation)` combinations.

```bash
# 1. (Optional) regenerate the grid
Rscript scripts/hpc/generate_params_grid.R

# 2. Smoke-test all methods first
Rscript scripts/hpc/smoke_test.R

# 3. Submit the array (40 jobs)
bash scripts/hpc/submit_benchmark.sh

# 4. Combine results once jobs complete
Rscript scripts/hpc/collect_results.R
```

Edit `submit_benchmark.sh` for your cluster's partition, account, and conda
environment, and edit `run_benchmark_job.R` for the paths to BEATRICE / PAINTOR
on your system.
