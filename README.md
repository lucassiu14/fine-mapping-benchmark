# Fine-Mapping Benchmarking Framework

Benchmarking framework for evaluating statistical fine-mapping methods using simulated genetic data. Designed to run on both personal computers and HPC clusters.

## Setup

### Prerequisites

- R (>= 4.1.0)
- Git

### Installation

```bash
git clone https://github.com/YOUR_USERNAME/fine-mapping-benchmark.git
cd fine-mapping-benchmark
```

Open R in the project directory and restore dependencies:

```r
renv::restore()
```

## Project Structure

```
fine-mapping-benchmark/
├── R/                    # Core functions
│   └── wrappers/         # Per-method wrapper functions
├── config/               # Simulation and method configuration
├── scripts/              # Entry-point scripts
├── results/              # [gitignored] Output from benchmark runs
├── figures/              # [gitignored] Generated figures
├── tests/                # Unit tests
├── DESCRIPTION           # Project metadata and dependencies
├── renv.lock             # Dependency lockfile
└── README.md
```

## Usage

```r
# 1. Set your compute backend
library(future)
plan(multisession, workers = 4)        # laptop
# plan(future.batchtools_slurm)        # HPC cluster

# 2. Run the benchmark
source("scripts/run_benchmark.R")
run_benchmark(
  config_file = "config/scenarios.yaml",
  n_replicates = 10,
  output_dir = "results/"
)

# 3. Generate figures
source("scripts/generate_figures.R")
generate_figures(results_dir = "results/", output_dir = "figures/")
```

## Adding a New Method

1. Create `R/wrappers/wrap_yourmethod.R` implementing `wrap_yourmethod(sumstats, params)`.
2. The function must return an `fm_result` object (see `R/utils.R`).
3. Register the method in `R/run_method.R`.
