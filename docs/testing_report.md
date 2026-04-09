# Fine-Mapping Benchmark — Function & Argument Test Report

**Generated:** 2026-04-09 19:47:52  
**R version:** 4.5.3  
**Platform:** aarch64-apple-darwin25.3.0  

**Results: 201 PASS / 0 FAIL / 2 SKIP**

---

## Table of contents

- [simulate_genotypes](#simulate-genotypes)
- [simulate_phenotypes](#simulate-phenotypes)
- [simulate_genotypes — save / output_dir](#simulate-genotypes-save-output-dir)
- [simulate_phenotypes — save / output_dir](#simulate-phenotypes-save-output-dir)
- [run_simulation](#run-simulation)
- [run_methods](#run-methods)
- [evaluate_methods](#evaluate-methods)
- [plot_results](#plot-results)
- [run_susie / run_susie_region — method arguments](#run-susie-run-susie-region-method-arguments)
- [run_abf / run_abf_region — method arguments](#run-abf-run-abf-region-method-arguments)
- [run_susie_inf / run_susie_inf_region — method arguments](#run-susie-inf-run-susie-inf-region-method-arguments)
- [run_carma / run_carma_region — method arguments](#run-carma-run-carma-region-method-arguments)
- [External wrappers — argument forwarding (FINEMAP, PAINTOR, BEATRICE, FUNMAP)](#-xternal-wrappers-argument-forwarding-)

## simulate_genotypes

**25 PASS / 0 FAIL / 0 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| n_regions = 1 returns list of length 1 | ✅ PASS | 0.72s |
| n_regions = 3 returns list of length 3 | ✅ PASS | 1.46s |
| n sets number of rows in X | ✅ PASS | 0.48s |
| p as scalar applied to all regions | ✅ PASS | 1.00s |
| p as vector sets different targets per region | ✅ PASS | 1.06s |
| p > 500 with bundled VCF warns and caps | ✅ PASS | 1.28s |
| vcf_files = NULL uses bundled example VCF | ✅ PASS | 0.67s |
| vcf_files = single path reused for all regions | ✅ PASS | 0.95s |
| vcf_files wrong length errors | ✅ PASS | 0.00s |
| vcf_files missing file errors | ✅ PASS | 0.00s |
| min_maf = 0 accepts all variants | ✅ PASS | 0.47s |
| min_maf = 0.1 (stricter filter) works | ✅ PASS | 0.49s |
| max_maf = 0.3 applies upper MAF filter | ✅ PASS | 0.45s |
| standardise = TRUE gives ~zero-mean columns | ✅ PASS | 0.45s |
| standardise = FALSE returns 0/1/2 coding | ✅ PASS | 0.44s |
| standardise = FALSE returns X_raw identical to X | ✅ PASS | 0.49s |
| genetic_map_dir = NULL (uses tempdir) works | ✅ PASS | 0.89s |
| genetic_map_dir = existing path caches maps | ✅ PASS | 0.46s |
| seed ensures reproducibility | ✅ PASS | 0.90s |
| seed = NULL accepted (no reproducibility required) | ✅ PASS | 0.49s |
| verbose = FALSE suppresses messages | ✅ PASS | 0.44s |
| verbose = TRUE prints region progress | ✅ PASS | 0.45s |
| return value has X, X_raw, n, p, maf, variant_ids, region_id, vcf_source | ✅ PASS | 0.00s |
| n_regions must be positive integer (error on 0) | ✅ PASS | 0.00s |
| p vector length mismatch errors | ✅ PASS | 0.00s |

## simulate_phenotypes

**44 PASS / 0 FAIL / 0 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| S = 1 (scalar) runs without error | ✅ PASS | 0.00s |
| S = 3 (scalar) selects 3 causal variants | ✅ PASS | 0.00s |
| S as vector (different S per region) | ✅ PASS | 0.00s |
| S vector wrong length errors | ✅ PASS | 0.00s |
| S > p errors | ✅ PASS | 0.00s |
| phi = 0.1 runs without error | ✅ PASS | 0.00s |
| phi = 0.8 runs without error | ✅ PASS | 0.00s |
| phi as vector per region | ✅ PASS | 0.00s |
| phi outside (0,1) errors | ✅ PASS | 0.00s |
| model = 'sparse' works | ✅ PASS | 0.00s |
| model = 'sparse_inf' works | ✅ PASS | 0.01s |
| model invalid string errors | ✅ PASS | 0.00s |
| p_causal = 0.2 (sparse_inf) partitions variance | ✅ PASS | 0.00s |
| p_causal = 1.0 (fully sparse, no inf component) | ✅ PASS | 0.00s |
| p_causal outside (0,1] errors | ✅ PASS | 0.00s |
| inf_model = 'beatrice' (noncausal variants only) | ✅ PASS | 0.00s |
| inf_model = 'susie_inf' (all variants) | ✅ PASS | 0.00s |
| inf_model invalid string errors | ✅ PASS | 0.00s |
| effect_distribution = 'normal' draws from N(0, effect_variance) | ✅ PASS | 0.00s |
| effect_distribution = 'equal' distributes variance equally | ✅ PASS | 0.00s |
| effect_distribution invalid errors | ✅ PASS | 0.00s |
| effect_variance = 0.1 accepted | ✅ PASS | 0.00s |
| effect_variance = 1.0 accepted | ✅ PASS | 0.00s |
| effect_variance <= 0 errors | ✅ PASS | 0.00s |
| annotations = 'none' (no annotation matrix) | ✅ PASS | 0.00s |
| annotations = 'binary' creates binary annotation matrix | ✅ PASS | 0.00s |
| annotations = 'continuous' creates continuous annotation matrix | ✅ PASS | 0.00s |
| annotations as user-supplied matrix | ✅ PASS | 0.00s |
| annotations invalid string errors | ✅ PASS | 0.00s |
| n_annotations = 1 works | ✅ PASS | 0.00s |
| n_annotations = 5 works | ✅ PASS | 0.00s |
| annotation_proportions = NULL (random proportions) | ✅ PASS | 0.00s |
| annotation_proportions scalar (same for all annotations) | ✅ PASS | 0.00s |
| annotation_proportions vector (per-annotation) | ✅ PASS | 0.00s |
| annotation_proportions vector wrong length errors | ✅ PASS | 0.00s |
| enrichment = NULL (random enrichments) | ✅ PASS | 0.00s |
| enrichment scalar (same for all annotations) | ✅ PASS | 0.00s |
| enrichment vector (per-annotation) | ✅ PASS | 0.00s |
| enrichment vector wrong length errors | ✅ PASS | 0.00s |
| seed ensures reproducibility | ✅ PASS | 0.00s |
| seed = NULL accepted | ✅ PASS | 0.00s |
| verbose = FALSE suppresses messages | ✅ PASS | 0.02s |
| return fields: y, z, beta_hat, se, LD, truth | ✅ PASS | 0.00s |
| truth fields: causal_indices, causal_effects, beta_true, pve, S, phi, model | ✅ PASS | 0.00s |

## simulate_genotypes — save / output_dir

**6 PASS / 0 FAIL / 0 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| save = FALSE writes no files | ✅ PASS | 0.51s |
| save = TRUE writes .rds file | ✅ PASS | 0.45s |
| saved .rds is readable and has correct structure | ✅ PASS | 0.46s |
| output_dir created if it does not exist | ✅ PASS | 0.45s |
| filename encodes n_regions, n, p, seed | ✅ PASS | 0.91s |
| seed = NULL gives 'noseed' tag in filename | ✅ PASS | 0.45s |

## simulate_phenotypes — save / output_dir

**5 PASS / 0 FAIL / 0 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| save = FALSE writes no files | ✅ PASS | 0.00s |
| save = TRUE writes .rds file | ✅ PASS | 0.00s |
| saved .rds has y, z, truth fields | ✅ PASS | 0.01s |
| output_dir created if it does not exist | ✅ PASS | 0.00s |
| filename encodes model, S, phi, seed | ✅ PASS | 0.00s |

## run_simulation

**31 PASS / 0 FAIL / 0 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| n_iter = 1 produces correct number of scenarios | ✅ PASS | 0.46s |
| n_iter = 3 produces correct scenario count | ✅ PASS | 0.48s |
| S vector sweeps correctly | ✅ PASS | 0.51s |
| phi vector sweeps correctly | ✅ PASS | 0.48s |
| model = 'sparse' runs without error | ✅ PASS | 0.48s |
| model = 'sparse_inf' sweeps p_causal | ✅ PASS | 0.47s |
| inf_model = 'beatrice' accepted in sparse_inf | ✅ PASS | 0.48s |
| inf_model = 'susie_inf' accepted in sparse_inf | ✅ PASS | 0.48s |
| effect_distribution = 'normal' recorded in params | ✅ PASS | 0.48s |
| effect_distribution = 'equal' works | ✅ PASS | 0.48s |
| effect_variance = 0.5 accepted | ✅ PASS | 0.48s |
| annotations = 'none' produces NULL annotation matrix | ✅ PASS | 0.47s |
| annotations = 'binary' with n_annotations = 2 | ✅ PASS | 0.49s |
| annotations = 'continuous' with n_annotations = 3 | ✅ PASS | 0.48s |
| annotation_proportions scalar passed through correctly | ✅ PASS | 0.47s |
| enrichment scalar passed through correctly | ✅ PASS | 0.51s |
| vcf_dir missing directory errors | ✅ PASS | 0.00s |
| vcf_files = single VCF path used for all regions | ✅ PASS | 0.96s |
| min_maf = 0.05 passed to simulate_genotypes | ✅ PASS | 0.49s |
| max_maf = 0.4 passed to simulate_genotypes | ✅ PASS | 0.48s |
| standardise = FALSE returns raw 0/1/2 genotypes | ✅ PASS | 0.48s |
| seed = 42 ensures reproducibility | ✅ PASS | 0.95s |
| save = TRUE writes .rds file | ✅ PASS | 0.49s |
| save = FALSE writes no files | ✅ PASS | 0.47s |
| output_dir is created if it does not exist | ✅ PASS | 0.50s |
| verbose = FALSE suppresses messages | ✅ PASS | 0.48s |
| return value has genotypes, scenarios, params | ✅ PASS | 0.00s |
| scenarios have correct fields | ✅ PASS | 0.00s |
| params records all key settings | ✅ PASS | 0.00s |
| n_iter must be a positive integer (error on 0) | ✅ PASS | 0.00s |
| phi outside (0,1) errors | ✅ PASS | 0.00s |

## run_methods

**20 PASS / 0 FAIL / 0 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| methods = 'susie' runs on SIM_MINI | ✅ PASS | 0.22s |
| methods = 'abf' runs on SIM_MINI | ✅ PASS | 0.00s |
| methods = 'susie_inf' runs on SIM_MINI | ✅ PASS | 0.27s |
| methods = 'carma' runs on SIM_MINI | ✅ PASS | 25.33s |
| multiple methods run together | ✅ PASS | 0.38s |
| method_args forwarded to susie (L and coverage) | ✅ PASS | 0.11s |
| method_args forwarded to abf (prior_variance, coverage) | ✅ PASS | 0.00s |
| method_args forwarded to susie_inf (L) | ✅ PASS | 0.16s |
| method_args forwarded to carma (rho.index) | ✅ PASS | 29.88s |
| unknown method name errors | ✅ PASS | 0.00s |
| method_args for non-run method warns | ✅ PASS | 0.00s |
| save = TRUE writes per-method .rds and run_metadata.rds | ✅ PASS | 0.01s |
| save = FALSE produces no files | ✅ PASS | 0.00s |
| verbose = FALSE suppresses messages | ✅ PASS | 0.01s |
| return value has per-method list with results, n_total, n_failed | ✅ PASS | 0.00s |
| each fit has pip, credible_sets, method, runtime_seconds | ✅ PASS | 0.00s |
| pip length equals n_snps | ✅ PASS | 0.00s |
| pip values in [0, 1] | ✅ PASS | 0.00s |
| failed fits return error field not NA pip | ✅ PASS | 0.00s |
| methods case-insensitive (SUSIE == susie) | ✅ PASS | 0.21s |

## evaluate_methods

**22 PASS / 0 FAIL / 0 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| basic evaluation returns named list per method | ✅ PASS | 0.00s |
| global stratum present for each method | ✅ PASS | 0.00s |
| by_S stratum present and named correctly | ✅ PASS | 0.00s |
| by_phi stratum present and named correctly | ✅ PASS | 0.00s |
| by_p_causal is NULL for sparse model | ✅ PASS | 0.00s |
| by_p_causal present for sparse_inf model | ✅ PASS | 0.57s |
| global auprc is numeric in [0, 1] | ✅ PASS | 0.00s |
| global cs_coverage is in [0, 1] or NA | ✅ PASS | 0.00s |
| global cs_power is in [0, 1] or NA | ✅ PASS | 0.00s |
| fdr_power_curve has required columns | ✅ PASS | 0.00s |
| pip_calibration has required columns | ✅ PASS | 0.00s |
| SE fields present when n_iter >= 2 | ✅ PASS | 0.00s |
| pip_thresholds custom (coarser) works | ✅ PASS | 0.02s |
| n_pip_cal_bins = 5 produces 5-row calibration table | ✅ PASS | 0.04s |
| n_pip_cal_bins = 20 produces 20-row calibration table | ✅ PASS | 0.04s |
| save = TRUE writes evaluation.rds and evaluation_summary.csv | ✅ PASS | 0.04s |
| save = FALSE writes no files | ✅ PASS | 0.03s |
| output_dir created if absent (save = TRUE) | ✅ PASS | 0.04s |
| verbose = FALSE suppresses messages | ✅ PASS | 0.07s |
| methods_evaluated field present in return value | ✅ PASS | 0.00s |
| simulation missing 'scenarios' errors | ✅ PASS | 0.00s |
| results missing 'methods_run' errors | ✅ PASS | 0.00s |

## plot_results

**13 PASS / 0 FAIL / 0 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| output_file explicit path writes PDF there | ✅ PASS | 1.52s |
| output_dir writes evaluation.pdf inside that directory | ✅ PASS | 1.04s |
| output_dir created if it does not exist | ✅ PASS | 1.08s |
| output_file takes precedence over output_dir | ✅ PASS | 1.10s |
| save = FALSE does not write any file | ✅ PASS | 0.00s |
| save = TRUE (default) writes PDF | ✅ PASS | 1.11s |
| methods = 'susie' only (subset of evaluated methods) | ✅ PASS | 1.09s |
| methods = c('susie', 'abf') includes both | ✅ PASS | 1.22s |
| methods subset to unknown method produces no-valid-method error | ✅ PASS | 0.00s |
| verbose = FALSE produces no messages | ✅ PASS | 1.23s |
| verbose = TRUE prints section messages | ✅ PASS | 1.32s |
| return value is the resolved output path (invisibly) | ✅ PASS | 1.37s |
| sparse_inf eval with by_p_causal section rendered | ✅ PASS | 1.90s |

## run_susie / run_susie_region — method arguments

**12 PASS / 0 FAIL / 0 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| L = 5 accepted | ✅ PASS | 0.01s |
| L = 1 (single-component) accepted | ✅ PASS | 0.00s |
| coverage = 0.5 accepted | ✅ PASS | 0.01s |
| coverage = 0.99 accepted | ✅ PASS | 0.01s |
| min_abs_corr = 0 (no purity filter) accepted | ✅ PASS | 0.01s |
| min_abs_corr = 0.8 (strict purity filter) accepted | ✅ PASS | 0.01s |
| max_iter = 50 accepted | ✅ PASS | 0.01s |
| estimate_residual_variance = FALSE accepted | ✅ PASS | 0.01s |
| estimate_prior_variance = FALSE accepted | ✅ PASS | 0.00s |
| prior_variance = 0.05 accepted | ✅ PASS | 0.01s |
| susie output has pip, credible_sets, method, runtime_seconds | ✅ PASS | 0.01s |
| susie pip sums approximately to number of credible sets | ✅ PASS | 0.01s |

## run_abf / run_abf_region — method arguments

**8 PASS / 0 FAIL / 0 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| prior_variance = 0.04 (default) accepted | ✅ PASS | 0.00s |
| prior_variance = 0.1 accepted | ✅ PASS | 0.00s |
| coverage = 0.5 accepted | ✅ PASS | 0.00s |
| coverage = 0.99 accepted | ✅ PASS | 0.00s |
| abf returns exactly one credible set | ✅ PASS | 0.00s |
| abf pip sums to 1 (normalised ABF) | ✅ PASS | 0.00s |
| abf additional contains log10_abf | ✅ PASS | 0.00s |
| larger prior_variance increases ABF magnitude | ✅ PASS | 0.00s |

## run_susie_inf / run_susie_inf_region — method arguments

**6 PASS / 0 FAIL / 0 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| L = 5 accepted | ✅ PASS | 0.02s |
| L = 1 accepted | ✅ PASS | 0.00s |
| coverage = 0.9 accepted | ✅ PASS | 0.01s |
| max_iter = 50 accepted | ✅ PASS | 0.01s |
| susie_inf output has pip, credible_sets, method | ✅ PASS | 0.01s |
| susie_inf pip in [0,1] | ✅ PASS | 0.02s |

## run_carma / run_carma_region — method arguments

**7 PASS / 0 FAIL / 0 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| rho.index = 0.95 (default) accepted | ✅ PASS | 0.18s |
| rho.index = 0.9 accepted | ✅ PASS | 0.25s |
| num.causal = 5 accepted | ✅ PASS | 0.17s |
| num.causal = 1 accepted | ✅ PASS | 0.14s |
| carma returns pip of correct length | ✅ PASS | 0.16s |
| carma pip in [0,1] | ✅ PASS | 0.16s |
| carma returns exactly one credible set (global) | ✅ PASS | 0.33s |

## External wrappers — argument forwarding (FINEMAP, PAINTOR, BEATRICE, FUNMAP)

**2 PASS / 0 FAIL / 2 SKIP**

| Test | Status | Notes |
|------|--------|-------|
| finemap: finemap_path arg recognised | ⏭️ SKIP | FINEMAP binary not available on this machine |
| paintor: paintor_path arg recognised | ⏭️ SKIP | PAINTOR binary not available / not compiled |
| beatrice: beatrice_dir and python args forwarded (graceful error if absent) | ✅ PASS | 0.02s |
| funmap: python arg forwarded (graceful error if absent) | ✅ PASS | 0.00s |

---

## Appendix — Function argument reference

### `simulate_genotypes()`

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `n_regions` | integer | 3 | Number of independent genomic regions |
| `n` | integer | 500 | Number of individuals |
| `p` | integer or vector | 200 | Target SNPs per region (scalar or per-region vector) |
| `vcf_files` | character vector or NULL | NULL | VCF files (one per region); NULL = bundled example |
| `min_maf` | numeric ∈ [0, 0.5] | 0.01 | Minimum MAF filter |
| `max_maf` | numeric or NA | NA | Maximum MAF filter; NA = no upper filter |
| `standardise` | logical | TRUE | Standardise genotypes to mean 0, variance 1 |
| `genetic_map_dir` | character or NULL | NULL | Cache directory for HapMap genetic maps |
| `seed` | integer or NULL | NULL | Random seed for reproducibility |
| `save` | logical | FALSE | Save genotype list as .rds to `output_dir` |
| `output_dir` | character | 'results' | Directory for saved output (created if absent) |
| `verbose` | logical | TRUE | Print progress messages |

### `simulate_phenotypes()`

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `genotypes` | list | — | Output from `simulate_genotypes()` |
| `S` | integer or vector | 1 | Causal variants per region (scalar or per-region vector) |
| `phi` | numeric or vector | 0.1 | PVE (scalar or per-region vector), must be in (0, 1) |
| `model` | 'sparse' or 'sparse_inf' | 'sparse' | Genetic architecture model |
| `p_causal` | numeric ∈ (0, 1] | 0.5 | Fraction of PVE from sparse component (sparse_inf only) |
| `inf_model` | 'beatrice' or 'susie_inf' | 'beatrice' | Infinitesimal component formulation (sparse_inf only) |
| `effect_distribution` | 'normal' or 'equal' | 'normal' | Effect size distribution |
| `effect_variance` | numeric > 0 | 0.36 | Variance for normal effect sizes |
| `annotations` | 'none', 'binary', 'continuous', or matrix | 'none' | Annotation mode |
| `n_annotations` | integer ≥ 1 | 3 | Number of annotation columns (for binary/continuous) |
| `annotation_proportions` | numeric, vector, or NULL | NULL | Proportion of 1s per binary annotation |
| `enrichment` | numeric, vector, or NULL | NULL | Fold-enrichment for annotation-guided selection |
| `seed` | integer or NULL | NULL | Random seed |
| `save` | logical | FALSE | Save phenotype list as .rds to `output_dir` |
| `output_dir` | character | 'results' | Directory for saved output (created if absent) |
| `verbose` | logical | TRUE | Print progress messages |

### `run_simulation()`

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `n_regions` | integer | 3 | Number of genomic regions |
| `n` | integer | 500 | Number of individuals |
| `p` | integer or vector | 200 | Target SNPs per region |
| `n_iter` | integer ≥ 1 | 5 | Replicates per parameter combination |
| `S` | integer vector | c(1,2,3,5) | Causal-variant values to sweep |
| `phi` | numeric vector ∈ (0,1) | c(0.1,0.2,0.4,0.6) | PVE values to sweep |
| `model` | 'sparse' or 'sparse_inf' | 'sparse' | Genetic architecture |
| `p_causal` | numeric vector ∈ (0,1] | c(0.1,0.2,0.4) | p_causal values to sweep (sparse_inf only) |
| `inf_model` | 'beatrice' or 'susie_inf' | 'beatrice' | Infinitesimal formulation (sparse_inf only) |
| `effect_distribution` | 'normal' or 'equal' | 'normal' | Effect size distribution |
| `effect_variance` | numeric > 0 | 0.36 | Normal effect variance |
| `annotations` | 'none', 'binary', 'continuous', or matrix | 'none' | Annotation mode |
| `n_annotations` | integer ≥ 1 | 3 | Number of annotation columns |
| `annotation_proportions` | numeric, vector, or NULL | NULL | Binary annotation proportions |
| `enrichment` | numeric, vector, or NULL | NULL | Annotation enrichment |
| `vcf_dir` | character or NULL | NULL | Directory of VCF files (from prepare_vcfs.R) |
| `vcf_files` | character vector or NULL | NULL | Explicit VCF paths (overrides vcf_dir) |
| `min_maf` | numeric | 0.01 | Minimum MAF |
| `max_maf` | numeric or NA | NA | Maximum MAF |
| `standardise` | logical | TRUE | Standardise genotypes |
| `seed` | integer or NULL | NULL | Master random seed |
| `save` | logical | FALSE | Save result as .rds |
| `output_dir` | character | 'results' | Output directory |
| `verbose` | logical | TRUE | Print progress |

### `run_methods()`

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `simulation` | list | — | Output of `run_simulation()` |
| `methods` | character vector | 'susie' | Method names to run (case-insensitive) |
| `method_args` | named list | list() | Per-method argument overrides |
| `save` | logical | FALSE | Save per-method .rds files |
| `output_dir` | character | 'results' | Output directory |
| `verbose` | logical | TRUE | Print progress |

**Supported methods and their key tuneable arguments (via `method_args`):**

| Method | Key arguments |
|--------|--------------|
| `susie` | `L`, `coverage`, `min_abs_corr`, `max_iter`, `estimate_residual_variance`, `estimate_prior_variance`, `prior_variance` |
| `susie_inf` | `L`, `coverage`, `max_iter` |
| `abf` | `prior_variance`, `coverage` |
| `carma` | `rho.index`, `num.causal` |
| `finemap` | `finemap_path`, `n_causal`, `prior_std` |
| `paintor` | `paintor_path`, `max_causal` |
| `beatrice` | `beatrice_dir`, `python`, `max_iter` |
| `funmap` | `python`, `L`, `max_iter` |

### `evaluate_methods()`

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `simulation` | list | — | Output of `run_simulation()` |
| `results` | list | — | Output of `run_methods()` |
| `pip_thresholds` | numeric vector | seq(0,1,by=0.005) | PIP thresholds for power/FDR curve |
| `n_pip_cal_bins` | integer | 10 | Equal-width bins for PIP calibration |
| `save` | logical | FALSE | Write evaluation.rds and evaluation_summary.csv |
| `output_dir` | character | 'results' | Output directory |
| `verbose` | logical | TRUE | Print progress |

### `plot_results()`

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `eval_out` | list | — | Output of `evaluate_methods()` |
| `output_file` | character or NULL | NULL | Full PDF path (overrides `output_dir` when set) |
| `output_dir` | character | 'results' | Directory to save `evaluation.pdf` when `output_file` is NULL |
| `save` | logical | TRUE | If FALSE, skip writing the PDF entirely |
| `methods` | character vector or NULL | NULL | Methods to include (NULL = all evaluated) |
| `verbose` | logical | TRUE | Print progress |

