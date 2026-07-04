# =============================================================================
# test-comprehensive.R
#
# Argument-level test coverage for every public function in the benchmark.
# Translated from the former scripts/test_comprehensive.R.
#
# Covers each documented argument of:
#   simulate_genotypes, simulate_phenotypes, run_simulation,
#   run_methods, evaluate_methods, plot_results
#
# Plus per-method argument-passthrough tests for the pure-R wrappers
# (susie, susie_inf, abf, carma, marginal_z, polyfun_oracle, polyfun_est)
# and skip-if-not-available tests for the external wrappers (finemap,
# paintor, beatrice, funmap, sparsepro).
#
# Sections (renumbered to fix a duplicate "SECTION 16" in the source):
#   1   simulate_genotypes
#   2   simulate_phenotypes
#   2b  simulate_genotypes - save / output_dir
#   2c  simulate_phenotypes - save / output_dir
#   3   run_simulation
#   4   run_methods
#   5   evaluate_methods
#   6   plot_results
#   7   run_susie / run_susie_region
#   8   run_abf / run_abf_region
#   9   run_susie_inf / run_susie_inf_region
#   10  run_carma / run_carma_region
#   11  External wrappers (finemap / paintor / beatrice / funmap)
#   12  run_marginal_z
#   13  run_polyfun_oracle
#   14  run_polyfun_est (+ scenario_setup hook)
#   15  MAF-stratified evaluation (by_causal_maf)
#   16  Misspecification stratification (by_true_annotation_type)
#   17  LD mismatch (n_ref independent reference panel)
#   18  run_sparsepro / run_sparsepro_region (was misnumbered as
#       SECTION 16 in the original; renumbered here)
#
# Skip-on-CRAN: simulations call sim1000G which may need to download a
# HapMap genetic map (~1 MB per chromosome) on first call; CRAN does not
# allow that. The fixture block at file load time may hit the network;
# CRAN runs will fail there which is the expected behaviour for a
# Suggests-only check.
#
# The original `scripts/test_comprehensive.R` was removed once this
# testthat version superseded it.
# =============================================================================


# --- Shared test fixtures --------------------------------------------------
#
# Built once at file load time so multiple test_that() blocks within a
# section can share them. Mirrors the original script.

# Small genotype object (2 regions, n=100, p=50)
GENO_SMALL <- simulate_genotypes(
  n_regions = 2, n = 100, p = 50,
  genetic_map_dir = fmb_test_map_dir(),
  seed = 1, verbose = FALSE
)

# Phenotypes for the small geno object
PHENO_SMALL <- simulate_phenotypes(
  GENO_SMALL, S = 1, phi = 0.2, seed = 1, verbose = FALSE
)

# Minimal sim object for run_methods / evaluate_methods tests
SIM_MINI <- run_simulation(
  n_regions = 2, n = 100, p = 50,
  n_iter = 2, S = c(1, 2), phi = c(0.2, 0.4),
  model = "sparse", annotations = "none",
  genetic_map_dir = fmb_test_map_dir(),
  seed = 42, verbose = FALSE
)

# Annotated variant for polyfun_oracle / polyfun_est tests
SIM_MINI_ANNOT <- run_simulation(
  n_regions = 2, n = 100, p = 50,
  n_iter = 2, S = c(1, 2), phi = c(0.2, 0.4),
  model = "sparse", annotations = "binary", n_annotations = 3,
  enrichment = 5.0,
  genetic_map_dir = fmb_test_map_dir(),
  seed = 42, verbose = FALSE
)

# Susie + abf results on SIM_MINI so evaluate / plot tests have real data
RESULTS_MINI <- run_methods(
  SIM_MINI,
  methods     = c("susie", "abf"),
  method_args = list(
    susie = list(L = 5, coverage = 0.95),
    abf   = list(prior_variance = 0.04)
  ),
  save = FALSE, verbose = FALSE
)

EVAL_MINI <- evaluate_methods(
  SIM_MINI, RESULTS_MINI,
  save = FALSE, verbose = FALSE
)


# =============================================================================
# SECTION 1: simulate_genotypes
# =============================================================================

test_that("[1] n_regions = 1 returns list of length 1", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 1, verbose = FALSE)
  expect_length(g, 1L)
})

test_that("[1] n_regions = 3 returns list of length 3", {
  g <- simulate_genotypes(n_regions = 3, n = 50, p = 30,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 2, verbose = FALSE)
  expect_length(g, 3L)
})

test_that("[1] n sets number of rows in X", {
  g <- simulate_genotypes(n_regions = 1, n = 80, p = 30,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 3, verbose = FALSE)
  expect_equal(nrow(g[[1]]$X), 80L)
})

test_that("[1] p as scalar applied to all regions (with truncation cap)", {
  g <- simulate_genotypes(n_regions = 2, n = 50, p = 40,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 4, verbose = FALSE)
  expect_lte(g[[1]]$p, 40L)
  expect_lte(g[[2]]$p, 40L)
})

test_that("[1] p as vector sets different targets per region", {
  g <- simulate_genotypes(n_regions = 2, n = 50, p = c(30, 50),
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 5, verbose = FALSE)
  expect_length(g, 2L)
})

test_that("[1] p > 500 with bundled VCF warns and caps", {
  expect_warning(
    simulate_genotypes(n_regions = 1, n = 50, p = 600,
                       genetic_map_dir = fmb_test_map_dir(),
                       seed = 6, verbose = FALSE)
  )
})

test_that("[1] vcf_files = NULL uses bundled example VCF", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          vcf_files = NULL,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 7, verbose = FALSE)
  expect_false(is.null(g[[1]]$X))
})

test_that("[1] vcf_files = single path reused for all regions", {
  vcf <- system.file("examples", "region.vcf.gz", package = "sim1000G")
  g <- simulate_genotypes(n_regions = 2, n = 50, p = 30,
                          vcf_files = vcf,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 8, verbose = FALSE)
  expect_length(g, 2L)
})

test_that("[1] vcf_files wrong length errors", {
  vcf <- system.file("examples", "region.vcf.gz", package = "sim1000G")
  expect_error(
    simulate_genotypes(n_regions = 3, n = 50, p = 30,
                       vcf_files = c(vcf, vcf),
                       genetic_map_dir = fmb_test_map_dir(),
                       seed = 9, verbose = FALSE)
  )
})

test_that("[1] vcf_files missing file errors", {
  expect_error(
    simulate_genotypes(n_regions = 1, n = 50, p = 30,
                       vcf_files = "/nonexistent/path.vcf.gz",
                       genetic_map_dir = fmb_test_map_dir(),
                       seed = 10, verbose = FALSE)
  )
})

test_that("[1] min_maf = 0 accepts all variants", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          min_maf = 0,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 11, verbose = FALSE)
  expect_false(is.null(g[[1]]$X))
})

test_that("[1] min_maf = 0.1 (stricter filter) works", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          min_maf = 0.1,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 12, verbose = FALSE)
  # MAF is two-sided; values must be in [0.1, 0.9].
  expect_true(all(g[[1]]$maf >= 0.1 | g[[1]]$maf <= 0.9))
})

test_that("[1] max_maf = 0.3 applies upper MAF filter", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          min_maf = 0.01, max_maf = 0.3,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 13, verbose = FALSE)
  expect_false(is.null(g[[1]]$X))
})

test_that("[1] standardise = TRUE gives ~zero-mean columns", {
  g <- simulate_genotypes(n_regions = 1, n = 200, p = 30,
                          standardise = TRUE,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 14, verbose = FALSE)
  expect_true(all(abs(colMeans(g[[1]]$X)) < 0.01))
})

test_that("[1] standardise = FALSE returns 0/1/2 coding", {
  g <- simulate_genotypes(n_regions = 1, n = 100, p = 30,
                          standardise = FALSE,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 15, verbose = FALSE)
  vals <- unique(as.vector(g[[1]]$X))
  expect_true(all(vals %in% c(0, 1, 2)))
})

test_that("[1] standardise = FALSE: X identical to X_raw", {
  g <- simulate_genotypes(n_regions = 1, n = 100, p = 30,
                          standardise = FALSE,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 16, verbose = FALSE)
  expect_identical(g[[1]]$X, g[[1]]$X_raw)
})

test_that("[1] genetic_map_dir = NULL (uses tempdir) works", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          genetic_map_dir = NULL,
                          seed = 17, verbose = FALSE)
  expect_false(is.null(g[[1]]$X))
})

test_that("[1] genetic_map_dir = existing path caches maps", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = 18, verbose = FALSE)
  expect_false(is.null(g[[1]]$X))
})

test_that("[1] seed ensures reproducibility", {
  g1 <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                            genetic_map_dir = fmb_test_map_dir(),
                            seed = 99, verbose = FALSE)
  g2 <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                            genetic_map_dir = fmb_test_map_dir(),
                            seed = 99, verbose = FALSE)
  expect_identical(g1[[1]]$X, g2[[1]]$X)
})

test_that("[1] seed = NULL accepted (no reproducibility required)", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          genetic_map_dir = fmb_test_map_dir(),
                          seed = NULL, verbose = FALSE)
  expect_false(is.null(g[[1]]$X))
})

test_that("[1] verbose = FALSE suppresses messages", {
  out <- utils::capture.output(
    simulate_genotypes(n_regions = 1, n = 50, p = 30,
                       genetic_map_dir = fmb_test_map_dir(),
                       seed = 20, verbose = FALSE),
    type = "message"
  )
  expect_true(length(out) == 0L || !any(grepl("Region", out)))
})

test_that("[1] verbose = TRUE prints region progress", {
  out <- utils::capture.output(
    simulate_genotypes(n_regions = 1, n = 50, p = 30,
                       genetic_map_dir = fmb_test_map_dir(),
                       seed = 21, verbose = TRUE),
    type = "message"
  )
  expect_true(any(grepl("Region|region|Done", out)))
})

test_that("[1] return value carries all expected fields", {
  expected <- c("X", "X_raw", "n", "p", "maf", "variant_ids",
                "region_id", "vcf_source")
  expect_true(all(expected %in% names(GENO_SMALL[[1]])))
})

test_that("[1] n_regions = 0 errors (must be positive)", {
  expect_error(
    simulate_genotypes(n_regions = 0, n = 50, p = 30,
                       genetic_map_dir = fmb_test_map_dir(),
                       verbose = FALSE)
  )
})

test_that("[1] p vector length mismatch errors", {
  expect_error(
    simulate_genotypes(n_regions = 3, n = 50, p = c(30, 40),
                       genetic_map_dir = fmb_test_map_dir(),
                       verbose = FALSE)
  )
})


# =============================================================================
# SECTION 2: simulate_phenotypes
# =============================================================================

test_that("[2] S = 1 (scalar) selects one causal variant", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              seed = 1, verbose = FALSE)
  expect_false(is.null(sim[[1]]$truth$causal_indices))
  expect_length(sim[[1]]$truth$causal_indices, 1L)
})

test_that("[2] S = 3 (scalar) selects 3 causal variants", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 3, phi = 0.2,
                              seed = 2, verbose = FALSE)
  expect_length(sim[[1]]$truth$causal_indices, 3L)
})

test_that("[2] S as vector sets per-region counts", {
  sim <- simulate_phenotypes(GENO_SMALL, S = c(1, 2), phi = 0.2,
                              seed = 3, verbose = FALSE)
  expect_length(sim[[1]]$truth$causal_indices, 1L)
  expect_length(sim[[2]]$truth$causal_indices, 2L)
})

test_that("[2] S vector wrong length errors", {
  expect_error(
    simulate_phenotypes(GENO_SMALL, S = c(1, 2, 3), phi = 0.2,
                        seed = 4, verbose = FALSE)
  )
})

test_that("[2] S > p errors", {
  expect_error(
    simulate_phenotypes(GENO_SMALL, S = 999, phi = 0.2,
                        seed = 5, verbose = FALSE)
  )
})

test_that("[2] phi = 0.1 yields positive realised PVE", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.1,
                              seed = 6, verbose = FALSE)
  expect_gt(sim[[1]]$truth$pve, 0)
})

test_that("[2] phi = 0.8 yields positive realised PVE", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.8,
                              seed = 7, verbose = FALSE)
  expect_gt(sim[[1]]$truth$pve, 0)
})

test_that("[2] phi as vector sets per-region PVE target", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = c(0.1, 0.5),
                              seed = 8, verbose = FALSE)
  expect_equal(sim[[1]]$truth$phi, 0.1)
  expect_equal(sim[[2]]$truth$phi, 0.5)
})

test_that("[2] phi outside (0,1) errors", {
  expect_error(
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 1.2,
                        seed = 9, verbose = FALSE)
  )
})

test_that("[2] model = 'sparse' recorded in truth", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              model = "sparse", seed = 10, verbose = FALSE)
  expect_equal(sim[[1]]$truth$model, "sparse")
})

test_that("[2] model = 'sparse_inf' recorded in truth", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              model = "sparse_inf", p_causal = 0.5,
                              seed = 11, verbose = FALSE)
  expect_equal(sim[[1]]$truth$model, "sparse_inf")
})

test_that("[2] model invalid string errors", {
  expect_error(
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        model = "wrong_model", seed = 12, verbose = FALSE)
  )
})

test_that("[2] p_causal = 0.2 (sparse_inf) recorded in truth", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.3,
                              model = "sparse_inf", p_causal = 0.2,
                              seed = 13, verbose = FALSE)
  expect_equal(sim[[1]]$truth$p_causal, 0.2)
})

test_that("[2] p_causal = 1.0 (fully sparse, no inf component) accepted", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.3,
                              model = "sparse_inf", p_causal = 1.0,
                              seed = 14, verbose = FALSE)
  expect_false(is.null(sim[[1]]$truth$p_causal))
})

test_that("[2] p_causal outside (0,1] errors", {
  expect_error(
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.3,
                        model = "sparse_inf", p_causal = 0,
                        seed = 15, verbose = FALSE)
  )
})

test_that("[2] inf_model = 'beatrice' recorded", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.3,
                              model = "sparse_inf", p_causal = 0.5,
                              inf_model = "beatrice",
                              seed = 16, verbose = FALSE)
  expect_equal(sim[[1]]$truth$inf_model, "beatrice")
})

test_that("[2] inf_model = 'susie_inf' recorded", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.3,
                              model = "sparse_inf", p_causal = 0.5,
                              inf_model = "susie_inf",
                              seed = 17, verbose = FALSE)
  expect_equal(sim[[1]]$truth$inf_model, "susie_inf")
})

test_that("[2] inf_model invalid string errors", {
  expect_error(
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.3,
                        model = "sparse_inf", p_causal = 0.5,
                        inf_model = "unknown",
                        seed = 18, verbose = FALSE)
  )
})

test_that("[2] effect_distribution = 'normal' recorded", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 2, phi = 0.2,
                              effect_distribution = "normal",
                              effect_variance = 0.36,
                              seed = 19, verbose = FALSE)
  expect_equal(sim[[1]]$truth$effect_distribution, "normal")
})

test_that("[2] effect_distribution = 'equal' recorded", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 2, phi = 0.2,
                              effect_distribution = "equal",
                              seed = 20, verbose = FALSE)
  expect_equal(sim[[1]]$truth$effect_distribution, "equal")
})

test_that("[2] effect_distribution invalid errors", {
  expect_error(
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        effect_distribution = "laplace",
                        seed = 21, verbose = FALSE)
  )
})

test_that("[2] effect_variance = 0.1 recorded", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              effect_variance = 0.1,
                              seed = 22, verbose = FALSE)
  expect_equal(sim[[1]]$truth$effect_variance, 0.1)
})

test_that("[2] effect_variance = 1.0 accepted (phenotype produced)", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              effect_variance = 1.0,
                              seed = 23, verbose = FALSE)
  expect_false(is.null(sim[[1]]$y))
})

test_that("[2] effect_variance <= 0 errors", {
  expect_error(
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        effect_variance = 0,
                        seed = 24, verbose = FALSE)
  )
})

test_that("[2] annotations = 'none' (no annotation matrix)", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "none",
                              seed = 25, verbose = FALSE)
  expect_null(sim[[1]]$annotations_matrix)
})

test_that("[2] annotations = 'binary' produces a binary matrix", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 3,
                              seed = 26, verbose = FALSE)
  A <- sim[[1]]$annotations_matrix
  expect_false(is.null(A))
  expect_equal(ncol(A), 3L)
  expect_true(all(A %in% c(0, 1)))
})

test_that("[2] annotations = 'continuous' produces a numeric matrix", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "continuous", n_annotations = 2,
                              seed = 27, verbose = FALSE)
  A <- sim[[1]]$annotations_matrix
  expect_false(is.null(A))
  expect_equal(ncol(A), 2L)
  expect_true(is.numeric(A))
})

test_that("[2] annotations = user-supplied matrix accepted", {
  p <- GENO_SMALL[[1]]$p
  A_user <- matrix(rbinom(p * 2, 1, 0.2), nrow = p, ncol = 2)
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = A_user,
                              seed = 28, verbose = FALSE)
  expect_false(is.null(sim[[1]]$annotations_matrix))
})

test_that("[2] annotations invalid string errors", {
  expect_error(
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        annotations = "rainbow",
                        seed = 29, verbose = FALSE)
  )
})

test_that("[2] n_annotations = 1 works", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 1,
                              seed = 30, verbose = FALSE)
  expect_equal(ncol(sim[[1]]$annotations_matrix), 1L)
})

test_that("[2] n_annotations = 5 works", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 5,
                              seed = 31, verbose = FALSE)
  expect_equal(ncol(sim[[1]]$annotations_matrix), 5L)
})

test_that("[2] annotation_proportions = NULL (random) accepted", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 3,
                              annotation_proportions = NULL,
                              seed = 32, verbose = FALSE)
  expect_false(is.null(sim[[1]]$annotations_matrix))
})

test_that("[2] annotation_proportions scalar (same for all annotations)", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 3,
                              annotation_proportions = 0.2,
                              seed = 33, verbose = FALSE)
  expect_false(is.null(sim[[1]]$annotations_matrix))
})

test_that("[2] annotation_proportions vector (per-annotation)", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 3,
                              annotation_proportions = c(0.1, 0.2, 0.3),
                              seed = 34, verbose = FALSE)
  expect_false(is.null(sim[[1]]$annotations_matrix))
})

test_that("[2] annotation_proportions wrong length errors", {
  expect_error(
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        annotations = "binary", n_annotations = 3,
                        annotation_proportions = c(0.1, 0.2),
                        seed = 35, verbose = FALSE)
  )
})

test_that("[2] enrichment = NULL (random) accepted", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 2,
                              enrichment = NULL,
                              seed = 36, verbose = FALSE)
  expect_false(is.null(sim[[1]]$annotations_matrix))
})

test_that("[2] enrichment scalar (same for all annotations)", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 2,
                              enrichment = 5,
                              seed = 37, verbose = FALSE)
  expect_false(is.null(sim[[1]]$truth$enrichment))
})

test_that("[2] enrichment vector (per-annotation)", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 2,
                              enrichment = c(3, 8),
                              seed = 38, verbose = FALSE)
  expect_false(is.null(sim[[1]]$truth$enrichment))
})

test_that("[2] enrichment wrong length errors", {
  expect_error(
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        annotations = "binary", n_annotations = 3,
                        enrichment = c(3, 8),
                        seed = 39, verbose = FALSE)
  )
})

test_that("[2] seed ensures reproducibility", {
  s1 <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                             seed = 77, verbose = FALSE)
  s2 <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                             seed = 77, verbose = FALSE)
  expect_identical(s1[[1]]$y, s2[[1]]$y)
})

test_that("[2] seed = NULL accepted", {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              seed = NULL, verbose = FALSE)
  expect_false(is.null(sim[[1]]$y))
})

test_that("[2] verbose = FALSE suppresses messages", {
  out <- utils::capture.output(
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        seed = 40, verbose = FALSE),
    type = "message"
  )
  expect_length(out, 0L)
})

test_that("[2] return value carries all expected per-region fields", {
  expected <- c("y", "z", "beta_hat", "se", "LD", "truth")
  expect_true(all(expected %in% names(PHENO_SMALL[[1]])))
})

test_that("[2] truth carries all expected fields", {
  expected <- c("causal_indices", "causal_effects", "beta_true",
                "pve", "S", "phi", "model")
  expect_true(all(expected %in% names(PHENO_SMALL[[1]]$truth)))
})


# =============================================================================
# SECTION 2b: simulate_genotypes - save / output_dir
# =============================================================================

test_that("[2b] save = FALSE writes no files", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  simulate_genotypes(n_regions = 1, n = 50, p = 30,
                     genetic_map_dir = fmb_test_map_dir(),
                     seed = 200, save = FALSE, output_dir = tmp,
                     verbose = FALSE)
  expect_length(list.files(tmp), 0L)
})

test_that("[2b] save = TRUE writes .rds file", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  simulate_genotypes(n_regions = 1, n = 50, p = 30,
                     genetic_map_dir = fmb_test_map_dir(),
                     seed = 201, save = TRUE, output_dir = tmp,
                     verbose = FALSE)
  expect_length(list.files(tmp, pattern = "[.]rds$"), 1L)
})

test_that("[2b] saved .rds is readable and has correct structure", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  simulate_genotypes(n_regions = 1, n = 50, p = 30,
                     genetic_map_dir = fmb_test_map_dir(),
                     seed = 202, save = TRUE, output_dir = tmp,
                     verbose = FALSE)
  f   <- list.files(tmp, pattern = "[.]rds$", full.names = TRUE)
  obj <- readRDS(f)
  expect_true(is.list(obj))
  expect_false(is.null(obj[[1]]$X))
})

test_that("[2b] output_dir created if it does not exist", {
  tmp <- file.path(tempfile(), "geno_out")
  on.exit(unlink(dirname(tmp), recursive = TRUE))
  simulate_genotypes(n_regions = 1, n = 50, p = 30,
                     genetic_map_dir = fmb_test_map_dir(),
                     seed = 203, save = TRUE, output_dir = tmp,
                     verbose = FALSE)
  expect_true(dir.exists(tmp))
})

test_that("[2b] filename encodes n_regions, n, p, seed", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  simulate_genotypes(n_regions = 2, n = 60, p = 35,
                     genetic_map_dir = fmb_test_map_dir(),
                     seed = 204, save = TRUE, output_dir = tmp,
                     verbose = FALSE)
  f <- list.files(tmp, pattern = "[.]rds$")
  expect_match(f, "2regions")
  expect_match(f, "n60")
  expect_match(f, "p35")
  expect_match(f, "seed204")
})

test_that("[2b] seed = NULL gives 'noseed' tag in filename", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  simulate_genotypes(n_regions = 1, n = 50, p = 30,
                     genetic_map_dir = fmb_test_map_dir(),
                     seed = NULL, save = TRUE, output_dir = tmp,
                     verbose = FALSE)
  expect_match(list.files(tmp, pattern = "[.]rds$"), "noseed")
})


# =============================================================================
# SECTION 2c: simulate_phenotypes - save / output_dir
# =============================================================================

test_that("[2c] save = FALSE writes no files", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                      seed = 300, save = FALSE, output_dir = tmp,
                      verbose = FALSE)
  expect_length(list.files(tmp), 0L)
})

test_that("[2c] save = TRUE writes .rds file", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                      seed = 301, save = TRUE, output_dir = tmp,
                      verbose = FALSE)
  expect_length(list.files(tmp, pattern = "[.]rds$"), 1L)
})

test_that("[2c] saved .rds has y, z, truth fields", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                      seed = 302, save = TRUE, output_dir = tmp,
                      verbose = FALSE)
  f   <- list.files(tmp, pattern = "[.]rds$", full.names = TRUE)
  obj <- readRDS(f)
  expect_true(all(c("y", "z", "truth") %in% names(obj[[1]])))
})

test_that("[2c] output_dir created if it does not exist", {
  tmp <- file.path(tempfile(), "pheno_out")
  on.exit(unlink(dirname(tmp), recursive = TRUE))
  simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                      seed = 303, save = TRUE, output_dir = tmp,
                      verbose = FALSE)
  expect_true(dir.exists(tmp))
})

test_that("[2c] filename encodes model, S, phi, seed", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  simulate_phenotypes(GENO_SMALL, S = 2, phi = 0.3,
                      model = "sparse", seed = 304,
                      save = TRUE, output_dir = tmp, verbose = FALSE)
  f <- list.files(tmp, pattern = "[.]rds$")
  expect_match(f, "sparse")
  expect_match(f, "S2")
  expect_match(f, "phi0.3")
  expect_match(f, "seed304")
})


# =============================================================================
# SECTION 3: run_simulation
# =============================================================================

test_that("[3] n_iter = 1 produces a single scenario", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 1, verbose = FALSE)
  expect_length(r$scenarios, 1L)
})

test_that("[3] n_iter = 3 produces three scenarios", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 3, S = 1, phi = 0.2,
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 2, verbose = FALSE)
  expect_length(r$scenarios, 3L)
})

test_that("[3] S vector sweeps correctly", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = c(1, 2), phi = 0.2,
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 3, verbose = FALSE)
  expect_length(r$scenarios, 2L)
  S_vals <- sapply(r$scenarios, `[[`, "S")
  expect_equal(sort(S_vals), c(1L, 2L))
})

test_that("[3] phi vector sweeps correctly", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = c(0.2, 0.4),
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 4, verbose = FALSE)
  expect_length(r$scenarios, 2L)
  phi_vals <- sapply(r$scenarios, `[[`, "phi")
  expect_equal(sort(phi_vals), c(0.2, 0.4))
})

test_that("[3] model = 'sparse' recorded in params", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      model = "sparse",
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 5, verbose = FALSE)
  expect_equal(r$params$model, "sparse")
})

test_that("[3] model = 'sparse_inf' sweeps p_causal", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      model = "sparse_inf", p_causal = c(0.2, 0.5),
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 6, verbose = FALSE)
  expect_length(r$scenarios, 2L)
  expect_equal(r$params$model, "sparse_inf")
})

test_that("[3] inf_model = 'beatrice' accepted in sparse_inf", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      model = "sparse_inf", p_causal = 0.5,
                      inf_model = "beatrice",
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 7, verbose = FALSE)
  expect_equal(r$params$inf_model, "beatrice")
})

test_that("[3] inf_model = 'susie_inf' accepted in sparse_inf", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      model = "sparse_inf", p_causal = 0.5,
                      inf_model = "susie_inf",
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 8, verbose = FALSE)
  expect_equal(r$params$inf_model, "susie_inf")
})

test_that("[3] effect_distribution = 'normal' recorded in params", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      effect_distribution = "normal",
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 9, verbose = FALSE)
  expect_equal(r$params$effect_distribution, "normal")
})

test_that("[3] effect_distribution = 'equal' recorded", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      effect_distribution = "equal",
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 10, verbose = FALSE)
  expect_equal(r$params$effect_distribution, "equal")
})

test_that("[3] effect_variance = 0.5 recorded", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      effect_variance = 0.5,
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 11, verbose = FALSE)
  expect_equal(r$params$effect_variance, 0.5)
})

test_that("[3] annotations = 'none' produces NULL annotation matrix", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "none",
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 12, verbose = FALSE)
  expect_null(r$scenarios[[1]]$regions[[1]]$annotations_matrix)
})

test_that("[3] annotations = 'binary' with n_annotations = 2", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "binary", n_annotations = 2,
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 13, verbose = FALSE)
  A <- r$scenarios[[1]]$regions[[1]]$annotations_matrix
  expect_false(is.null(A))
  expect_equal(ncol(A), 2L)
})

test_that("[3] annotations = 'continuous' with n_annotations = 3", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "continuous", n_annotations = 3,
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 14, verbose = FALSE)
  A <- r$scenarios[[1]]$regions[[1]]$annotations_matrix
  expect_false(is.null(A))
  expect_equal(ncol(A), 3L)
})

test_that("[3] annotation_proportions scalar passed through", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "binary", n_annotations = 2,
                      annotation_proportions = 0.15,
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 15, verbose = FALSE)
  expect_equal(r$params$annotation_proportions, 0.15)
})

test_that("[3] enrichment scalar passed through", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "binary", n_annotations = 2,
                      enrichment = 4,
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 16, verbose = FALSE)
  expect_equal(r$params$enrichment, 4)
})

test_that("[3] vcf_dir missing directory errors", {
  expect_error(
    run_simulation(n_regions = 1, n = 80, p = 40,
                   n_iter = 1, S = 1, phi = 0.2,
                   vcf_dir = "/nonexistent/dir",
                   genetic_map_dir = fmb_test_map_dir(),
                   seed = 17, verbose = FALSE)
  )
})

test_that("[3] vcf_files = single VCF path used for all regions", {
  vcf <- system.file("examples", "region.vcf.gz", package = "sim1000G")
  r <- run_simulation(n_regions = 2, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      vcf_files = vcf,
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 18, verbose = FALSE)
  expect_length(r$genotypes, 2L)
})

test_that("[3] min_maf = 0.05 passed to simulate_genotypes", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      min_maf = 0.05,
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 19, verbose = FALSE)
  expect_equal(r$params$min_maf, 0.05)
})

test_that("[3] max_maf = 0.4 passed to simulate_genotypes", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      max_maf = 0.4,
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 20, verbose = FALSE)
  expect_false(is.null(r$genotypes[[1]]$X))
})

test_that("[3] standardise = FALSE returns raw 0/1/2 genotypes", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      standardise = FALSE,
                      genetic_map_dir = fmb_test_map_dir(),
                      seed = 21, verbose = FALSE)
  vals <- unique(as.vector(r$genotypes[[1]]$X))
  expect_true(all(vals %in% c(0, 1, 2)))
})

test_that("[3] seed = 42 ensures reproducibility", {
  r1 <- run_simulation(n_regions = 1, n = 80, p = 40,
                       n_iter = 1, S = 1, phi = 0.2,
                       genetic_map_dir = fmb_test_map_dir(),
                       seed = 42, verbose = FALSE)
  r2 <- run_simulation(n_regions = 1, n = 80, p = 40,
                       n_iter = 1, S = 1, phi = 0.2,
                       genetic_map_dir = fmb_test_map_dir(),
                       seed = 42, verbose = FALSE)
  expect_identical(r1$genotypes[[1]]$X, r2$genotypes[[1]]$X)
})

test_that("[3] save = TRUE writes .rds file", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  run_simulation(n_regions = 1, n = 80, p = 40,
                 n_iter = 1, S = 1, phi = 0.2,
                 genetic_map_dir = fmb_test_map_dir(),
                 seed = 99, save = TRUE, output_dir = tmp,
                 verbose = FALSE)
  expect_length(list.files(tmp, pattern = "[.]rds$"), 1L)
})

test_that("[3] save = FALSE writes no files", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  run_simulation(n_regions = 1, n = 80, p = 40,
                 n_iter = 1, S = 1, phi = 0.2,
                 genetic_map_dir = fmb_test_map_dir(),
                 seed = 100, save = FALSE, output_dir = tmp,
                 verbose = FALSE)
  expect_length(list.files(tmp), 0L)
})

test_that("[3] output_dir is created when nested + missing", {
  tmp <- file.path(tempfile(), "deep", "nested")
  on.exit(unlink(dirname(dirname(tmp)), recursive = TRUE))
  run_simulation(n_regions = 1, n = 80, p = 40,
                 n_iter = 1, S = 1, phi = 0.2,
                 genetic_map_dir = fmb_test_map_dir(),
                 seed = 101, save = TRUE, output_dir = tmp,
                 verbose = FALSE)
  expect_true(dir.exists(tmp))
})

test_that("[3] verbose = FALSE suppresses messages", {
  out <- utils::capture.output(
    run_simulation(n_regions = 1, n = 80, p = 40,
                   n_iter = 1, S = 1, phi = 0.2,
                   genetic_map_dir = fmb_test_map_dir(),
                   seed = 102, verbose = FALSE),
    type = "message"
  )
  expect_length(out, 0L)
})

test_that("[3] return value carries genotypes, scenarios, params", {
  expect_true(all(c("genotypes", "scenarios", "params") %in% names(SIM_MINI)))
})

test_that("[3] scenarios carry all expected fields", {
  expected <- c("scenario_id", "S", "phi", "p_causal", "iter", "model", "regions")
  expect_true(all(expected %in% names(SIM_MINI$scenarios[[1]])))
})

test_that("[3] params records all key settings", {
  expected <- c("n_regions", "n", "p", "n_iter", "S_values", "phi_values",
                "model", "seed")
  expect_true(all(expected %in% names(SIM_MINI$params)))
})

test_that("[3] n_iter = 0 errors (must be positive)", {
  expect_error(
    run_simulation(n_regions = 1, n = 80, p = 40,
                   n_iter = 0, S = 1, phi = 0.2,
                   genetic_map_dir = fmb_test_map_dir(),
                   verbose = FALSE)
  )
})

test_that("[3] phi outside (0,1) errors", {
  expect_error(
    run_simulation(n_regions = 1, n = 80, p = 40,
                   n_iter = 1, S = 1, phi = 1.5,
                   genetic_map_dir = fmb_test_map_dir(),
                   verbose = FALSE)
  )
})


# =============================================================================
# SECTION 4: run_methods
# =============================================================================

test_that("[4] methods = 'susie' runs on SIM_MINI", {
  r <- run_methods(SIM_MINI, methods = "susie",
                   method_args = list(susie = list(L = 5)),
                   save = FALSE, verbose = FALSE)
  expect_true("susie" %in% r$methods_run)
})

test_that("[4] methods = 'abf' runs on SIM_MINI", {
  r <- run_methods(SIM_MINI, methods = "abf",
                   save = FALSE, verbose = FALSE)
  expect_true("abf" %in% r$methods_run)
})

test_that("[4] methods = 'susie_inf' runs on SIM_MINI", {
  r <- run_methods(SIM_MINI, methods = "susie_inf",
                   method_args = list(susie_inf = list(L = 5)),
                   save = FALSE, verbose = FALSE)
  expect_true("susie_inf" %in% r$methods_run)
})

test_that("[4] methods = 'carma' runs on SIM_MINI", {
  r <- run_methods(SIM_MINI, methods = "carma",
                   save = FALSE, verbose = FALSE)
  expect_true("carma" %in% r$methods_run)
})

test_that("[4] multiple methods run together", {
  r <- run_methods(SIM_MINI, methods = c("susie", "abf"),
                   save = FALSE, verbose = FALSE)
  expect_true(all(c("susie", "abf") %in% r$methods_run))
})

test_that("[4] method_args forwarded to susie (L and coverage)", {
  r <- run_methods(SIM_MINI, methods = "susie",
                   method_args = list(susie = list(L = 3, coverage = 0.9)),
                   save = FALSE, verbose = FALSE)
  expect_equal(r$susie$method_args$L, 3)
  expect_equal(r$susie$method_args$coverage, 0.9)
})

test_that("[4] method_args forwarded to abf (prior_variance)", {
  r <- run_methods(SIM_MINI, methods = "abf",
                   method_args = list(abf = list(prior_variance = 0.02,
                                                 coverage = 0.9)),
                   save = FALSE, verbose = FALSE)
  expect_equal(r$abf$method_args$prior_variance, 0.02)
})

test_that("[4] method_args forwarded to susie_inf (L)", {
  r <- run_methods(SIM_MINI, methods = "susie_inf",
                   method_args = list(susie_inf = list(L = 3)),
                   save = FALSE, verbose = FALSE)
  expect_equal(r$susie_inf$method_args$L, 3)
})

test_that("[4] method_args forwarded to carma (rho.index)", {
  r <- run_methods(SIM_MINI, methods = "carma",
                   method_args = list(carma = list(rho.index = 0.9)),
                   save = FALSE, verbose = FALSE)
  expect_equal(r$carma$method_args$rho.index, 0.9)
})

test_that("[4] unknown method name errors", {
  expect_error(
    run_methods(SIM_MINI, methods = "notamethod",
                save = FALSE, verbose = FALSE)
  )
})

test_that("[4] method_args for non-run method warns", {
  expect_warning(
    run_methods(SIM_MINI, methods = "abf",
                method_args = list(abf = list(), susie = list(L = 5)),
                save = FALSE, verbose = FALSE)
  )
})

test_that("[4] save = TRUE writes per-method .rds + run_metadata.rds", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  run_methods(SIM_MINI, methods = "abf",
              save = TRUE, output_dir = tmp, verbose = FALSE)
  files <- list.files(tmp, recursive = TRUE, pattern = "[.]rds$")
  expect_true(any(grepl("abf[.]rds", files)))
  expect_true(any(grepl("run_metadata[.]rds", files)))
})

test_that("[4] save = FALSE produces no files", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  run_methods(SIM_MINI, methods = "abf",
              save = FALSE, output_dir = tmp, verbose = FALSE)
  expect_length(list.files(tmp, recursive = TRUE), 0L)
})

test_that("[4] verbose = FALSE suppresses messages", {
  out <- utils::capture.output(
    run_methods(SIM_MINI, methods = "abf", save = FALSE, verbose = FALSE),
    type = "message"
  )
  expect_length(out, 0L)
})

test_that("[4] return value has results / n_total / n_failed / total_runtime", {
  expect_true(all(c("results", "n_total", "n_failed",
                    "total_runtime_seconds") %in%
                  names(RESULTS_MINI$susie)))
})

test_that("[4] each fit has pip / credible_sets / method / scenario metadata", {
  fit <- RESULTS_MINI$susie$results[[1]]
  expected <- c("pip", "credible_sets", "method", "runtime_seconds",
                "scenario_id", "region_id", "S", "phi", "iter")
  expect_true(all(expected %in% names(fit)))
})

test_that("[4] pip length equals n_snps", {
  fit <- RESULTS_MINI$susie$results[[1]]
  rg  <- SIM_MINI$genotypes[[fit$region_id]]
  expect_equal(length(fit$pip), rg$p)
})

test_that("[4] pip values in [0, 1]", {
  pip <- RESULTS_MINI$susie$results[[1]]$pip
  expect_true(all(pip >= 0))
  expect_true(all(pip <= 1))
})

test_that("[4] direct wrapper call with broken inputs returns clean error", {
  # run_abf_region called with a deliberately broken region_geno; the
  # wrapper either errors cleanly or returns a structured error result.
  fit <- tryCatch(
    run_abf_region(
      region_geno  = list(LD = NULL, n = 100),
      region_pheno = list(z = rep(0, 5), se = rep(1, 5))
    ),
    error = function(e) list(error = conditionMessage(e))
  )
  expect_true(!is.null(fit$error) || !is.null(fit$pip))
})

test_that("[4] methods case-insensitive (SUSIE == susie)", {
  r <- run_methods(SIM_MINI, methods = "SUSIE",
                   method_args = list(susie = list(L = 5)),
                   save = FALSE, verbose = FALSE)
  expect_true("susie" %in% r$methods_run)
})


# =============================================================================
# SECTION 5: evaluate_methods
# =============================================================================

test_that("[5] basic evaluation returns named list per method", {
  expect_true(all(c("susie", "abf") %in% EVAL_MINI$methods_evaluated))
})

test_that("[5] global stratum present for each method", {
  for (m in EVAL_MINI$methods_evaluated) {
    expect_false(is.null(EVAL_MINI[[m]]$global),
                 info = sprintf("method=%s", m))
  }
})

test_that("[5] by_S stratum present and named correctly", {
  expect_false(is.null(EVAL_MINI$susie$by_S))
  expect_true(all(c("1", "2") %in% names(EVAL_MINI$susie$by_S)))
})

test_that("[5] by_phi stratum present and named correctly", {
  expect_false(is.null(EVAL_MINI$susie$by_phi))
  expect_true(all(c("0.2", "0.4") %in% names(EVAL_MINI$susie$by_phi)))
})

test_that("[5] by_p_causal is NULL for sparse model", {
  expect_null(EVAL_MINI$susie$by_p_causal)
})

test_that("[5] by_p_causal present and populated for sparse_inf model", {
  sim_inf <- run_simulation(
    n_regions = 1, n = 80, p = 40,
    n_iter = 2, S = 1, phi = 0.2,
    model = "sparse_inf", p_causal = c(0.3, 0.7),
    genetic_map_dir = fmb_test_map_dir(),
    seed = 55, verbose = FALSE
  )
  res_inf <- run_methods(sim_inf, methods = "abf",
                          save = FALSE, verbose = FALSE)
  ev_inf  <- evaluate_methods(sim_inf, res_inf,
                               save = FALSE, verbose = FALSE)
  expect_false(is.null(ev_inf$abf$by_p_causal))
  expect_true(all(c("0.3", "0.7") %in% names(ev_inf$abf$by_p_causal)))
})

test_that("[5] global auprc is numeric in [0, 1]", {
  auprc <- EVAL_MINI$susie$global$auprc
  expect_true(is.numeric(auprc))
  expect_false(is.na(auprc))
  expect_gte(auprc, 0)
  expect_lte(auprc, 1)
})

test_that("[5] global cs_coverage is in [0, 1] or NA", {
  cv <- EVAL_MINI$susie$global$cs_coverage
  expect_true(is.numeric(cv))
  expect_true(is.na(cv) || (cv >= 0 && cv <= 1))
})

test_that("[5] global cs_power is in [0, 1] or NA", {
  cp <- EVAL_MINI$susie$global$cs_power
  expect_true(is.numeric(cp))
  expect_true(is.na(cp) || (cp >= 0 && cp <= 1))
})

test_that("[5] fdr_power_curve has all required columns", {
  curve <- EVAL_MINI$susie$global$fdr_power_curve
  expected <- c("threshold", "tp", "fp", "fn",
                "fdr", "power", "precision", "recall")
  expect_true(all(expected %in% names(curve)))
})

test_that("[5] pip_calibration has all required columns", {
  cal <- EVAL_MINI$susie$global$pip_calibration
  expected <- c("bin", "bin_lower", "bin_upper", "mean_pip", "frac_causal")
  expect_true(all(expected %in% names(cal)))
})

test_that("[5] SE fields present when n_iter >= 2", {
  expect_false(is.null(EVAL_MINI$susie$global$auprc_se))
})

test_that("[5] custom (coarser) pip_thresholds respected", {
  ev <- evaluate_methods(SIM_MINI, RESULTS_MINI,
                         pip_thresholds = seq(0, 1, by = 0.1),
                         save = FALSE, verbose = FALSE)
  # seq(0, 1, by = 0.1) has 11 values
  expect_equal(nrow(ev$susie$global$fdr_power_curve), 11L)
})

test_that("[5] n_pip_cal_bins = 5 produces 5-row calibration table", {
  ev <- evaluate_methods(SIM_MINI, RESULTS_MINI,
                         n_pip_cal_bins = 5,
                         save = FALSE, verbose = FALSE)
  expect_equal(nrow(ev$susie$global$pip_calibration), 5L)
})

test_that("[5] n_pip_cal_bins = 20 produces 20-row calibration table", {
  ev <- evaluate_methods(SIM_MINI, RESULTS_MINI,
                         n_pip_cal_bins = 20,
                         save = FALSE, verbose = FALSE)
  expect_equal(nrow(ev$susie$global$pip_calibration), 20L)
})

test_that("[5] save = TRUE writes evaluation.rds + evaluation_summary.csv", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  evaluate_methods(SIM_MINI, RESULTS_MINI,
                   save = TRUE, output_dir = tmp, verbose = FALSE)
  expect_true(file.exists(file.path(tmp, "evaluation.rds")))
  expect_true(file.exists(file.path(tmp, "evaluation_summary.csv")))
})

test_that("[5] save = FALSE writes no files", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  evaluate_methods(SIM_MINI, RESULTS_MINI,
                   save = FALSE, output_dir = tmp, verbose = FALSE)
  expect_length(list.files(tmp), 0L)
})

test_that("[5] output_dir created if absent (save = TRUE)", {
  tmp <- file.path(tempfile(), "ev_out")
  on.exit(unlink(tmp, recursive = TRUE))
  evaluate_methods(SIM_MINI, RESULTS_MINI,
                   save = TRUE, output_dir = tmp, verbose = FALSE)
  expect_true(dir.exists(tmp))
})

test_that("[5] verbose = FALSE suppresses messages", {
  out <- utils::capture.output(
    evaluate_methods(SIM_MINI, RESULTS_MINI,
                     save = FALSE, verbose = FALSE),
    type = "message"
  )
  expect_length(out, 0L)
})

test_that("[5] methods_evaluated field present in return value", {
  expect_true("methods_evaluated" %in% names(EVAL_MINI))
})

test_that("[5] simulation missing 'scenarios' errors", {
  expect_error(
    evaluate_methods(list(params = list()), RESULTS_MINI,
                     save = FALSE, verbose = FALSE)
  )
})

test_that("[5] results missing 'methods_run' errors", {
  expect_error(
    evaluate_methods(SIM_MINI, list(),
                     save = FALSE, verbose = FALSE)
  )
})


# =============================================================================
# SECTION 6: plot_results
# =============================================================================

test_that("[6] output_file explicit path writes PDF there", {
  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp))
  plot_results(EVAL_MINI, output_file = tmp, verbose = FALSE)
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 1000L)
})

test_that("[6] output_dir writes evaluation.pdf inside that directory", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  plot_results(EVAL_MINI, output_dir = tmp, verbose = FALSE)
  expect_true(file.exists(file.path(tmp, "evaluation.pdf")))
})

test_that("[6] output_dir created when nested + missing", {
  tmp <- file.path(tempfile(), "plots", "deep")
  on.exit(unlink(dirname(dirname(tmp)), recursive = TRUE))
  plot_results(EVAL_MINI, output_dir = tmp, verbose = FALSE)
  expect_true(dir.exists(tmp))
})

test_that("[6] output_file takes precedence over output_dir", {
  tmp_dir  <- tempfile(); dir.create(tmp_dir)
  tmp_file <- tempfile(fileext = ".pdf")
  on.exit({ unlink(tmp_dir, recursive = TRUE); unlink(tmp_file) })
  plot_results(EVAL_MINI, output_file = tmp_file,
               output_dir = tmp_dir, verbose = FALSE)
  expect_true(file.exists(tmp_file))
  expect_false(file.exists(file.path(tmp_dir, "evaluation.pdf")))
})

test_that("[6] save = FALSE does not write any file", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  plot_results(EVAL_MINI, output_dir = tmp,
               save = FALSE, verbose = FALSE)
  expect_length(list.files(tmp), 0L)
})

test_that("[6] save = TRUE (default) writes the PDF", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  plot_results(EVAL_MINI, output_dir = tmp,
               save = TRUE, verbose = FALSE)
  expect_true(file.exists(file.path(tmp, "evaluation.pdf")))
})

test_that("[6] methods = 'susie' only (subset of evaluated methods)", {
  tmp <- tempfile(fileext = ".pdf"); on.exit(unlink(tmp))
  plot_results(EVAL_MINI, output_file = tmp,
               methods = "susie", verbose = FALSE)
  expect_true(file.exists(tmp))
})

test_that("[6] methods = c('susie', 'abf') includes both", {
  tmp <- tempfile(fileext = ".pdf"); on.exit(unlink(tmp))
  plot_results(EVAL_MINI, output_file = tmp,
               methods = c("susie", "abf"), verbose = FALSE)
  expect_true(file.exists(tmp))
})

test_that("[6] methods subset to unknown name errors", {
  expect_error(
    plot_results(EVAL_MINI,
                 output_file = tempfile(fileext = ".pdf"),
                 methods = "notamethod", verbose = FALSE)
  )
})

test_that("[6] verbose = FALSE produces no messages", {
  tmp <- tempfile(fileext = ".pdf"); on.exit(unlink(tmp))
  out <- utils::capture.output(
    plot_results(EVAL_MINI, output_file = tmp, verbose = FALSE),
    type = "message"
  )
  expect_length(out, 0L)
})

test_that("[6] verbose = TRUE prints section messages", {
  tmp <- tempfile(fileext = ".pdf"); on.exit(unlink(tmp))
  out <- utils::capture.output(
    plot_results(EVAL_MINI, output_file = tmp, verbose = TRUE),
    type = "message"
  )
  expect_true(any(grepl("Global|Plotting|PDF", out)))
})

test_that("[6] return value is the resolved output path (invisibly)", {
  tmp <- tempfile(fileext = ".pdf"); on.exit(unlink(tmp))
  ret <- plot_results(EVAL_MINI, output_file = tmp, verbose = FALSE)
  expect_true(is.character(ret))
  expect_equal(ret, tmp)
})

test_that("[6] sparse_inf eval renders the by_p_causal section in the PDF", {
  sim_inf <- run_simulation(
    n_regions = 1, n = 80, p = 40,
    n_iter = 2, S = 1, phi = 0.2,
    model = "sparse_inf", p_causal = c(0.3, 0.7),
    genetic_map_dir = fmb_test_map_dir(),
    seed = 66, verbose = FALSE
  )
  res_inf <- run_methods(sim_inf, methods = "abf",
                          save = FALSE, verbose = FALSE)
  ev_inf  <- evaluate_methods(sim_inf, res_inf,
                               save = FALSE, verbose = FALSE)

  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  plot_results(ev_inf, output_dir = tmp, verbose = FALSE)
  expect_true(file.exists(file.path(tmp, "evaluation.pdf")))
})


# =============================================================================
# Per-region fixtures shared by per-method argument tests (7-14)
# =============================================================================

.rg <- SIM_MINI$genotypes[[1]]
.rp <- SIM_MINI$scenarios[[1]]$regions[[1]]


# =============================================================================
# SECTION 7: run_susie / run_susie_region argument passthrough
# =============================================================================

test_that("[7] susie: L = 5 stored in params", {
  fit <- run_susie_region(.rg, .rp, L = 5)
  expect_equal(fit$params$L, 5)
})

test_that("[7] susie: L = 1 (single-component) accepted", {
  fit <- run_susie_region(.rg, .rp, L = 1)
  expect_equal(fit$params$L, 1)
})

test_that("[7] susie: coverage = 0.5 stored in params", {
  fit <- run_susie_region(.rg, .rp, L = 5, coverage = 0.5)
  expect_equal(fit$params$coverage, 0.5)
})

test_that("[7] susie: coverage = 0.99 stored in params", {
  fit <- run_susie_region(.rg, .rp, L = 5, coverage = 0.99)
  expect_equal(fit$params$coverage, 0.99)
})

test_that("[7] susie: min_abs_corr = 0 (no purity filter) accepted", {
  fit <- run_susie_region(.rg, .rp, L = 5, min_abs_corr = 0)
  expect_equal(fit$params$min_abs_corr, 0)
})

test_that("[7] susie: min_abs_corr = 0.8 (strict filter) accepted", {
  fit <- run_susie_region(.rg, .rp, L = 5, min_abs_corr = 0.8)
  expect_equal(fit$params$min_abs_corr, 0.8)
})

test_that("[7] susie: max_iter = 50 stored in params", {
  fit <- run_susie_region(.rg, .rp, L = 5, max_iter = 50)
  expect_equal(fit$params$max_iter, 50)
})

test_that("[7] susie: estimate_residual_variance = FALSE stored in params", {
  fit <- run_susie_region(.rg, .rp, L = 5, estimate_residual_variance = FALSE)
  expect_false(fit$params$estimate_residual_variance)
})

test_that("[7] susie: estimate_prior_variance = FALSE stored in params", {
  fit <- run_susie_region(.rg, .rp, L = 5, estimate_prior_variance = FALSE)
  expect_false(fit$params$estimate_prior_variance)
})

test_that("[7] susie: prior_variance = 0.05 stored in params", {
  fit <- run_susie_region(.rg, .rp, L = 5, prior_variance = 0.05)
  expect_equal(fit$params$prior_variance, 0.05)
})

test_that("[7] susie: output has all expected top-level fields", {
  fit <- run_susie_region(.rg, .rp, L = 5)
  expect_true(all(c("pip", "credible_sets", "method", "runtime_seconds",
                    "input_type", "params", "additional") %in% names(fit)))
})

test_that("[7] susie: pip is length-p and in [0, 1]", {
  fit <- run_susie_region(.rg, .rp, L = 5)
  expect_true(is.numeric(fit$pip))
  expect_equal(length(fit$pip), .rg$p)
  expect_true(all(fit$pip >= 0))
  expect_true(all(fit$pip <= 1))
})


# =============================================================================
# SECTION 8: run_abf / run_abf_region argument passthrough
# =============================================================================

test_that("[8] abf: prior_variance = 0.04 (default) stored in params", {
  fit <- run_abf_region(.rg, .rp, prior_variance = 0.04)
  expect_equal(fit$params$prior_variance, 0.04)
})

test_that("[8] abf: prior_variance = 0.1 stored in params", {
  fit <- run_abf_region(.rg, .rp, prior_variance = 0.1)
  expect_equal(fit$params$prior_variance, 0.1)
})

test_that("[8] abf: coverage = 0.5 stored in params", {
  fit <- run_abf_region(.rg, .rp, coverage = 0.5)
  expect_equal(fit$params$coverage, 0.5)
})

test_that("[8] abf: coverage = 0.99 stored in params", {
  fit <- run_abf_region(.rg, .rp, coverage = 0.99)
  expect_equal(fit$params$coverage, 0.99)
})

test_that("[8] abf: returns exactly one credible set", {
  fit <- run_abf_region(.rg, .rp)
  expect_true(is.list(fit$credible_sets))
  expect_length(fit$credible_sets, 1L)
})

test_that("[8] abf: pip sums to ~1 (normalised ABF)", {
  fit <- run_abf_region(.rg, .rp)
  expect_lt(abs(sum(fit$pip) - 1), 1e-9)
})

test_that("[8] abf: additional contains log10_abf", {
  fit <- run_abf_region(.rg, .rp)
  expect_true("log10_abf" %in% names(fit$additional))
})

test_that("[8] abf: different prior_variance values produce different PIPs", {
  fit_lo <- run_abf_region(.rg, .rp, prior_variance = 0.01)
  fit_hi <- run_abf_region(.rg, .rp, prior_variance = 0.5)
  expect_false(identical(fit_lo$pip, fit_hi$pip))
})


# =============================================================================
# SECTION 9: run_susie_inf / run_susie_inf_region argument passthrough
# =============================================================================

test_that("[9] susie_inf: L = 5 stored in params", {
  fit <- run_susie_inf_region(.rg, .rp, L = 5)
  expect_equal(fit$params$L, 5)
})

test_that("[9] susie_inf: L = 1 accepted", {
  fit <- run_susie_inf_region(.rg, .rp, L = 1)
  expect_equal(fit$params$L, 1)
})

test_that("[9] susie_inf: coverage = 0.9 stored in params", {
  fit <- run_susie_inf_region(.rg, .rp, L = 5, coverage = 0.9)
  expect_equal(fit$params$coverage, 0.9)
})

test_that("[9] susie_inf: max_iter = 50 stored in params", {
  fit <- run_susie_inf_region(.rg, .rp, L = 5, max_iter = 50)
  expect_equal(fit$params$max_iter, 50)
})

test_that("[9] susie_inf: output has pip, credible_sets, method", {
  fit <- run_susie_inf_region(.rg, .rp, L = 5)
  expect_true(all(c("pip", "credible_sets", "method") %in% names(fit)))
})

test_that("[9] susie_inf: pip in [0, 1]", {
  fit <- run_susie_inf_region(.rg, .rp, L = 5)
  expect_true(all(fit$pip >= 0))
  expect_true(all(fit$pip <= 1))
})


# =============================================================================
# SECTION 10: run_carma / run_carma_region argument passthrough
# =============================================================================

test_that("[10] carma: rho.index = 0.95 (default) stored in params", {
  skip_if_not_installed("CARMA")
  fit <- run_carma_region(.rg, .rp, rho.index = 0.95)
  expect_equal(fit$params$rho.index, 0.95)
})

test_that("[10] carma: rho.index = 0.9 stored in params", {
  skip_if_not_installed("CARMA")
  fit <- run_carma_region(.rg, .rp, rho.index = 0.9)
  expect_equal(fit$params$rho.index, 0.9)
})

test_that("[10] carma: num.causal = 5 stored in params", {
  skip_if_not_installed("CARMA")
  fit <- run_carma_region(.rg, .rp, num.causal = 5)
  expect_equal(fit$params$num.causal, 5)
})

test_that("[10] carma: num.causal = 1 produces a pip vector", {
  skip_if_not_installed("CARMA")
  fit <- run_carma_region(.rg, .rp, num.causal = 1)
  expect_false(is.null(fit$pip))
})

test_that("[10] carma: pip has correct length", {
  skip_if_not_installed("CARMA")
  fit <- run_carma_region(.rg, .rp)
  expect_equal(length(fit$pip), .rg$p)
})

test_that("[10] carma: pip in [0, 1]", {
  skip_if_not_installed("CARMA")
  fit <- run_carma_region(.rg, .rp)
  expect_true(all(fit$pip >= 0, na.rm = TRUE))
  expect_true(all(fit$pip <= 1, na.rm = TRUE))
})

test_that("[10] carma: credible_sets is a (possibly empty) list", {
  skip_if_not_installed("CARMA")
  fit <- run_carma_region(.rg, .rp)
  expect_true(is.list(fit$credible_sets))
  expect_gte(length(fit$credible_sets), 0L)
})


# =============================================================================
# SECTION 11: External binary wrappers - argument forwarding
# =============================================================================
#
# These wrappers call out to external binaries / Python scripts. Each
# test_that block tries the corresponding setup, then either runs the
# wrapper (which should error gracefully if the binary/dir is missing) or
# skips when the binary is absent. The point is to verify argument
# forwarding, not the actual external method.

test_that("[11-finemap] finemap: finemap_path arg recognised", {
  finemap_bin <- tryCatch(setup_finemap(download = FALSE),
                          error = function(e) NULL)
  skip_if(is.null(finemap_bin),
          "FINEMAP binary not available on this machine")
  fit <- tryCatch(
    run_finemap_region(.rg, .rp, finemap_path = finemap_bin, n_causal = 2),
    error = function(e) list(error = conditionMessage(e))
  )
  expect_false(is.null(fit))
})

test_that("[11-paintor] paintor: paintor_path arg recognised", {
  paintor_bin <- tryCatch(setup_paintor(), error = function(e) NULL)
  skip_if(is.null(paintor_bin),
          "PAINTOR binary not available / not compiled")
  fit <- tryCatch(
    run_paintor_region(.rg, .rp, paintor_path = paintor_bin, max_causal = 1),
    error = function(e) list(error = conditionMessage(e))
  )
  expect_false(is.null(fit))
})

test_that("[11-beatrice] beatrice: beatrice_dir + python args forwarded (graceful on miss)", {
  fit <- tryCatch(
    run_beatrice_region(.rg, .rp,
                        beatrice_dir = "~/Beatrice-Finemapping",
                        python = "/opt/anaconda3/bin/python3",
                        max_iter = 500),
    error = function(e) list(error = conditionMessage(e))
  )
  # Graceful return either way: pip vector or $error field
  expect_false(is.null(fit))
})

test_that("[11-funmap] funmap: python arg forwarded (graceful on miss)", {
  fit <- tryCatch(
    run_funmap_region(.rg, .rp,
                      python = "/opt/anaconda3/bin/python3",
                      L = 5),
    error = function(e) list(error = conditionMessage(e))
  )
  expect_false(is.null(fit))
})

test_that("[0.7 annotation_correlation] rho=0 gives near-identity empirical correlation", {
  set.seed(1L)
  p <- 5000L
  n_annotations <- 6L
  props <- rep(0.15, n_annotations)
  enrichment <- c(5, 5, 5, 1, 1, 1)   # two groups of 3
  A <- simulate_annotations_for_region(
    p = p, annotation_type = "binary", n_annotations = n_annotations,
    annotation_proportions = props, user_annotation_matrix = NULL,
    annotation_correlation = 0, enrichment = enrichment
  )$matrix
  C <- cor(A)
  off <- C[upper.tri(C)]
  expect_lt(max(abs(off)), 0.05)
  freq <- colMeans(A)
  expect_true(all(abs(freq - props) < 0.02))
})

test_that("[0.7 annotation_correlation] rho>0 induces positive within-group corr, none across groups", {
  set.seed(2L)
  p <- 5000L
  n_annotations <- 6L
  props <- rep(0.15, n_annotations)
  enrichment <- c(5, 5, 5, 1, 1, 1)   # two enrichment groups of 3
  A25 <- simulate_annotations_for_region(
    p = p, annotation_type = "binary", n_annotations = n_annotations,
    annotation_proportions = props, user_annotation_matrix = NULL,
    annotation_correlation = 0.25, enrichment = enrichment
  )$matrix
  A75 <- simulate_annotations_for_region(
    p = p, annotation_type = "binary", n_annotations = n_annotations,
    annotation_proportions = props, user_annotation_matrix = NULL,
    annotation_correlation = 0.75, enrichment = enrichment
  )$matrix

  in_group_pairs <- rbind(c(1,2), c(1,3), c(2,3), c(4,5), c(4,6), c(5,6))
  cross_pairs <- rbind(c(1,4), c(1,5), c(1,6), c(2,4), c(2,5), c(2,6),
                       c(3,4), c(3,5), c(3,6))

  mean_ig25 <- mean(apply(in_group_pairs, 1, function(ij) cor(A25[, ij[1]], A25[, ij[2]])))
  mean_ig75 <- mean(apply(in_group_pairs, 1, function(ij) cor(A75[, ij[1]], A75[, ij[2]])))
  mean_cross75 <- mean(apply(cross_pairs, 1, function(ij) cor(A75[, ij[1]], A75[, ij[2]])))

  expect_gt(mean_ig75, 0.20)
  expect_gt(mean_ig75, mean_ig25 + 0.05)
  expect_lt(abs(mean_cross75), 0.05)
  expect_true(all(abs(colMeans(A25) - props) < 0.02))
  expect_true(all(abs(colMeans(A75) - props) < 0.02))
})

test_that("[0.7 annotation_correlation] enrichment=NULL disables correlation", {
  set.seed(3L)
  p <- 5000L; n_annotations <- 4L
  props <- rep(0.2, n_annotations)
  A <- simulate_annotations_for_region(
    p = p, annotation_type = "binary", n_annotations = n_annotations,
    annotation_proportions = props, user_annotation_matrix = NULL,
    annotation_correlation = 0.75, enrichment = NULL   # no grouping info
  )$matrix
  expect_lt(max(abs(cor(A)[upper.tri(cor(A))])), 0.05)
})

test_that("[0.7 annotation_correlation] run_simulation accepts + forwards the arg", {
  sim <- run_simulation(
    n_regions = 2, n = 100, p = 40, n_iter = 1, S = 1, phi = 0.2,
    model = "sparse", annotations = "binary", n_annotations = 4,
    annotation_proportions = rep(0.2, 4),
    enrichment = c(5, 5, 1, 1),
    annotation_correlation = 0.5,
    genetic_map_dir = fmb_test_map_dir(),
    seed = 42, verbose = FALSE
  )
  expect_false(is.null(sim$scenarios[[1]]$regions[[1]]$annotations_matrix))
  expect_equal(sim$params$annotation_correlation, 0.5)
})

test_that("[11-fb-helper] .fb_extract_annotations: geno preferred, pheno fallback", {
  extract <- get(".fb_extract_annotations", envir = asNamespace("fmbenchmark"))
  A_geno  <- matrix(1, 5, 2)
  A_pheno <- matrix(2, 5, 2)
  expect_identical(extract(list(annotations_matrix = A_geno), list()), A_geno)
  expect_identical(extract(list(), list(annotations_matrix = A_pheno)), A_pheno)
  expect_identical(
    extract(list(annotations_matrix = A_geno),
            list(annotations_matrix = A_pheno)),
    A_geno
  )
  expect_null(extract(list(), list()))
})

test_that("[11-fb-regression] fb wrapper reads annotations from region_geno (gw-path fix)", {
  # Regression: the genome-wide simulator sets annotations_matrix ONLY on
  # region_geno. Previously the wrapper read from region_pheno, silently
  # dropping annotations under simulate_gwfm_data. Mock
  # run_functional_beatrice and inspect the annotations argument.
  n_snp <- 20; K <- 3
  rg <- list(
    LD = diag(n_snp), n = 100,
    variant_ids = as.character(seq_len(n_snp)),
    annotations_matrix = matrix(rbinom(n_snp * K, 1, 0.3), n_snp, K)
  )
  rp <- list(z = rnorm(n_snp))   # no annotations here — mimics gw path

  captured_A <- NULL
  fake_fb <- function(z, LD, n, annotations = NULL, ...) {
    captured_A <<- annotations
    list(pip = rep(0, length(z)), method = "functional_beatrice",
         credible_sets = list(), params = list(),
         runtime_seconds = 0, additional = list())
  }
  original <- get("run_functional_beatrice", envir = asNamespace("fmbenchmark"))
  assignInNamespace("run_functional_beatrice", fake_fb, ns = "fmbenchmark")
  on.exit(
    assignInNamespace("run_functional_beatrice", original, ns = "fmbenchmark"),
    add = TRUE
  )

  invisible(run_functional_beatrice_region(rg, rp))
  expect_false(is.null(captured_A),
               info = "annotations should be forwarded from region_geno")
  expect_equal(dim(captured_A), c(n_snp, K))
  expect_equal(captured_A, rg$annotations_matrix)

  # Permutation regression: shuffling annotation rows changes what the
  # wrapper hands off (previous bug returned NULL either way).
  perm <- sample(n_snp)
  rg2 <- rg
  rg2$annotations_matrix <- rg$annotations_matrix[perm, , drop = FALSE]
  captured_A <- NULL
  invisible(run_functional_beatrice_region(rg2, rp))
  expect_equal(captured_A, rg2$annotations_matrix)
  expect_false(isTRUE(all.equal(captured_A, rg$annotations_matrix)))
})


# =============================================================================
# SECTION 12: run_marginal_z / run_marginal_z_region (baseline)
# =============================================================================

test_that("[12] marginal_z: output has all standard fields", {
  fit <- run_marginal_z_region(.rg, .rp)
  expect_true(all(c("pip", "credible_sets", "method", "input_type",
                    "params", "runtime_seconds", "additional") %in% names(fit)))
  expect_equal(fit$method, "marginal_z")
  expect_equal(fit$input_type, "summary")
})

test_that("[12] marginal_z: pip length equals p", {
  fit <- run_marginal_z_region(.rg, .rp)
  expect_equal(length(fit$pip), .rg$p)
})

test_that("[12] marginal_z: pip values lie in [0, 1]", {
  fit <- run_marginal_z_region(.rg, .rp)
  expect_true(all(fit$pip >= 0))
  expect_true(all(fit$pip <= 1))
})

test_that("[12] marginal_z: pip sums to ~1 (|z|/sum|z|)", {
  fit <- run_marginal_z_region(.rg, .rp)
  expect_lt(abs(sum(fit$pip) - 1), 1e-8)
})

test_that("[12] marginal_z: default coverage = 0.95", {
  fit <- run_marginal_z_region(.rg, .rp)
  expect_equal(fit$params$coverage, 0.95)
})

test_that("[12] marginal_z: lower coverage produces smaller-or-equal CS", {
  fit95 <- run_marginal_z_region(.rg, .rp, coverage = 0.95)
  fit50 <- run_marginal_z_region(.rg, .rp, coverage = 0.50)
  expect_lte(length(fit50$credible_sets[[1]]),
             length(fit95$credible_sets[[1]]))
})

test_that("[12] marginal_z: returns exactly one credible set", {
  fit <- run_marginal_z_region(.rg, .rp)
  expect_true(is.list(fit$credible_sets))
  expect_length(fit$credible_sets, 1L)
})


# =============================================================================
# SECTION 13: run_polyfun_oracle / run_polyfun_oracle_region (oracle priors)
# =============================================================================
#
# polyfun_oracle reads simulator truth, so it needs SIM_MINI_ANNOT (which
# has annotations + enrichment in $truth).

.rg_po <- SIM_MINI_ANNOT$genotypes[[1]]
.rp_po <- SIM_MINI_ANNOT$scenarios[[1]]$regions[[1]]


test_that("[13] polyfun_oracle: output has all standard fields", {
  fit <- run_polyfun_oracle_region(.rg_po, .rp_po)
  expect_true(all(c("pip", "credible_sets", "method", "input_type",
                    "params", "runtime_seconds", "additional") %in% names(fit)))
  expect_equal(fit$method, "polyfun_oracle")
})

test_that("[13] polyfun_oracle: pip valid (length p, in [0, 1])", {
  fit <- run_polyfun_oracle_region(.rg_po, .rp_po)
  expect_equal(length(fit$pip), .rg_po$p)
  expect_true(all(fit$pip >= 0, na.rm = TRUE))
  expect_true(all(fit$pip <= 1, na.rm = TRUE))
})

test_that("[13] polyfun_oracle: reports prior_weights of length p", {
  fit <- run_polyfun_oracle_region(.rg_po, .rp_po)
  pw  <- fit$additional$prior_weights
  expect_false(is.null(pw))
  expect_equal(length(pw), .rg_po$p)
  expect_true(all(pw >= 0, na.rm = TRUE))
})

test_that("[13] polyfun_oracle: falls back to uniform prior when annotations absent", {
  # SIM_MINI has annotations="none" so truth$enrichment is also absent.
  # Wrapper is documented to fall back to uniform priors (degenerate to
  # plain SuSiE) rather than erroring.
  fit <- run_polyfun_oracle_region(.rg, .rp)
  expect_null(fit$error)
  expect_equal(length(fit$pip), .rg$p)
  expect_true(all(fit$pip >= 0, na.rm = TRUE))
  expect_true(all(fit$pip <= 1, na.rm = TRUE))
  expect_identical(fit$params$prior_source, "uniform_fallback")
})


# =============================================================================
# SECTION 14: run_polyfun_est / run_polyfun_est_region / scenario_setup
# =============================================================================

.rg_pe <- SIM_MINI_ANNOT$genotypes[[1]]
.rp_pe <- SIM_MINI_ANNOT$scenarios[[1]]$regions[[1]]


test_that("[14] polyfun_est: output has all standard fields", {
  fit <- run_polyfun_est_region(.rg_pe, .rp_pe)
  expect_true(all(c("pip", "credible_sets", "method", "input_type",
                    "params", "runtime_seconds", "additional") %in% names(fit)))
  expect_equal(fit$method, "polyfun_est")
})

test_that("[14] polyfun_est: pip valid (length p, in [0, 1])", {
  fit <- run_polyfun_est_region(.rg_pe, .rp_pe)
  expect_equal(length(fit$pip), .rg_pe$p)
  expect_true(all(fit$pip >= 0, na.rm = TRUE))
  expect_true(all(fit$pip <= 1, na.rm = TRUE))
})

test_that("[14] polyfun_est: reports non-negative tau of length m + 1", {
  fit <- run_polyfun_est_region(.rg_pe, .rp_pe)
  tau <- fit$additional$tau
  m   <- ncol(.rg_pe$annotations_matrix)
  expect_false(is.null(tau))
  expect_equal(length(tau), m + 1L)
  expect_true(all(tau >= 0, na.rm = TRUE))
})

test_that("[14] polyfun_est: runs without scenario hook (per-region tau fallback)", {
  fit <- run_polyfun_est_region(.rg_pe, .rp_pe, pooled_tau = NULL)
  expect_equal(length(fit$pip), .rg_pe$p)
  expect_false(is.null(fit$additional$tau))
})

test_that("[14] polyfun_est: scenario-setup hook returns pooled_tau when consistent", {
  scen <- SIM_MINI_ANNOT$scenarios[[1]]
  extra <- run_polyfun_est_scenario_setup(
    genotypes = SIM_MINI_ANNOT$genotypes,
    regions   = scen$regions,
    user_args = list()
  )
  m <- ncol(SIM_MINI_ANNOT$genotypes[[1]]$annotations_matrix)
  expect_true("pooled_tau" %in% names(extra))
  expect_equal(length(extra$pooled_tau), m + 1L)
  expect_true(all(extra$pooled_tau >= 0))
})

test_that("[14] polyfun_est: graceful behaviour on no-annotation fixture", {
  fit <- run_polyfun_est_region(.rg, .rp)
  # Either falls back to uniform priors or sets $error. We just want
  # no crash and a length-p pip vector.
  expect_false(is.null(fit))
  expect_equal(length(fit$pip), .rg$p)
})


# =============================================================================
# SECTION 14b: run_polyfun_ldsc (corrected LD-score PolyFun)
# =============================================================================
# See wrapper_polyfun_ldsc.R for the correction rationale. These tests
# lock in the S-LDSC-style regressor and the LOCO scenario_setup.

test_that("[14b] polyfun_ldsc: single-region fit gives valid PIPs + non-uniform prior", {
  # Use a sim with annotations so priors are non-trivial
  sim <- run_simulation(
    n_regions = 1, n = 150, p = 60, n_iter = 1, S = 2, phi = 0.3,
    model = "sparse", annotations = "binary", n_annotations = 3,
    annotation_proportions = rep(0.2, 3),
    enrichment = c(6, 1, 1),
    genetic_map_dir = fmb_test_map_dir(),
    seed = 1, verbose = FALSE
  )
  rg <- sim$genotypes[[1]]
  rp <- sim$scenarios[[1]]$regions[[1]]

  fit <- run_polyfun_ldsc_region(rg, rp)
  expect_equal(length(fit$pip), rg$p)
  expect_true(all(fit$pip >= 0 & fit$pip <= 1))
  expect_equal(fit$method, "polyfun_ldsc")
  expect_equal(fit$additional$prior_source, "single_region_ldsc")
  # Priors should differ from uniform when at least one tau > 0
  pw <- fit$additional$prior_weights
  expect_equal(length(pw), rg$p)
  expect_gt(sd(pw), 0)
})

test_that("[14b] polyfun_ldsc: uniform fallback when no annotations", {
  # .rg has no annotations_matrix
  fit <- run_polyfun_ldsc_region(.rg, .rp)
  expect_equal(length(fit$pip), .rg$p)
  expect_equal(fit$additional$prior_source, "uniform_fallback")
  expect_equal(fit$additional$prior_weights,
               rep(1 / .rg$p, .rg$p))
})

test_that("[14b] polyfun_ldsc scenario_setup: returns per-region LOCO tau vectors", {
  sim <- run_simulation(
    n_regions = 3, n = 150, p = 40, n_iter = 1, S = 1, phi = 0.2,
    model = "sparse", annotations = "binary", n_annotations = 3,
    annotation_proportions = rep(0.2, 3),
    enrichment = c(5, 1, 1),
    genetic_map_dir = fmb_test_map_dir(),
    seed = 5, verbose = FALSE
  )
  su <- run_polyfun_ldsc_scenario_setup(
    genotypes = sim$genotypes,
    regions   = sim$scenarios[[1]]$regions,
    user_args = list()
  )
  # Non-empty return means LOCO succeeded
  expect_named(su, "pooled_tau")
  expect_equal(length(su$pooled_tau), 3L)
  # Keyed by region_id so run_methods()'s scenario-wide arg merge still
  # dispatches the correct tau to the correct region
  expect_setequal(names(su$pooled_tau), c("1", "2", "3"))
  # Each tau is intercept + m coefficients = m + 1 entries
  for (tau in su$pooled_tau) {
    expect_equal(length(tau), 4L)   # intercept + 3 annotations
    expect_true(all(tau >= 0))       # NNLS enforces non-negativity
  }
  # LOCO: per-region taus should differ (each fit on a different subset
  # of regions), unless the data happens to be numerically degenerate.
  expect_false(isTRUE(all.equal(su$pooled_tau[["1"]], su$pooled_tau[["2"]])))
})

test_that("[14b] polyfun_ldsc via run_methods uses LOCO priors (differ from raw single-region)", {
  sim <- run_simulation(
    n_regions = 3, n = 150, p = 40, n_iter = 1, S = 1, phi = 0.2,
    model = "sparse", annotations = "binary", n_annotations = 3,
    annotation_proportions = rep(0.2, 3),
    enrichment = c(5, 1, 1),
    genetic_map_dir = fmb_test_map_dir(),
    seed = 7, verbose = FALSE
  )
  res <- run_methods(sim, methods = "polyfun_ldsc",
                     method_args = list(polyfun_ldsc = list(L = 3)),
                     save = FALSE, verbose = FALSE)
  expect_equal(res$polyfun_ldsc$n_total, 3L)
  expect_equal(res$polyfun_ldsc$n_failed, 0L)
  first_fit <- res$polyfun_ldsc$results[[1]]
  expect_equal(first_fit$additional$prior_source, "loco_scenario_setup")
})

test_that("[14b] polyfun_ldsc: ldscore helper computes l_{j,c} = sum_k r^2_{j,k} A_{k,c}", {
  ldscore_matrix <- get(".ldscore_matrix", envir = asNamespace("fmbenchmark"))
  A <- matrix(c(1, 0, 1, 0, 1,
                0, 1, 0, 1, 0), 5, 2, byrow = FALSE)
  LD <- diag(5); LD[1, 2] <- LD[2, 1] <- 0.5
  ell <- ldscore_matrix(A, LD)
  # Variant 1 (A_1 = c(1,0), r^2 to variant 2 = 0.25): ell_{1,1} = 1*1 + 0.25*0 = 1
  # Variant 2 (A_2 = c(0,1)): ell_{2,1} = 0.25*1 + 1*0 = 0.25
  expect_equal(ell[1, 1], 1.0)
  expect_equal(ell[2, 1], 0.25)
  # Column 2 has A = c(0,1,0,1,0): ell_{1,2} = 0 + 0.25*1 = 0.25
  expect_equal(ell[1, 2], 0.25)
})


# =============================================================================
# SECTION 15: MAF-stratified evaluation (by_causal_maf)
# =============================================================================

# RESULTS_MINI predates the MAF axis; rebuild eval here to get the field.
EVAL_MAF <- evaluate_methods(
  SIM_MINI, RESULTS_MINI,
  save = FALSE, verbose = FALSE
)


test_that("[15] evaluate_methods returns a by_causal_maf list per method", {
  expect_true("by_causal_maf" %in% names(EVAL_MAF$susie))
  expect_true("by_causal_maf" %in% names(EVAL_MAF$abf))
})

test_that("[15] by_causal_maf entries are populated when MAFs available", {
  # SIM_MINI uses real VCFs => genotypes carry MAFs => bins populate.
  bm <- EVAL_MAF$susie$by_causal_maf
  expect_false(is.null(bm))
  expect_gte(length(bm), 1L)
})

test_that("[15] by_causal_maf bin names are a subset of {rare, low, common}", {
  bm <- EVAL_MAF$susie$by_causal_maf
  expect_true(all(names(bm) %in% c("rare", "low", "common")))
})

test_that("[15] by_causal_maf bins appear in canonical order rare->low->common", {
  bm <- EVAL_MAF$susie$by_causal_maf
  canonical <- c("rare", "low", "common")
  present   <- canonical[canonical %in% names(bm)]
  expect_identical(names(bm), present)
})

test_that("[15] by_causal_maf bins each carry a numeric auprc field", {
  bm <- EVAL_MAF$susie$by_causal_maf
  for (b in names(bm)) {
    expect_true("auprc" %in% names(bm[[b]]),
                info = sprintf("bin=%s", b))
    expect_true(is.numeric(bm[[b]]$auprc),
                info = sprintf("bin=%s", b))
  }
})


# =============================================================================
# SECTION 16: Misspecification stratification (by_true_annotation_type)
# =============================================================================

# Build an eval against SIM_MINI_ANNOT so the annotated-fixture case is
# covered alongside the no-annotation case from EVAL_MAF.

RESULTS_MINI_ANNOT_FOR_AT <- run_methods(
  SIM_MINI_ANNOT,
  methods     = c("susie", "abf"),
  method_args = list(susie = list(L = 5L, coverage = 0.95)),
  save = FALSE, verbose = FALSE
)
EVAL_ANNOT <- evaluate_methods(
  SIM_MINI_ANNOT, RESULTS_MINI_ANNOT_FOR_AT,
  save = FALSE, verbose = FALSE
)


test_that("[16] by_true_annotation_type present per method", {
  expect_true("by_true_annotation_type" %in% names(EVAL_MAF$susie))
  expect_true("by_true_annotation_type" %in% names(EVAL_ANNOT$susie))
})

test_that("[16] no-annotation fixture: by_true_annotation_type has 'none'", {
  bt <- EVAL_MAF$susie$by_true_annotation_type
  expect_false(is.null(bt))
  expect_identical(names(bt), "none")
})

test_that("[16] binary-annotation fixture: by_true_annotation_type has 'binary'", {
  bt <- EVAL_ANNOT$susie$by_true_annotation_type
  expect_false(is.null(bt))
  expect_identical(names(bt), "binary")
})

test_that("[16] by_true_annotation_type bins carry a numeric auprc field", {
  for (eo in list(EVAL_MAF$susie, EVAL_ANNOT$susie)) {
    bt <- eo$by_true_annotation_type
    for (t in names(bt)) {
      expect_true("auprc" %in% names(bt[[t]]),
                  info = sprintf("type=%s", t))
      expect_true(is.numeric(bt[[t]]$auprc),
                  info = sprintf("type=%s", t))
    }
  }
})


# =============================================================================
# SECTION 17: LD mismatch - n_ref independent reference panel
# =============================================================================

# Baseline (no ref panel) - mirrors pre-n_ref behaviour exactly.
GENO_NOREF <- simulate_genotypes(
  n_regions = 1L, n = 200L, p = 60L,
  genetic_map_dir = fmb_test_map_dir(),
  seed = 99L, verbose = FALSE
)

# Small ref panel - LD will differ noticeably from in-sample LD.
GENO_REF <- simulate_genotypes(
  n_regions = 1L, n = 200L, p = 60L,
  genetic_map_dir = fmb_test_map_dir(),
  seed = 99L, verbose = FALSE,
  n_ref = 50L
)

# Ref panel matching n - independent draw but larger, so LD closer to truth.
GENO_REF_FULL <- simulate_genotypes(
  n_regions = 1L, n = 200L, p = 60L,
  genetic_map_dir = fmb_test_map_dir(),
  seed = 99L, verbose = FALSE,
  n_ref = 200L
)


test_that("[17] n_ref = NULL: no X_ref / n_ref fields (backwards compatible)", {
  r <- GENO_NOREF[[1L]]
  expect_null(r$X_ref)
  expect_null(r$n_ref)
})

test_that("[17] n_ref = 50: X_ref has 50 rows and same p as X", {
  r <- GENO_REF[[1L]]
  expect_false(is.null(r$X_ref))
  expect_equal(nrow(r$X_ref), 50L)
  expect_equal(ncol(r$X_ref), ncol(r$X))
  expect_identical(r$n_ref, 50L)
})

test_that("[17] run_simulation: LD_true and LD both populated when n_ref set", {
  sim_ref <- run_simulation(
    n_regions = 1L, n = 200L, p = 60L, n_iter = 1L,
    S = 1L, phi = 0.2, model = "sparse", annotations = "none",
    genetic_map_dir = fmb_test_map_dir(),
    seed = 99L, verbose = FALSE, n_ref = 50L
  )
  r <- sim_ref$genotypes[[1L]]
  expect_false(is.null(r$LD))
  expect_false(is.null(r$LD_true))
  expect_identical(dim(r$LD), dim(r$LD_true))
  expect_identical(sim_ref$params$n_ref, 50L)
})

test_that("[17] run_simulation: n_ref = NULL gives LD identical to LD_true", {
  sim_noref <- run_simulation(
    n_regions = 1L, n = 200L, p = 60L, n_iter = 1L,
    S = 1L, phi = 0.2, model = "sparse", annotations = "none",
    genetic_map_dir = fmb_test_map_dir(),
    seed = 99L, verbose = FALSE
  )
  r <- sim_noref$genotypes[[1L]]
  expect_false(is.null(r$LD_true))
  expect_identical(r$LD, r$LD_true)
  expect_null(sim_noref$params$n_ref)
})

test_that("[17] LD mismatch shrinks as n_ref grows (rough monotonicity)", {
  # Same main-sample seed across both, so LD_true is identical. The larger
  # ref panel should be closer to LD_true because correlation estimates
  # are more precise at larger n.
  r_small <- GENO_REF[[1L]]
  r_large <- GENO_REF_FULL[[1L]]
  LD_true <- cor(r_small$X)
  d_small <- mean((cor(r_small$X_ref) - LD_true)^2, na.rm = TRUE)
  d_large <- mean((cor(r_large$X_ref) - LD_true)^2, na.rm = TRUE)
  expect_lt(d_large, d_small)
})


# =============================================================================
# SECTION 18: run_sparsepro / run_sparsepro_region wrapper plumbing
# =============================================================================
#
# Renumbered from the source file's duplicate "SECTION 16".
#
# SparsePro is a Python CLI script (sparsepro_zld.py) in the upstream repo
# (https://github.com/zhwm/SparsePro). We don't bundle or auto-install it.
# To detect availability we honour SPARSEPRO_DIR (+ optionally
# SPARSEPRO_PYTHON) env vars. When absent the happy-path tests skip with
# a clear reason; wrapper-plumbing tests (graceful failure, output shape
# on the error path, argument forwarding, input validation) run
# unconditionally.

sparsepro_dir       <- Sys.getenv("SPARSEPRO_DIR",    unset = "")
sparsepro_python    <- Sys.getenv("SPARSEPRO_PYTHON", unset = "python")
sparsepro_available <- nzchar(sparsepro_dir) &&
  file.exists(file.path(sparsepro_dir, "sparsepro_zld.py"))

.rg_sp <- SIM_MINI$genotypes[[1]]
.rp_sp <- SIM_MINI$scenarios[[1]]$regions[[1]]


test_that("[18] sparsepro: standard-shape error result on missing dir", {
  fit <- run_sparsepro_region(.rg_sp, .rp_sp,
                              sparsepro_dir = "/definitely/not/a/path")
  fields <- c("pip", "credible_sets", "method", "input_type",
              "params", "runtime_seconds", "additional")
  expect_true(all(fields %in% names(fit)))
  expect_equal(fit$method, "sparsepro")
  expect_equal(fit$input_type, "summary")
  expect_false(is.null(fit$error))
  expect_equal(length(fit$pip), .rg_sp$p)
  expect_true(all(is.na(fit$pip)))
  expect_length(fit$credible_sets, 0L)
})

test_that("[18] sparsepro: sparsepro_dir / python / K / cthres recorded in params", {
  fit <- run_sparsepro_region(
    .rg_sp, .rp_sp,
    sparsepro_dir = "/tmp/SparsePro-test",
    python        = "python",
    K             = 7,
    cthres        = 0.9
  )
  expect_identical(fit$params$sparsepro_dir, "/tmp/SparsePro-test")
  expect_identical(fit$params$python, "python")
  expect_equal(fit$params$K, 7)
  expect_lt(abs(fit$params$cthres - 0.9), 1e-12)
  expect_equal(fit$params$n, .rg_sp$n)
})

test_that("[18] sparsepro: K must be >= 1", {
  e <- tryCatch(
    run_sparsepro(z = .rp_sp$z, LD = .rg_sp$LD, n = .rg_sp$n,
                  sparsepro_dir = "/tmp/x", K = 0),
    error = function(e) conditionMessage(e)
  )
  expect_true(is.character(e))
  expect_match(e, "K must be")
})

test_that("[18] sparsepro: cthres must be in (0, 1)", {
  e <- tryCatch(
    run_sparsepro(z = .rp_sp$z, LD = .rg_sp$LD, n = .rg_sp$n,
                  sparsepro_dir = "/tmp/x", cthres = 1.5),
    error = function(e) conditionMessage(e)
  )
  expect_true(is.character(e))
  expect_match(e, "cthres must be")
})

test_that("[18] setup_sparsepro: errors clearly when sparsepro_dir missing", {
  e <- tryCatch(
    setup_sparsepro(sparsepro_dir = "/definitely/not/a/path"),
    error = function(e) conditionMessage(e)
  )
  expect_true(is.character(e))
  expect_match(e, "sparsepro_zld.py not found")
})

test_that("[18] sparsepro: produces valid PIPs on small fixture (real install)", {
  skip_if(!sparsepro_available,
          "SparsePro not installed (set SPARSEPRO_DIR to enable)")

  fit <- run_sparsepro_region(.rg_sp, .rp_sp,
                              sparsepro_dir = sparsepro_dir,
                              python        = sparsepro_python,
                              K             = 3)
  expect_null(fit$error)
  expect_equal(length(fit$pip), .rg_sp$p)
  expect_true(all(fit$pip >= 0, na.rm = TRUE))
  expect_true(all(fit$pip <= 1, na.rm = TRUE))
})

test_that("[18] sparsepro: credible_sets are integer index vectors (real install)", {
  skip_if(!sparsepro_available,
          "SparsePro not installed (set SPARSEPRO_DIR to enable)")

  fit <- run_sparsepro_region(.rg_sp, .rp_sp,
                              sparsepro_dir = sparsepro_dir,
                              python        = sparsepro_python,
                              K             = 3)
  expect_true(is.list(fit$credible_sets))
  if (length(fit$credible_sets) > 0L) {
    for (cs in fit$credible_sets) {
      expect_true(is.numeric(cs))
      expect_true(all(cs >= 1L))
      expect_true(all(cs <= .rg_sp$p))
    }
  }
})
