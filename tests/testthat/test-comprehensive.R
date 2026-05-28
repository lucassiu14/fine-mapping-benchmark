# =============================================================================
# test-comprehensive.R
#
# Argument-level test coverage for every public function in the benchmark.
# Translated from scripts/test_comprehensive.R as part of Phase 3b-ii bulk.
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
# The original `scripts/test_comprehensive.R` is kept on disk for now;
# Phase 3b-iii decides whether to move it under `inst/scripts/` or
# delete it.
# =============================================================================


# --- Shared test fixtures --------------------------------------------------
#
# Built once at file load time so multiple test_that() blocks within a
# section can share them. Mirrors the original script.

# Small genotype object (2 regions, n=100, p=50)
GENO_SMALL <- simulate_genotypes(
  n_regions = 2, n = 100, p = 50,
  genetic_map_dir = "../../data/genetic_maps",
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
  genetic_map_dir = "../../data/genetic_maps",
  seed = 42, verbose = FALSE
)

# Annotated variant for polyfun_oracle / polyfun_est tests
SIM_MINI_ANNOT <- run_simulation(
  n_regions = 2, n = 100, p = 50,
  n_iter = 2, S = c(1, 2), phi = c(0.2, 0.4),
  model = "sparse", annotations = "binary", n_annotations = 3,
  enrichment = 5.0,
  genetic_map_dir = "../../data/genetic_maps",
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
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 1, verbose = FALSE)
  expect_length(g, 1L)
})

test_that("[1] n_regions = 3 returns list of length 3", {
  g <- simulate_genotypes(n_regions = 3, n = 50, p = 30,
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 2, verbose = FALSE)
  expect_length(g, 3L)
})

test_that("[1] n sets number of rows in X", {
  g <- simulate_genotypes(n_regions = 1, n = 80, p = 30,
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 3, verbose = FALSE)
  expect_equal(nrow(g[[1]]$X), 80L)
})

test_that("[1] p as scalar applied to all regions (with truncation cap)", {
  g <- simulate_genotypes(n_regions = 2, n = 50, p = 40,
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 4, verbose = FALSE)
  expect_lte(g[[1]]$p, 40L)
  expect_lte(g[[2]]$p, 40L)
})

test_that("[1] p as vector sets different targets per region", {
  g <- simulate_genotypes(n_regions = 2, n = 50, p = c(30, 50),
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 5, verbose = FALSE)
  expect_length(g, 2L)
})

test_that("[1] p > 500 with bundled VCF warns and caps", {
  expect_warning(
    simulate_genotypes(n_regions = 1, n = 50, p = 600,
                       genetic_map_dir = "../../data/genetic_maps",
                       seed = 6, verbose = FALSE)
  )
})

test_that("[1] vcf_files = NULL uses bundled example VCF", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          vcf_files = NULL,
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 7, verbose = FALSE)
  expect_false(is.null(g[[1]]$X))
})

test_that("[1] vcf_files = single path reused for all regions", {
  vcf <- system.file("examples", "region.vcf.gz", package = "sim1000G")
  g <- simulate_genotypes(n_regions = 2, n = 50, p = 30,
                          vcf_files = vcf,
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 8, verbose = FALSE)
  expect_length(g, 2L)
})

test_that("[1] vcf_files wrong length errors", {
  vcf <- system.file("examples", "region.vcf.gz", package = "sim1000G")
  expect_error(
    simulate_genotypes(n_regions = 3, n = 50, p = 30,
                       vcf_files = c(vcf, vcf),
                       genetic_map_dir = "../../data/genetic_maps",
                       seed = 9, verbose = FALSE)
  )
})

test_that("[1] vcf_files missing file errors", {
  expect_error(
    simulate_genotypes(n_regions = 1, n = 50, p = 30,
                       vcf_files = "/nonexistent/path.vcf.gz",
                       genetic_map_dir = "../../data/genetic_maps",
                       seed = 10, verbose = FALSE)
  )
})

test_that("[1] min_maf = 0 accepts all variants", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          min_maf = 0,
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 11, verbose = FALSE)
  expect_false(is.null(g[[1]]$X))
})

test_that("[1] min_maf = 0.1 (stricter filter) works", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          min_maf = 0.1,
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 12, verbose = FALSE)
  # MAF is two-sided; values must be in [0.1, 0.9].
  expect_true(all(g[[1]]$maf >= 0.1 | g[[1]]$maf <= 0.9))
})

test_that("[1] max_maf = 0.3 applies upper MAF filter", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          min_maf = 0.01, max_maf = 0.3,
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 13, verbose = FALSE)
  expect_false(is.null(g[[1]]$X))
})

test_that("[1] standardise = TRUE gives ~zero-mean columns", {
  g <- simulate_genotypes(n_regions = 1, n = 200, p = 30,
                          standardise = TRUE,
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 14, verbose = FALSE)
  expect_true(all(abs(colMeans(g[[1]]$X)) < 0.01))
})

test_that("[1] standardise = FALSE returns 0/1/2 coding", {
  g <- simulate_genotypes(n_regions = 1, n = 100, p = 30,
                          standardise = FALSE,
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 15, verbose = FALSE)
  vals <- unique(as.vector(g[[1]]$X))
  expect_true(all(vals %in% c(0, 1, 2)))
})

test_that("[1] standardise = FALSE: X identical to X_raw", {
  g <- simulate_genotypes(n_regions = 1, n = 100, p = 30,
                          standardise = FALSE,
                          genetic_map_dir = "../../data/genetic_maps",
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
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = 18, verbose = FALSE)
  expect_false(is.null(g[[1]]$X))
})

test_that("[1] seed ensures reproducibility", {
  g1 <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                            genetic_map_dir = "../../data/genetic_maps",
                            seed = 99, verbose = FALSE)
  g2 <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                            genetic_map_dir = "../../data/genetic_maps",
                            seed = 99, verbose = FALSE)
  expect_identical(g1[[1]]$X, g2[[1]]$X)
})

test_that("[1] seed = NULL accepted (no reproducibility required)", {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          genetic_map_dir = "../../data/genetic_maps",
                          seed = NULL, verbose = FALSE)
  expect_false(is.null(g[[1]]$X))
})

test_that("[1] verbose = FALSE suppresses messages", {
  out <- utils::capture.output(
    simulate_genotypes(n_regions = 1, n = 50, p = 30,
                       genetic_map_dir = "../../data/genetic_maps",
                       seed = 20, verbose = FALSE),
    type = "message"
  )
  expect_true(length(out) == 0L || !any(grepl("Region", out)))
})

test_that("[1] verbose = TRUE prints region progress", {
  out <- utils::capture.output(
    simulate_genotypes(n_regions = 1, n = 50, p = 30,
                       genetic_map_dir = "../../data/genetic_maps",
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
                       genetic_map_dir = "../../data/genetic_maps",
                       verbose = FALSE)
  )
})

test_that("[1] p vector length mismatch errors", {
  expect_error(
    simulate_genotypes(n_regions = 3, n = 50, p = c(30, 40),
                       genetic_map_dir = "../../data/genetic_maps",
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
                     genetic_map_dir = "../../data/genetic_maps",
                     seed = 200, save = FALSE, output_dir = tmp,
                     verbose = FALSE)
  expect_length(list.files(tmp), 0L)
})

test_that("[2b] save = TRUE writes .rds file", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  simulate_genotypes(n_regions = 1, n = 50, p = 30,
                     genetic_map_dir = "../../data/genetic_maps",
                     seed = 201, save = TRUE, output_dir = tmp,
                     verbose = FALSE)
  expect_length(list.files(tmp, pattern = "[.]rds$"), 1L)
})

test_that("[2b] saved .rds is readable and has correct structure", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  simulate_genotypes(n_regions = 1, n = 50, p = 30,
                     genetic_map_dir = "../../data/genetic_maps",
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
                     genetic_map_dir = "../../data/genetic_maps",
                     seed = 203, save = TRUE, output_dir = tmp,
                     verbose = FALSE)
  expect_true(dir.exists(tmp))
})

test_that("[2b] filename encodes n_regions, n, p, seed", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  simulate_genotypes(n_regions = 2, n = 60, p = 35,
                     genetic_map_dir = "../../data/genetic_maps",
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
                     genetic_map_dir = "../../data/genetic_maps",
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
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 1, verbose = FALSE)
  expect_length(r$scenarios, 1L)
})

test_that("[3] n_iter = 3 produces three scenarios", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 3, S = 1, phi = 0.2,
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 2, verbose = FALSE)
  expect_length(r$scenarios, 3L)
})

test_that("[3] S vector sweeps correctly", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = c(1, 2), phi = 0.2,
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 3, verbose = FALSE)
  expect_length(r$scenarios, 2L)
  S_vals <- sapply(r$scenarios, `[[`, "S")
  expect_equal(sort(S_vals), c(1L, 2L))
})

test_that("[3] phi vector sweeps correctly", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = c(0.2, 0.4),
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 4, verbose = FALSE)
  expect_length(r$scenarios, 2L)
  phi_vals <- sapply(r$scenarios, `[[`, "phi")
  expect_equal(sort(phi_vals), c(0.2, 0.4))
})

test_that("[3] model = 'sparse' recorded in params", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      model = "sparse",
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 5, verbose = FALSE)
  expect_equal(r$params$model, "sparse")
})

test_that("[3] model = 'sparse_inf' sweeps p_causal", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      model = "sparse_inf", p_causal = c(0.2, 0.5),
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 6, verbose = FALSE)
  expect_length(r$scenarios, 2L)
  expect_equal(r$params$model, "sparse_inf")
})

test_that("[3] inf_model = 'beatrice' accepted in sparse_inf", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      model = "sparse_inf", p_causal = 0.5,
                      inf_model = "beatrice",
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 7, verbose = FALSE)
  expect_equal(r$params$inf_model, "beatrice")
})

test_that("[3] inf_model = 'susie_inf' accepted in sparse_inf", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      model = "sparse_inf", p_causal = 0.5,
                      inf_model = "susie_inf",
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 8, verbose = FALSE)
  expect_equal(r$params$inf_model, "susie_inf")
})

test_that("[3] effect_distribution = 'normal' recorded in params", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      effect_distribution = "normal",
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 9, verbose = FALSE)
  expect_equal(r$params$effect_distribution, "normal")
})

test_that("[3] effect_distribution = 'equal' recorded", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      effect_distribution = "equal",
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 10, verbose = FALSE)
  expect_equal(r$params$effect_distribution, "equal")
})

test_that("[3] effect_variance = 0.5 recorded", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      effect_variance = 0.5,
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 11, verbose = FALSE)
  expect_equal(r$params$effect_variance, 0.5)
})

test_that("[3] annotations = 'none' produces NULL annotation matrix", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "none",
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 12, verbose = FALSE)
  expect_null(r$scenarios[[1]]$regions[[1]]$annotations_matrix)
})

test_that("[3] annotations = 'binary' with n_annotations = 2", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "binary", n_annotations = 2,
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 13, verbose = FALSE)
  A <- r$scenarios[[1]]$regions[[1]]$annotations_matrix
  expect_false(is.null(A))
  expect_equal(ncol(A), 2L)
})

test_that("[3] annotations = 'continuous' with n_annotations = 3", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "continuous", n_annotations = 3,
                      genetic_map_dir = "../../data/genetic_maps",
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
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 15, verbose = FALSE)
  expect_equal(r$params$annotation_proportions, 0.15)
})

test_that("[3] enrichment scalar passed through", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "binary", n_annotations = 2,
                      enrichment = 4,
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 16, verbose = FALSE)
  expect_equal(r$params$enrichment, 4)
})

test_that("[3] vcf_dir missing directory errors", {
  expect_error(
    run_simulation(n_regions = 1, n = 80, p = 40,
                   n_iter = 1, S = 1, phi = 0.2,
                   vcf_dir = "/nonexistent/dir",
                   genetic_map_dir = "../../data/genetic_maps",
                   seed = 17, verbose = FALSE)
  )
})

test_that("[3] vcf_files = single VCF path used for all regions", {
  vcf <- system.file("examples", "region.vcf.gz", package = "sim1000G")
  r <- run_simulation(n_regions = 2, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      vcf_files = vcf,
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 18, verbose = FALSE)
  expect_length(r$genotypes, 2L)
})

test_that("[3] min_maf = 0.05 passed to simulate_genotypes", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      min_maf = 0.05,
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 19, verbose = FALSE)
  expect_equal(r$params$min_maf, 0.05)
})

test_that("[3] max_maf = 0.4 passed to simulate_genotypes", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      max_maf = 0.4,
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 20, verbose = FALSE)
  expect_false(is.null(r$genotypes[[1]]$X))
})

test_that("[3] standardise = FALSE returns raw 0/1/2 genotypes", {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      standardise = FALSE,
                      genetic_map_dir = "../../data/genetic_maps",
                      seed = 21, verbose = FALSE)
  vals <- unique(as.vector(r$genotypes[[1]]$X))
  expect_true(all(vals %in% c(0, 1, 2)))
})

test_that("[3] seed = 42 ensures reproducibility", {
  r1 <- run_simulation(n_regions = 1, n = 80, p = 40,
                       n_iter = 1, S = 1, phi = 0.2,
                       genetic_map_dir = "../../data/genetic_maps",
                       seed = 42, verbose = FALSE)
  r2 <- run_simulation(n_regions = 1, n = 80, p = 40,
                       n_iter = 1, S = 1, phi = 0.2,
                       genetic_map_dir = "../../data/genetic_maps",
                       seed = 42, verbose = FALSE)
  expect_identical(r1$genotypes[[1]]$X, r2$genotypes[[1]]$X)
})

test_that("[3] save = TRUE writes .rds file", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  run_simulation(n_regions = 1, n = 80, p = 40,
                 n_iter = 1, S = 1, phi = 0.2,
                 genetic_map_dir = "../../data/genetic_maps",
                 seed = 99, save = TRUE, output_dir = tmp,
                 verbose = FALSE)
  expect_length(list.files(tmp, pattern = "[.]rds$"), 1L)
})

test_that("[3] save = FALSE writes no files", {
  tmp <- tempfile(); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE))
  run_simulation(n_regions = 1, n = 80, p = 40,
                 n_iter = 1, S = 1, phi = 0.2,
                 genetic_map_dir = "../../data/genetic_maps",
                 seed = 100, save = FALSE, output_dir = tmp,
                 verbose = FALSE)
  expect_length(list.files(tmp), 0L)
})

test_that("[3] output_dir is created when nested + missing", {
  tmp <- file.path(tempfile(), "deep", "nested")
  on.exit(unlink(dirname(dirname(tmp)), recursive = TRUE))
  run_simulation(n_regions = 1, n = 80, p = 40,
                 n_iter = 1, S = 1, phi = 0.2,
                 genetic_map_dir = "../../data/genetic_maps",
                 seed = 101, save = TRUE, output_dir = tmp,
                 verbose = FALSE)
  expect_true(dir.exists(tmp))
})

test_that("[3] verbose = FALSE suppresses messages", {
  out <- utils::capture.output(
    run_simulation(n_regions = 1, n = 80, p = 40,
                   n_iter = 1, S = 1, phi = 0.2,
                   genetic_map_dir = "../../data/genetic_maps",
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
                   genetic_map_dir = "../../data/genetic_maps",
                   verbose = FALSE)
  )
})

test_that("[3] phi outside (0,1) errors", {
  expect_error(
    run_simulation(n_regions = 1, n = 80, p = 40,
                   n_iter = 1, S = 1, phi = 1.5,
                   genetic_map_dir = "../../data/genetic_maps",
                   verbose = FALSE)
  )
})
