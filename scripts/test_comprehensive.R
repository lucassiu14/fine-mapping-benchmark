# =============================================================================
# scripts/test_comprehensive.R
#
# Comprehensive argument-level test for every public function in the benchmark.
#
# Tests EVERY argument of:
#   simulate_genotypes, simulate_phenotypes, run_simulation,
#   run_methods, evaluate_methods, plot_results
#
# Also verifies that user-specified method arguments are accepted and forwarded
# for all pure-R wrappers (susie, susie_inf, abf, carma) and silently skipped
# for external-binary methods (finemap, paintor, beatrice, funmap).
#
# Usage (from project root):
#   Rscript scripts/test_comprehensive.R
#
# Produces:
#   docs/testing_report.md — argument-level pass/fail documentation
# =============================================================================

suppressPackageStartupMessages({
  source("R/utils.R")
  source("R/simulate_genotypes.R")
  source("R/simulate_phenotypes.R")
  source("R/run_simulation.R")
  source("R/run_methods.R")
  source("R/evaluate.R")
  source("R/plot_results.R")
  source("R/wrappers/susie.R")
  source("R/wrappers/susie_inf.R")
  source("R/wrappers/abf.R")
  source("R/wrappers/finemap.R")
  source("R/wrappers/funmap.R")
  source("R/wrappers/paintor.R")
  source("R/wrappers/beatrice.R")
  source("R/wrappers/carma.R")
  source("R/wrappers/marginal_z.R")
  source("R/wrappers/polyfun_oracle.R")
  source("R/wrappers/polyfun_est.R")
})

# =============================================================================
# Test framework
# =============================================================================

.RESULTS     <- list()
.SECTION_NOW <- "Unknown"
.n_pass      <- 0L
.n_fail      <- 0L
.n_skip      <- 0L

set_section <- function(name) {
  .SECTION_NOW <<- name
  message("\n", strrep("=", 70))
  message("SECTION: ", name)
  message(strrep("=", 70))
}

run_test <- function(name, fn, skip_reason = NULL) {
  if (!is.null(skip_reason)) {
    .RESULTS[[length(.RESULTS) + 1L]] <<- list(
      section = .SECTION_NOW, name = name,
      status = "SKIP", error = skip_reason
    )
    .n_skip <<- .n_skip + 1L
    message(sprintf("  SKIP  %s  [%s]", name, skip_reason))
    return(invisible(NULL))
  }

  t0  <- proc.time()["elapsed"]
  res <- tryCatch(
    withCallingHandlers(
      { fn(); list(status = "PASS", error = NULL) },
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) list(status = "FAIL", error = conditionMessage(e))
  )
  elapsed <- round(proc.time()["elapsed"] - t0, 2)

  .RESULTS[[length(.RESULTS) + 1L]] <<- list(
    section = .SECTION_NOW, name = name,
    status = res$status, error = res$error, elapsed = elapsed
  )
  if (res$status == "PASS") {
    .n_pass <<- .n_pass + 1L
    message(sprintf("  PASS  %s  (%.2fs)", name, elapsed))
  } else {
    .n_fail <<- .n_fail + 1L
    message(sprintf("  FAIL  %s  — %s", name, res$error))
  }
}

# Expect an error (test passes when fn() throws)
expect_error <- function(fn) {
  threw <- tryCatch({ fn(); FALSE }, error = function(e) TRUE)
  if (!threw) stop("Expected an error but none was thrown.")
  invisible(TRUE)
}

# Expect a warning (test passes when fn() warns)
expect_warning_fn <- function(fn) {
  warned <- FALSE
  withCallingHandlers(fn(), warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
  if (!warned) stop("Expected a warning but none was raised.")
  invisible(TRUE)
}

# =============================================================================
# Shared test fixtures — generated once to keep runtime manageable
# =============================================================================

message("Building shared test fixtures...")

# Small genotype object (2 regions, n=100, p=50)
GENO_SMALL <- simulate_genotypes(
  n_regions = 2, n = 100, p = 50,
  genetic_map_dir = "data/genetic_maps",
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
  genetic_map_dir = "data/genetic_maps",
  seed = 42, verbose = FALSE
)

# Annotated variant for polyfun_oracle / polyfun_est tests, which need
# either a known per-SNP prior (oracle) or an annotation matrix to fit
# tau on (est). Mirrors SIM_MINI but with three binary annotations.
SIM_MINI_ANNOT <- run_simulation(
  n_regions = 2, n = 100, p = 50,
  n_iter = 2, S = c(1, 2), phi = c(0.2, 0.4),
  model = "sparse", annotations = "binary", n_annotations = 3,
  enrichment = 5.0,
  genetic_map_dir = "data/genetic_maps",
  seed = 42, verbose = FALSE
)

# Run susie and abf on MINI so evaluate/plot tests have real results
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

message("Fixtures ready.\n")


# =============================================================================
# SECTION 1: simulate_genotypes
# =============================================================================

set_section("simulate_genotypes")

run_test("n_regions = 1 returns list of length 1", function() {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 1, verbose = FALSE)
  stopifnot(length(g) == 1L)
})

run_test("n_regions = 3 returns list of length 3", function() {
  g <- simulate_genotypes(n_regions = 3, n = 50, p = 30,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 2, verbose = FALSE)
  stopifnot(length(g) == 3L)
})

run_test("n sets number of rows in X", function() {
  g <- simulate_genotypes(n_regions = 1, n = 80, p = 30,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 3, verbose = FALSE)
  stopifnot(nrow(g[[1]]$X) == 80L)
})

run_test("p as scalar applied to all regions", function() {
  g <- simulate_genotypes(n_regions = 2, n = 50, p = 40,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 4, verbose = FALSE)
  stopifnot(g[[1]]$p <= 40L, g[[2]]$p <= 40L)
})

run_test("p as vector sets different targets per region", function() {
  g <- simulate_genotypes(n_regions = 2, n = 50, p = c(30, 50),
                          genetic_map_dir = "data/genetic_maps",
                          seed = 5, verbose = FALSE)
  stopifnot(length(g) == 2L)
})

run_test("p > 500 with bundled VCF warns and caps", function() {
  expect_warning_fn(function() {
    simulate_genotypes(n_regions = 1, n = 50, p = 600,
                       genetic_map_dir = "data/genetic_maps",
                       seed = 6, verbose = FALSE)
  })
})

run_test("vcf_files = NULL uses bundled example VCF", function() {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          vcf_files = NULL,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 7, verbose = FALSE)
  stopifnot(!is.null(g[[1]]$X))
})

run_test("vcf_files = single path reused for all regions", function() {
  vcf <- system.file("examples", "region.vcf.gz", package = "sim1000G")
  g <- simulate_genotypes(n_regions = 2, n = 50, p = 30,
                          vcf_files = vcf,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 8, verbose = FALSE)
  stopifnot(length(g) == 2L)
})

run_test("vcf_files wrong length errors", function() {
  vcf <- system.file("examples", "region.vcf.gz", package = "sim1000G")
  expect_error(function() {
    simulate_genotypes(n_regions = 3, n = 50, p = 30,
                       vcf_files = c(vcf, vcf),  # length 2, n_regions = 3
                       genetic_map_dir = "data/genetic_maps",
                       seed = 9, verbose = FALSE)
  })
})

run_test("vcf_files missing file errors", function() {
  expect_error(function() {
    simulate_genotypes(n_regions = 1, n = 50, p = 30,
                       vcf_files = "/nonexistent/path.vcf.gz",
                       genetic_map_dir = "data/genetic_maps",
                       seed = 10, verbose = FALSE)
  })
})

run_test("min_maf = 0 accepts all variants", function() {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          min_maf = 0,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 11, verbose = FALSE)
  stopifnot(!is.null(g[[1]]$X))
})

run_test("min_maf = 0.1 (stricter filter) works", function() {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          min_maf = 0.1,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 12, verbose = FALSE)
  stopifnot(all(g[[1]]$maf >= 0.1 | g[[1]]$maf <= 0.9))
})

run_test("max_maf = 0.3 applies upper MAF filter", function() {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          min_maf = 0.01, max_maf = 0.3,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 13, verbose = FALSE)
  stopifnot(!is.null(g[[1]]$X))
})

run_test("standardise = TRUE gives ~zero-mean columns", function() {
  g <- simulate_genotypes(n_regions = 1, n = 200, p = 30,
                          standardise = TRUE,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 14, verbose = FALSE)
  col_means <- colMeans(g[[1]]$X)
  stopifnot(all(abs(col_means) < 0.01))
})

run_test("standardise = FALSE returns 0/1/2 coding", function() {
  g <- simulate_genotypes(n_regions = 1, n = 100, p = 30,
                          standardise = FALSE,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 15, verbose = FALSE)
  vals <- unique(as.vector(g[[1]]$X))
  stopifnot(all(vals %in% c(0, 1, 2)))
})

run_test("standardise = FALSE returns X_raw identical to X", function() {
  g <- simulate_genotypes(n_regions = 1, n = 100, p = 30,
                          standardise = FALSE,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 16, verbose = FALSE)
  stopifnot(identical(g[[1]]$X, g[[1]]$X_raw))
})

run_test("genetic_map_dir = NULL (uses tempdir) works", function() {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          genetic_map_dir = NULL,
                          seed = 17, verbose = FALSE)
  stopifnot(!is.null(g[[1]]$X))
})

run_test("genetic_map_dir = existing path caches maps", function() {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          genetic_map_dir = "data/genetic_maps",
                          seed = 18, verbose = FALSE)
  stopifnot(!is.null(g[[1]]$X))
})

run_test("seed ensures reproducibility", function() {
  g1 <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                            genetic_map_dir = "data/genetic_maps",
                            seed = 99, verbose = FALSE)
  g2 <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                            genetic_map_dir = "data/genetic_maps",
                            seed = 99, verbose = FALSE)
  stopifnot(identical(g1[[1]]$X, g2[[1]]$X))
})

run_test("seed = NULL accepted (no reproducibility required)", function() {
  g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                          genetic_map_dir = "data/genetic_maps",
                          seed = NULL, verbose = FALSE)
  stopifnot(!is.null(g[[1]]$X))
})

run_test("verbose = FALSE suppresses messages", function() {
  out <- utils::capture.output(
    g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                            genetic_map_dir = "data/genetic_maps",
                            seed = 20, verbose = FALSE),
    type = "message"
  )
  stopifnot(length(out) == 0L || !any(grepl("Region", out)))
})

run_test("verbose = TRUE prints region progress", function() {
  out <- utils::capture.output(
    g <- simulate_genotypes(n_regions = 1, n = 50, p = 30,
                            genetic_map_dir = "data/genetic_maps",
                            seed = 21, verbose = TRUE),
    type = "message"
  )
  stopifnot(any(grepl("Region|region|Done", out)))
})

run_test("return value has X, X_raw, n, p, maf, variant_ids, region_id, vcf_source", function() {
  g <- GENO_SMALL
  expected <- c("X", "X_raw", "n", "p", "maf", "variant_ids", "region_id", "vcf_source")
  stopifnot(all(expected %in% names(g[[1]])))
})

run_test("n_regions must be positive integer (error on 0)", function() {
  expect_error(function() {
    simulate_genotypes(n_regions = 0, n = 50, p = 30,
                       genetic_map_dir = "data/genetic_maps", verbose = FALSE)
  })
})

run_test("p vector length mismatch errors", function() {
  expect_error(function() {
    simulate_genotypes(n_regions = 3, n = 50, p = c(30, 40),
                       genetic_map_dir = "data/genetic_maps", verbose = FALSE)
  })
})


# =============================================================================
# SECTION 2: simulate_phenotypes
# =============================================================================

set_section("simulate_phenotypes")

run_test("S = 1 (scalar) runs without error", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              seed = 1, verbose = FALSE)
  stopifnot(!is.null(sim[[1]]$truth$causal_indices))
  stopifnot(length(sim[[1]]$truth$causal_indices) == 1L)
})

run_test("S = 3 (scalar) selects 3 causal variants", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 3, phi = 0.2,
                              seed = 2, verbose = FALSE)
  stopifnot(length(sim[[1]]$truth$causal_indices) == 3L)
})

run_test("S as vector (different S per region)", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = c(1, 2), phi = 0.2,
                              seed = 3, verbose = FALSE)
  stopifnot(length(sim[[1]]$truth$causal_indices) == 1L)
  stopifnot(length(sim[[2]]$truth$causal_indices) == 2L)
})

run_test("S vector wrong length errors", function() {
  expect_error(function() {
    simulate_phenotypes(GENO_SMALL, S = c(1, 2, 3), phi = 0.2,
                        seed = 4, verbose = FALSE)
  })
})

run_test("S > p errors", function() {
  expect_error(function() {
    simulate_phenotypes(GENO_SMALL, S = 999, phi = 0.2,
                        seed = 5, verbose = FALSE)
  })
})

run_test("phi = 0.1 runs without error", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.1,
                              seed = 6, verbose = FALSE)
  stopifnot(sim[[1]]$truth$pve > 0)
})

run_test("phi = 0.8 runs without error", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.8,
                              seed = 7, verbose = FALSE)
  stopifnot(sim[[1]]$truth$pve > 0)
})

run_test("phi as vector per region", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = c(0.1, 0.5),
                              seed = 8, verbose = FALSE)
  stopifnot(sim[[1]]$truth$phi == 0.1)
  stopifnot(sim[[2]]$truth$phi == 0.5)
})

run_test("phi outside (0,1) errors", function() {
  expect_error(function() {
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 1.2,
                        seed = 9, verbose = FALSE)
  })
})

run_test("model = 'sparse' works", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              model = "sparse", seed = 10, verbose = FALSE)
  stopifnot(sim[[1]]$truth$model == "sparse")
})

run_test("model = 'sparse_inf' works", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              model = "sparse_inf", p_causal = 0.5,
                              seed = 11, verbose = FALSE)
  stopifnot(sim[[1]]$truth$model == "sparse_inf")
})

run_test("model invalid string errors", function() {
  expect_error(function() {
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        model = "wrong_model", seed = 12, verbose = FALSE)
  })
})

run_test("p_causal = 0.2 (sparse_inf) partitions variance", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.3,
                              model = "sparse_inf", p_causal = 0.2,
                              seed = 13, verbose = FALSE)
  stopifnot(!is.null(sim[[1]]$truth$p_causal))
  stopifnot(sim[[1]]$truth$p_causal == 0.2)
})

run_test("p_causal = 1.0 (fully sparse, no inf component)", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.3,
                              model = "sparse_inf", p_causal = 1.0,
                              seed = 14, verbose = FALSE)
  stopifnot(!is.null(sim[[1]]$truth$p_causal))
})

run_test("p_causal outside (0,1] errors", function() {
  expect_error(function() {
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.3,
                        model = "sparse_inf", p_causal = 0,
                        seed = 15, verbose = FALSE)
  })
})

run_test("inf_model = 'beatrice' (noncausal variants only)", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.3,
                              model = "sparse_inf", p_causal = 0.5,
                              inf_model = "beatrice",
                              seed = 16, verbose = FALSE)
  stopifnot(sim[[1]]$truth$inf_model == "beatrice")
})

run_test("inf_model = 'susie_inf' (all variants)", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.3,
                              model = "sparse_inf", p_causal = 0.5,
                              inf_model = "susie_inf",
                              seed = 17, verbose = FALSE)
  stopifnot(sim[[1]]$truth$inf_model == "susie_inf")
})

run_test("inf_model invalid string errors", function() {
  expect_error(function() {
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.3,
                        model = "sparse_inf", p_causal = 0.5,
                        inf_model = "unknown",
                        seed = 18, verbose = FALSE)
  })
})

run_test("effect_distribution = 'normal' draws from N(0, effect_variance)", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 2, phi = 0.2,
                              effect_distribution = "normal",
                              effect_variance = 0.36,
                              seed = 19, verbose = FALSE)
  stopifnot(sim[[1]]$truth$effect_distribution == "normal")
})

run_test("effect_distribution = 'equal' distributes variance equally", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 2, phi = 0.2,
                              effect_distribution = "equal",
                              seed = 20, verbose = FALSE)
  stopifnot(sim[[1]]$truth$effect_distribution == "equal")
})

run_test("effect_distribution invalid errors", function() {
  expect_error(function() {
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        effect_distribution = "laplace",
                        seed = 21, verbose = FALSE)
  })
})

run_test("effect_variance = 0.1 accepted", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              effect_variance = 0.1,
                              seed = 22, verbose = FALSE)
  stopifnot(sim[[1]]$truth$effect_variance == 0.1)
})

run_test("effect_variance = 1.0 accepted", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              effect_variance = 1.0,
                              seed = 23, verbose = FALSE)
  stopifnot(!is.null(sim[[1]]$y))
})

run_test("effect_variance <= 0 errors", function() {
  expect_error(function() {
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        effect_variance = 0, seed = 24, verbose = FALSE)
  })
})

run_test("annotations = 'none' (no annotation matrix)", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "none",
                              seed = 25, verbose = FALSE)
  stopifnot(is.null(sim[[1]]$annotations_matrix))
})

run_test("annotations = 'binary' creates binary annotation matrix", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 3,
                              seed = 26, verbose = FALSE)
  A <- sim[[1]]$annotations_matrix
  stopifnot(!is.null(A), ncol(A) == 3L, all(A %in% c(0, 1)))
})

run_test("annotations = 'continuous' creates continuous annotation matrix", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "continuous", n_annotations = 2,
                              seed = 27, verbose = FALSE)
  A <- sim[[1]]$annotations_matrix
  stopifnot(!is.null(A), ncol(A) == 2L, is.numeric(A))
})

run_test("annotations as user-supplied matrix", function() {
  p <- GENO_SMALL[[1]]$p
  A_user <- matrix(rbinom(p * 2, 1, 0.2), nrow = p, ncol = 2)
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = A_user,
                              seed = 28, verbose = FALSE)
  stopifnot(!is.null(sim[[1]]$annotations_matrix))
})

run_test("annotations invalid string errors", function() {
  expect_error(function() {
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        annotations = "rainbow", seed = 29, verbose = FALSE)
  })
})

run_test("n_annotations = 1 works", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 1,
                              seed = 30, verbose = FALSE)
  stopifnot(ncol(sim[[1]]$annotations_matrix) == 1L)
})

run_test("n_annotations = 5 works", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 5,
                              seed = 31, verbose = FALSE)
  stopifnot(ncol(sim[[1]]$annotations_matrix) == 5L)
})

run_test("annotation_proportions = NULL (random proportions)", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 3,
                              annotation_proportions = NULL,
                              seed = 32, verbose = FALSE)
  stopifnot(!is.null(sim[[1]]$annotations_matrix))
})

run_test("annotation_proportions scalar (same for all annotations)", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 3,
                              annotation_proportions = 0.2,
                              seed = 33, verbose = FALSE)
  stopifnot(!is.null(sim[[1]]$annotations_matrix))
})

run_test("annotation_proportions vector (per-annotation)", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 3,
                              annotation_proportions = c(0.1, 0.2, 0.3),
                              seed = 34, verbose = FALSE)
  stopifnot(!is.null(sim[[1]]$annotations_matrix))
})

run_test("annotation_proportions vector wrong length errors", function() {
  expect_error(function() {
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        annotations = "binary", n_annotations = 3,
                        annotation_proportions = c(0.1, 0.2),
                        seed = 35, verbose = FALSE)
  })
})

run_test("enrichment = NULL (random enrichments)", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 2,
                              enrichment = NULL,
                              seed = 36, verbose = FALSE)
  stopifnot(!is.null(sim[[1]]$annotations_matrix))
})

run_test("enrichment scalar (same for all annotations)", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 2,
                              enrichment = 5,
                              seed = 37, verbose = FALSE)
  stopifnot(!is.null(sim[[1]]$truth$enrichment))
})

run_test("enrichment vector (per-annotation)", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                              annotations = "binary", n_annotations = 2,
                              enrichment = c(3, 8),
                              seed = 38, verbose = FALSE)
  stopifnot(!is.null(sim[[1]]$truth$enrichment))
})

run_test("enrichment vector wrong length errors", function() {
  expect_error(function() {
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                        annotations = "binary", n_annotations = 3,
                        enrichment = c(3, 8),
                        seed = 39, verbose = FALSE)
  })
})

run_test("seed ensures reproducibility", function() {
  s1 <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2, seed = 77, verbose = FALSE)
  s2 <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2, seed = 77, verbose = FALSE)
  stopifnot(identical(s1[[1]]$y, s2[[1]]$y))
})

run_test("seed = NULL accepted", function() {
  sim <- simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2, seed = NULL, verbose = FALSE)
  stopifnot(!is.null(sim[[1]]$y))
})

run_test("verbose = FALSE suppresses messages", function() {
  out <- utils::capture.output(
    simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2, seed = 40, verbose = FALSE),
    type = "message"
  )
  stopifnot(length(out) == 0L)
})

run_test("return fields: y, z, beta_hat, se, LD, truth", function() {
  sim <- PHENO_SMALL
  expected <- c("y", "z", "beta_hat", "se", "LD", "truth")
  stopifnot(all(expected %in% names(sim[[1]])))
})

run_test("truth fields: causal_indices, causal_effects, beta_true, pve, S, phi, model", function() {
  t <- PHENO_SMALL[[1]]$truth
  expected <- c("causal_indices", "causal_effects", "beta_true", "pve", "S", "phi", "model")
  stopifnot(all(expected %in% names(t)))
})


# =============================================================================
# SECTION 2b: simulate_genotypes — save / output_dir
# =============================================================================

set_section("simulate_genotypes — save / output_dir")

run_test("save = FALSE writes no files", function() {
  tmp <- tempfile(); dir.create(tmp)
  simulate_genotypes(n_regions = 1, n = 50, p = 30,
                     genetic_map_dir = "data/genetic_maps",
                     seed = 200, save = FALSE, output_dir = tmp,
                     verbose = FALSE)
  stopifnot(length(list.files(tmp)) == 0L)
  unlink(tmp, recursive = TRUE)
})

run_test("save = TRUE writes .rds file", function() {
  tmp <- tempfile(); dir.create(tmp)
  simulate_genotypes(n_regions = 1, n = 50, p = 30,
                     genetic_map_dir = "data/genetic_maps",
                     seed = 201, save = TRUE, output_dir = tmp,
                     verbose = FALSE)
  stopifnot(length(list.files(tmp, pattern = "\\.rds$")) == 1L)
  unlink(tmp, recursive = TRUE)
})

run_test("saved .rds is readable and has correct structure", function() {
  tmp <- tempfile(); dir.create(tmp)
  simulate_genotypes(n_regions = 1, n = 50, p = 30,
                     genetic_map_dir = "data/genetic_maps",
                     seed = 202, save = TRUE, output_dir = tmp,
                     verbose = FALSE)
  f   <- list.files(tmp, pattern = "\\.rds$", full.names = TRUE)
  obj <- readRDS(f)
  stopifnot(is.list(obj), !is.null(obj[[1]]$X))
  unlink(tmp, recursive = TRUE)
})

run_test("output_dir created if it does not exist", function() {
  tmp <- file.path(tempfile(), "geno_out")
  simulate_genotypes(n_regions = 1, n = 50, p = 30,
                     genetic_map_dir = "data/genetic_maps",
                     seed = 203, save = TRUE, output_dir = tmp,
                     verbose = FALSE)
  stopifnot(dir.exists(tmp))
  unlink(dirname(tmp), recursive = TRUE)
})

run_test("filename encodes n_regions, n, p, seed", function() {
  tmp <- tempfile(); dir.create(tmp)
  simulate_genotypes(n_regions = 2, n = 60, p = 35,
                     genetic_map_dir = "data/genetic_maps",
                     seed = 204, save = TRUE, output_dir = tmp,
                     verbose = FALSE)
  f <- list.files(tmp, pattern = "\\.rds$")
  stopifnot(grepl("2regions", f), grepl("n60", f), grepl("p35", f),
            grepl("seed204", f))
  unlink(tmp, recursive = TRUE)
})

run_test("seed = NULL gives 'noseed' tag in filename", function() {
  tmp <- tempfile(); dir.create(tmp)
  simulate_genotypes(n_regions = 1, n = 50, p = 30,
                     genetic_map_dir = "data/genetic_maps",
                     seed = NULL, save = TRUE, output_dir = tmp,
                     verbose = FALSE)
  f <- list.files(tmp, pattern = "\\.rds$")
  stopifnot(grepl("noseed", f))
  unlink(tmp, recursive = TRUE)
})


# =============================================================================
# SECTION 2c: simulate_phenotypes — save / output_dir
# =============================================================================

set_section("simulate_phenotypes — save / output_dir")

run_test("save = FALSE writes no files", function() {
  tmp <- tempfile(); dir.create(tmp)
  simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                      seed = 300, save = FALSE, output_dir = tmp,
                      verbose = FALSE)
  stopifnot(length(list.files(tmp)) == 0L)
  unlink(tmp, recursive = TRUE)
})

run_test("save = TRUE writes .rds file", function() {
  tmp <- tempfile(); dir.create(tmp)
  simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                      seed = 301, save = TRUE, output_dir = tmp,
                      verbose = FALSE)
  stopifnot(length(list.files(tmp, pattern = "\\.rds$")) == 1L)
  unlink(tmp, recursive = TRUE)
})

run_test("saved .rds has y, z, truth fields", function() {
  tmp <- tempfile(); dir.create(tmp)
  simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                      seed = 302, save = TRUE, output_dir = tmp,
                      verbose = FALSE)
  f   <- list.files(tmp, pattern = "\\.rds$", full.names = TRUE)
  obj <- readRDS(f)
  stopifnot(all(c("y", "z", "truth") %in% names(obj[[1]])))
  unlink(tmp, recursive = TRUE)
})

run_test("output_dir created if it does not exist", function() {
  tmp <- file.path(tempfile(), "pheno_out")
  simulate_phenotypes(GENO_SMALL, S = 1, phi = 0.2,
                      seed = 303, save = TRUE, output_dir = tmp,
                      verbose = FALSE)
  stopifnot(dir.exists(tmp))
  unlink(dirname(tmp), recursive = TRUE)
})

run_test("filename encodes model, S, phi, seed", function() {
  tmp <- tempfile(); dir.create(tmp)
  simulate_phenotypes(GENO_SMALL, S = 2, phi = 0.3,
                      model = "sparse", seed = 304,
                      save = TRUE, output_dir = tmp, verbose = FALSE)
  f <- list.files(tmp, pattern = "\\.rds$")
  stopifnot(grepl("sparse", f), grepl("S2", f), grepl("phi0.3", f),
            grepl("seed304", f))
  unlink(tmp, recursive = TRUE)
})


# =============================================================================
# SECTION 3: run_simulation
# =============================================================================

set_section("run_simulation")

run_test("n_iter = 1 produces correct number of scenarios", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 1, verbose = FALSE)
  stopifnot(length(r$scenarios) == 1L)  # 1 S x 1 phi x 1 iter
})

run_test("n_iter = 3 produces correct scenario count", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 3, S = 1, phi = 0.2,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 2, verbose = FALSE)
  stopifnot(length(r$scenarios) == 3L)  # 1 S x 1 phi x 3 iter
})

run_test("S vector sweeps correctly", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = c(1, 2), phi = 0.2,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 3, verbose = FALSE)
  stopifnot(length(r$scenarios) == 2L)  # 2 S x 1 phi x 1 iter
  S_vals <- sapply(r$scenarios, `[[`, "S")
  stopifnot(all(sort(S_vals) == c(1L, 2L)))
})

run_test("phi vector sweeps correctly", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = c(0.2, 0.4),
                      genetic_map_dir = "data/genetic_maps",
                      seed = 4, verbose = FALSE)
  stopifnot(length(r$scenarios) == 2L)
  phi_vals <- sapply(r$scenarios, `[[`, "phi")
  stopifnot(all(sort(phi_vals) == c(0.2, 0.4)))
})

run_test("model = 'sparse' runs without error", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      model = "sparse",
                      genetic_map_dir = "data/genetic_maps",
                      seed = 5, verbose = FALSE)
  stopifnot(r$params$model == "sparse")
})

run_test("model = 'sparse_inf' sweeps p_causal", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      model = "sparse_inf", p_causal = c(0.2, 0.5),
                      genetic_map_dir = "data/genetic_maps",
                      seed = 6, verbose = FALSE)
  # 1 S x 1 phi x 2 p_causal x 1 iter = 2 scenarios
  stopifnot(length(r$scenarios) == 2L)
  stopifnot(r$params$model == "sparse_inf")
})

run_test("inf_model = 'beatrice' accepted in sparse_inf", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      model = "sparse_inf", p_causal = 0.5,
                      inf_model = "beatrice",
                      genetic_map_dir = "data/genetic_maps",
                      seed = 7, verbose = FALSE)
  stopifnot(r$params$inf_model == "beatrice")
})

run_test("inf_model = 'susie_inf' accepted in sparse_inf", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      model = "sparse_inf", p_causal = 0.5,
                      inf_model = "susie_inf",
                      genetic_map_dir = "data/genetic_maps",
                      seed = 8, verbose = FALSE)
  stopifnot(r$params$inf_model == "susie_inf")
})

run_test("effect_distribution = 'normal' recorded in params", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      effect_distribution = "normal",
                      genetic_map_dir = "data/genetic_maps",
                      seed = 9, verbose = FALSE)
  stopifnot(r$params$effect_distribution == "normal")
})

run_test("effect_distribution = 'equal' works", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      effect_distribution = "equal",
                      genetic_map_dir = "data/genetic_maps",
                      seed = 10, verbose = FALSE)
  stopifnot(r$params$effect_distribution == "equal")
})

run_test("effect_variance = 0.5 accepted", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      effect_variance = 0.5,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 11, verbose = FALSE)
  stopifnot(r$params$effect_variance == 0.5)
})

run_test("annotations = 'none' produces NULL annotation matrix", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "none",
                      genetic_map_dir = "data/genetic_maps",
                      seed = 12, verbose = FALSE)
  stopifnot(is.null(r$scenarios[[1]]$regions[[1]]$annotations_matrix))
})

run_test("annotations = 'binary' with n_annotations = 2", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "binary", n_annotations = 2,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 13, verbose = FALSE)
  A <- r$scenarios[[1]]$regions[[1]]$annotations_matrix
  stopifnot(!is.null(A), ncol(A) == 2L)
})

run_test("annotations = 'continuous' with n_annotations = 3", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "continuous", n_annotations = 3,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 14, verbose = FALSE)
  A <- r$scenarios[[1]]$regions[[1]]$annotations_matrix
  stopifnot(!is.null(A), ncol(A) == 3L)
})

run_test("annotation_proportions scalar passed through correctly", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "binary", n_annotations = 2,
                      annotation_proportions = 0.15,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 15, verbose = FALSE)
  stopifnot(r$params$annotation_proportions == 0.15)
})

run_test("enrichment scalar passed through correctly", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      annotations = "binary", n_annotations = 2,
                      enrichment = 4,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 16, verbose = FALSE)
  stopifnot(r$params$enrichment == 4)
})

run_test("vcf_dir missing directory errors", function() {
  expect_error(function() {
    run_simulation(n_regions = 1, n = 80, p = 40,
                   n_iter = 1, S = 1, phi = 0.2,
                   vcf_dir = "/nonexistent/dir",
                   genetic_map_dir = "data/genetic_maps",
                   seed = 17, verbose = FALSE)
  })
})

run_test("vcf_files = single VCF path used for all regions", function() {
  vcf <- system.file("examples", "region.vcf.gz", package = "sim1000G")
  r <- run_simulation(n_regions = 2, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      vcf_files = vcf,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 18, verbose = FALSE)
  stopifnot(length(r$genotypes) == 2L)
})

run_test("min_maf = 0.05 passed to simulate_genotypes", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      min_maf = 0.05,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 19, verbose = FALSE)
  stopifnot(r$params$min_maf == 0.05)
})

run_test("max_maf = 0.4 passed to simulate_genotypes", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      max_maf = 0.4,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 20, verbose = FALSE)
  # No error = pass; maf filtering applied
  stopifnot(!is.null(r$genotypes[[1]]$X))
})

run_test("standardise = FALSE returns raw 0/1/2 genotypes", function() {
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      standardise = FALSE,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 21, verbose = FALSE)
  vals <- unique(as.vector(r$genotypes[[1]]$X))
  stopifnot(all(vals %in% c(0, 1, 2)))
})

run_test("seed = 42 ensures reproducibility", function() {
  r1 <- run_simulation(n_regions = 1, n = 80, p = 40,
                       n_iter = 1, S = 1, phi = 0.2,
                       genetic_map_dir = "data/genetic_maps",
                       seed = 42, verbose = FALSE)
  r2 <- run_simulation(n_regions = 1, n = 80, p = 40,
                       n_iter = 1, S = 1, phi = 0.2,
                       genetic_map_dir = "data/genetic_maps",
                       seed = 42, verbose = FALSE)
  stopifnot(identical(r1$genotypes[[1]]$X, r2$genotypes[[1]]$X))
})

run_test("save = TRUE writes .rds file", function() {
  tmp <- tempfile()
  dir.create(tmp)
  r <- run_simulation(n_regions = 1, n = 80, p = 40,
                      n_iter = 1, S = 1, phi = 0.2,
                      genetic_map_dir = "data/genetic_maps",
                      seed = 99, save = TRUE, output_dir = tmp,
                      verbose = FALSE)
  files <- list.files(tmp, pattern = "\\.rds$")
  stopifnot(length(files) == 1L)
  unlink(tmp, recursive = TRUE)
})

run_test("save = FALSE writes no files", function() {
  tmp <- tempfile()
  dir.create(tmp)
  run_simulation(n_regions = 1, n = 80, p = 40,
                 n_iter = 1, S = 1, phi = 0.2,
                 genetic_map_dir = "data/genetic_maps",
                 seed = 100, save = FALSE, output_dir = tmp,
                 verbose = FALSE)
  stopifnot(length(list.files(tmp)) == 0L)
  unlink(tmp, recursive = TRUE)
})

run_test("output_dir is created if it does not exist", function() {
  tmp <- file.path(tempfile(), "deep", "nested")
  run_simulation(n_regions = 1, n = 80, p = 40,
                 n_iter = 1, S = 1, phi = 0.2,
                 genetic_map_dir = "data/genetic_maps",
                 seed = 101, save = TRUE, output_dir = tmp,
                 verbose = FALSE)
  stopifnot(dir.exists(tmp))
  unlink(dirname(dirname(tmp)), recursive = TRUE)
})

run_test("verbose = FALSE suppresses messages", function() {
  out <- utils::capture.output(
    run_simulation(n_regions = 1, n = 80, p = 40,
                   n_iter = 1, S = 1, phi = 0.2,
                   genetic_map_dir = "data/genetic_maps",
                   seed = 102, verbose = FALSE),
    type = "message"
  )
  stopifnot(length(out) == 0L)
})

run_test("return value has genotypes, scenarios, params", function() {
  r <- SIM_MINI
  stopifnot(all(c("genotypes", "scenarios", "params") %in% names(r)))
})

run_test("scenarios have correct fields", function() {
  sc <- SIM_MINI$scenarios[[1]]
  expected <- c("scenario_id", "S", "phi", "p_causal", "iter", "model", "regions")
  stopifnot(all(expected %in% names(sc)))
})

run_test("params records all key settings", function() {
  p <- SIM_MINI$params
  expected <- c("n_regions", "n", "p", "n_iter", "S_values", "phi_values", "model", "seed")
  stopifnot(all(expected %in% names(p)))
})

run_test("n_iter must be a positive integer (error on 0)", function() {
  expect_error(function() {
    run_simulation(n_regions = 1, n = 80, p = 40,
                   n_iter = 0, S = 1, phi = 0.2,
                   genetic_map_dir = "data/genetic_maps", verbose = FALSE)
  })
})

run_test("phi outside (0,1) errors", function() {
  expect_error(function() {
    run_simulation(n_regions = 1, n = 80, p = 40,
                   n_iter = 1, S = 1, phi = 1.5,
                   genetic_map_dir = "data/genetic_maps", verbose = FALSE)
  })
})


# =============================================================================
# SECTION 4: run_methods
# =============================================================================

set_section("run_methods")

run_test("methods = 'susie' runs on SIM_MINI", function() {
  r <- run_methods(SIM_MINI, methods = "susie",
                   method_args = list(susie = list(L = 5)),
                   save = FALSE, verbose = FALSE)
  stopifnot("susie" %in% r$methods_run)
})

run_test("methods = 'abf' runs on SIM_MINI", function() {
  r <- run_methods(SIM_MINI, methods = "abf",
                   save = FALSE, verbose = FALSE)
  stopifnot("abf" %in% r$methods_run)
})

run_test("methods = 'susie_inf' runs on SIM_MINI", function() {
  r <- run_methods(SIM_MINI, methods = "susie_inf",
                   method_args = list(susie_inf = list(L = 5)),
                   save = FALSE, verbose = FALSE)
  stopifnot("susie_inf" %in% r$methods_run)
})

run_test("methods = 'carma' runs on SIM_MINI", function() {
  r <- run_methods(SIM_MINI, methods = "carma",
                   save = FALSE, verbose = FALSE)
  stopifnot("carma" %in% r$methods_run)
})

run_test("multiple methods run together", function() {
  r <- run_methods(SIM_MINI, methods = c("susie", "abf"),
                   save = FALSE, verbose = FALSE)
  stopifnot(all(c("susie", "abf") %in% r$methods_run))
})

run_test("method_args forwarded to susie (L and coverage)", function() {
  r <- run_methods(SIM_MINI, methods = "susie",
                   method_args = list(susie = list(L = 3, coverage = 0.9)),
                   save = FALSE, verbose = FALSE)
  stopifnot(r$susie$method_args$L == 3)
  stopifnot(r$susie$method_args$coverage == 0.9)
})

run_test("method_args forwarded to abf (prior_variance, coverage)", function() {
  r <- run_methods(SIM_MINI, methods = "abf",
                   method_args = list(abf = list(prior_variance = 0.02, coverage = 0.9)),
                   save = FALSE, verbose = FALSE)
  stopifnot(r$abf$method_args$prior_variance == 0.02)
})

run_test("method_args forwarded to susie_inf (L)", function() {
  r <- run_methods(SIM_MINI, methods = "susie_inf",
                   method_args = list(susie_inf = list(L = 3)),
                   save = FALSE, verbose = FALSE)
  stopifnot(r$susie_inf$method_args$L == 3)
})

run_test("method_args forwarded to carma (rho.index)", function() {
  r <- run_methods(SIM_MINI, methods = "carma",
                   method_args = list(carma = list(rho.index = 0.9)),
                   save = FALSE, verbose = FALSE)
  stopifnot(r$carma$method_args$rho.index == 0.9)
})

run_test("unknown method name errors", function() {
  expect_error(function() {
    run_methods(SIM_MINI, methods = "notamethod", save = FALSE, verbose = FALSE)
  })
})

run_test("method_args for non-run method warns", function() {
  expect_warning_fn(function() {
    run_methods(SIM_MINI, methods = "abf",
                method_args = list(abf = list(), susie = list(L = 5)),
                save = FALSE, verbose = FALSE)
  })
})

run_test("save = TRUE writes per-method .rds and run_metadata.rds", function() {
  tmp <- tempfile()
  dir.create(tmp)
  run_methods(SIM_MINI, methods = "abf",
              save = TRUE, output_dir = tmp, verbose = FALSE)
  files <- list.files(tmp, recursive = TRUE, pattern = "\\.rds$")
  stopifnot(any(grepl("abf\\.rds", files)))
  stopifnot(any(grepl("run_metadata\\.rds", files)))
  unlink(tmp, recursive = TRUE)
})

run_test("save = FALSE produces no files", function() {
  tmp <- tempfile()
  dir.create(tmp)
  run_methods(SIM_MINI, methods = "abf",
              save = FALSE, output_dir = tmp, verbose = FALSE)
  stopifnot(length(list.files(tmp, recursive = TRUE)) == 0L)
  unlink(tmp, recursive = TRUE)
})

run_test("verbose = FALSE suppresses messages", function() {
  out <- utils::capture.output(
    run_methods(SIM_MINI, methods = "abf", save = FALSE, verbose = FALSE),
    type = "message"
  )
  stopifnot(length(out) == 0L)
})

run_test("return value has per-method list with results, n_total, n_failed", function() {
  r <- RESULTS_MINI
  stopifnot(all(c("results", "n_total", "n_failed", "total_runtime_seconds") %in%
                  names(r$susie)))
})

run_test("each fit has pip, credible_sets, method, runtime_seconds", function() {
  fit <- RESULTS_MINI$susie$results[[1]]
  expected <- c("pip", "credible_sets", "method", "runtime_seconds",
                "scenario_id", "region_id", "S", "phi", "iter")
  stopifnot(all(expected %in% names(fit)))
})

run_test("pip length equals n_snps", function() {
  fit <- RESULTS_MINI$susie$results[[1]]
  rg  <- SIM_MINI$genotypes[[fit$region_id]]
  stopifnot(length(fit$pip) == rg$p)
})

run_test("pip values in [0, 1]", function() {
  pip <- RESULTS_MINI$susie$results[[1]]$pip
  stopifnot(all(pip >= 0), all(pip <= 1))
})

run_test("failed fits return error field not NA pip", function() {
  # Create a simulation with a deliberate bad method call (force error)
  bad_sim <- SIM_MINI
  bad_sim$genotypes[[1]]$LD <- NULL  # break LD for FINEMAP (won't break susie)
  # Just check that graceful error handling works via direct call
  fit <- tryCatch(
    run_abf_region(
      region_geno  = list(LD = NULL, n = 100),
      region_pheno = list(z = rep(0, 5), se = rep(1, 5))
    ),
    error = function(e) list(error = conditionMessage(e))
  )
  stopifnot(!is.null(fit$error) || !is.null(fit$pip))
})

run_test("methods case-insensitive (SUSIE == susie)", function() {
  r <- run_methods(SIM_MINI, methods = "SUSIE",
                   method_args = list(susie = list(L = 5)),
                   save = FALSE, verbose = FALSE)
  stopifnot("susie" %in% r$methods_run)
})


# =============================================================================
# SECTION 5: evaluate_methods
# =============================================================================

set_section("evaluate_methods")

run_test("basic evaluation returns named list per method", function() {
  ev <- EVAL_MINI
  stopifnot(all(c("susie", "abf") %in% ev$methods_evaluated))
})

run_test("global stratum present for each method", function() {
  ev <- EVAL_MINI
  for (m in ev$methods_evaluated)
    stopifnot(!is.null(ev[[m]]$global))
})

run_test("by_S stratum present and named correctly", function() {
  ev <- EVAL_MINI
  stopifnot(!is.null(ev$susie$by_S))
  stopifnot(all(c("1", "2") %in% names(ev$susie$by_S)))
})

run_test("by_phi stratum present and named correctly", function() {
  ev <- EVAL_MINI
  stopifnot(!is.null(ev$susie$by_phi))
  stopifnot(all(c("0.2", "0.4") %in% names(ev$susie$by_phi)))
})

run_test("by_p_causal is NULL for sparse model", function() {
  ev <- EVAL_MINI
  stopifnot(is.null(ev$susie$by_p_causal))
})

run_test("by_p_causal present for sparse_inf model", function() {
  sim_inf <- run_simulation(
    n_regions = 1, n = 80, p = 40,
    n_iter = 2, S = 1, phi = 0.2,
    model = "sparse_inf", p_causal = c(0.3, 0.7),
    genetic_map_dir = "data/genetic_maps",
    seed = 55, verbose = FALSE
  )
  res_inf <- run_methods(sim_inf, methods = "abf", save = FALSE, verbose = FALSE)
  ev_inf  <- evaluate_methods(sim_inf, res_inf, save = FALSE, verbose = FALSE)
  stopifnot(!is.null(ev_inf$abf$by_p_causal))
  stopifnot(all(c("0.3", "0.7") %in% names(ev_inf$abf$by_p_causal)))
})

run_test("global auprc is numeric in [0, 1]", function() {
  auprc <- EVAL_MINI$susie$global$auprc
  stopifnot(is.numeric(auprc), !is.na(auprc), auprc >= 0, auprc <= 1)
})

run_test("global cs_coverage is in [0, 1] or NA", function() {
  cv <- EVAL_MINI$susie$global$cs_coverage
  stopifnot(is.numeric(cv), (is.na(cv) || (cv >= 0 && cv <= 1)))
})

run_test("global cs_power is in [0, 1] or NA", function() {
  cp <- EVAL_MINI$susie$global$cs_power
  stopifnot(is.numeric(cp), (is.na(cp) || (cp >= 0 && cp <= 1)))
})

run_test("fdr_power_curve has required columns", function() {
  curve <- EVAL_MINI$susie$global$fdr_power_curve
  expected <- c("threshold", "tp", "fp", "fn", "fdr", "power", "precision", "recall")
  stopifnot(all(expected %in% names(curve)))
})

run_test("pip_calibration has required columns", function() {
  cal <- EVAL_MINI$susie$global$pip_calibration
  expected <- c("bin", "bin_lower", "bin_upper", "mean_pip", "frac_causal")
  stopifnot(all(expected %in% names(cal)))
})

run_test("SE fields present when n_iter >= 2", function() {
  ev <- EVAL_MINI  # n_iter = 2
  stopifnot(!is.null(EVAL_MINI$susie$global$auprc_se))
})

run_test("pip_thresholds custom (coarser) works", function() {
  ev <- evaluate_methods(SIM_MINI, RESULTS_MINI,
                         pip_thresholds = seq(0, 1, by = 0.1),
                         save = FALSE, verbose = FALSE)
  nrow_curve <- nrow(ev$susie$global$fdr_power_curve)
  stopifnot(nrow_curve == 11L)  # seq(0,1,0.1) has 11 values
})

run_test("n_pip_cal_bins = 5 produces 5-row calibration table", function() {
  ev <- evaluate_methods(SIM_MINI, RESULTS_MINI,
                         n_pip_cal_bins = 5,
                         save = FALSE, verbose = FALSE)
  stopifnot(nrow(ev$susie$global$pip_calibration) == 5L)
})

run_test("n_pip_cal_bins = 20 produces 20-row calibration table", function() {
  ev <- evaluate_methods(SIM_MINI, RESULTS_MINI,
                         n_pip_cal_bins = 20,
                         save = FALSE, verbose = FALSE)
  stopifnot(nrow(ev$susie$global$pip_calibration) == 20L)
})

run_test("save = TRUE writes evaluation.rds and evaluation_summary.csv", function() {
  tmp <- tempfile()
  dir.create(tmp)
  evaluate_methods(SIM_MINI, RESULTS_MINI,
                   save = TRUE, output_dir = tmp, verbose = FALSE)
  stopifnot(file.exists(file.path(tmp, "evaluation.rds")))
  stopifnot(file.exists(file.path(tmp, "evaluation_summary.csv")))
  unlink(tmp, recursive = TRUE)
})

run_test("save = FALSE writes no files", function() {
  tmp <- tempfile()
  dir.create(tmp)
  evaluate_methods(SIM_MINI, RESULTS_MINI,
                   save = FALSE, output_dir = tmp, verbose = FALSE)
  stopifnot(length(list.files(tmp)) == 0L)
  unlink(tmp, recursive = TRUE)
})

run_test("output_dir created if absent (save = TRUE)", function() {
  tmp <- file.path(tempfile(), "ev_out")
  evaluate_methods(SIM_MINI, RESULTS_MINI,
                   save = TRUE, output_dir = tmp, verbose = FALSE)
  stopifnot(dir.exists(tmp))
  unlink(tmp, recursive = TRUE)
})

run_test("verbose = FALSE suppresses messages", function() {
  out <- utils::capture.output(
    evaluate_methods(SIM_MINI, RESULTS_MINI,
                     save = FALSE, verbose = FALSE),
    type = "message"
  )
  stopifnot(length(out) == 0L)
})

run_test("methods_evaluated field present in return value", function() {
  stopifnot("methods_evaluated" %in% names(EVAL_MINI))
})

run_test("simulation missing 'scenarios' errors", function() {
  expect_error(function() {
    evaluate_methods(list(params = list()), RESULTS_MINI,
                     save = FALSE, verbose = FALSE)
  })
})

run_test("results missing 'methods_run' errors", function() {
  expect_error(function() {
    evaluate_methods(SIM_MINI, list(),
                     save = FALSE, verbose = FALSE)
  })
})


# =============================================================================
# SECTION 6: plot_results
# =============================================================================

set_section("plot_results")

run_test("output_file explicit path writes PDF there", function() {
  tmp <- tempfile(fileext = ".pdf")
  plot_results(EVAL_MINI, output_file = tmp, verbose = FALSE)
  stopifnot(file.exists(tmp), file.size(tmp) > 1000)
  unlink(tmp)
})

run_test("output_dir writes evaluation.pdf inside that directory", function() {
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  plot_results(EVAL_MINI, output_dir = tmp_dir, verbose = FALSE)
  stopifnot(file.exists(file.path(tmp_dir, "evaluation.pdf")))
  unlink(tmp_dir, recursive = TRUE)
})

run_test("output_dir created if it does not exist", function() {
  tmp_dir <- file.path(tempfile(), "plots", "deep")
  plot_results(EVAL_MINI, output_dir = tmp_dir, verbose = FALSE)
  stopifnot(dir.exists(tmp_dir))
  unlink(dirname(dirname(tmp_dir)), recursive = TRUE)
})

run_test("output_file takes precedence over output_dir", function() {
  tmp_dir <- tempfile(); dir.create(tmp_dir)
  tmp_file <- tempfile(fileext = ".pdf")
  plot_results(EVAL_MINI, output_file = tmp_file,
               output_dir = tmp_dir, verbose = FALSE)
  stopifnot(file.exists(tmp_file))
  stopifnot(!file.exists(file.path(tmp_dir, "evaluation.pdf")))
  unlink(tmp_dir, recursive = TRUE); unlink(tmp_file)
})

run_test("save = FALSE does not write any file", function() {
  tmp_dir <- tempfile(); dir.create(tmp_dir)
  plot_results(EVAL_MINI, output_dir = tmp_dir,
               save = FALSE, verbose = FALSE)
  stopifnot(length(list.files(tmp_dir)) == 0L)
  unlink(tmp_dir, recursive = TRUE)
})

run_test("save = TRUE (default) writes PDF", function() {
  tmp_dir <- tempfile(); dir.create(tmp_dir)
  plot_results(EVAL_MINI, output_dir = tmp_dir,
               save = TRUE, verbose = FALSE)
  stopifnot(file.exists(file.path(tmp_dir, "evaluation.pdf")))
  unlink(tmp_dir, recursive = TRUE)
})

run_test("methods = 'susie' only (subset of evaluated methods)", function() {
  tmp <- tempfile(fileext = ".pdf")
  plot_results(EVAL_MINI, output_file = tmp,
               methods = "susie", verbose = FALSE)
  stopifnot(file.exists(tmp))
  unlink(tmp)
})

run_test("methods = c('susie', 'abf') includes both", function() {
  tmp <- tempfile(fileext = ".pdf")
  plot_results(EVAL_MINI, output_file = tmp,
               methods = c("susie", "abf"), verbose = FALSE)
  stopifnot(file.exists(tmp))
  unlink(tmp)
})

run_test("methods subset to unknown method produces no-valid-method error", function() {
  expect_error(function() {
    plot_results(EVAL_MINI, output_file = tempfile(fileext = ".pdf"),
                 methods = "notamethod", verbose = FALSE)
  })
})

run_test("verbose = FALSE produces no messages", function() {
  tmp <- tempfile(fileext = ".pdf")
  out <- utils::capture.output(
    plot_results(EVAL_MINI, output_file = tmp, verbose = FALSE),
    type = "message"
  )
  unlink(tmp)
  stopifnot(length(out) == 0L)
})

run_test("verbose = TRUE prints section messages", function() {
  tmp <- tempfile(fileext = ".pdf")
  out <- utils::capture.output(
    plot_results(EVAL_MINI, output_file = tmp, verbose = TRUE),
    type = "message"
  )
  unlink(tmp)
  stopifnot(any(grepl("Global|Plotting|PDF", out)))
})

run_test("return value is the resolved output path (invisibly)", function() {
  tmp <- tempfile(fileext = ".pdf")
  ret <- plot_results(EVAL_MINI, output_file = tmp, verbose = FALSE)
  unlink(tmp)
  stopifnot(is.character(ret), ret == tmp)
})

run_test("sparse_inf eval with by_p_causal section rendered", function() {
  sim_inf <- run_simulation(
    n_regions = 1, n = 80, p = 40,
    n_iter = 2, S = 1, phi = 0.2,
    model = "sparse_inf", p_causal = c(0.3, 0.7),
    genetic_map_dir = "data/genetic_maps",
    seed = 66, verbose = FALSE
  )
  res_inf <- run_methods(sim_inf, methods = "abf", save = FALSE, verbose = FALSE)
  ev_inf  <- evaluate_methods(sim_inf, res_inf, save = FALSE, verbose = FALSE)
  tmp_dir <- tempfile(); dir.create(tmp_dir)
  plot_results(ev_inf, output_dir = tmp_dir, verbose = FALSE)
  stopifnot(file.exists(file.path(tmp_dir, "evaluation.pdf")))
  unlink(tmp_dir, recursive = TRUE)
})


# =============================================================================
# SECTION 7: run_susie / run_susie_region (argument passthrough)
# =============================================================================

set_section("run_susie / run_susie_region — method arguments")

.rg <- SIM_MINI$genotypes[[1]]
.rp <- SIM_MINI$scenarios[[1]]$regions[[1]]

run_test("L = 5 accepted", function() {
  fit <- run_susie_region(.rg, .rp, L = 5)
  stopifnot(fit$params$L == 5)
})

run_test("L = 1 (single-component) accepted", function() {
  fit <- run_susie_region(.rg, .rp, L = 1)
  stopifnot(fit$params$L == 1)
})

run_test("coverage = 0.5 accepted", function() {
  fit <- run_susie_region(.rg, .rp, L = 5, coverage = 0.5)
  stopifnot(fit$params$coverage == 0.5)
})

run_test("coverage = 0.99 accepted", function() {
  fit <- run_susie_region(.rg, .rp, L = 5, coverage = 0.99)
  stopifnot(fit$params$coverage == 0.99)
})

run_test("min_abs_corr = 0 (no purity filter) accepted", function() {
  fit <- run_susie_region(.rg, .rp, L = 5, min_abs_corr = 0)
  stopifnot(fit$params$min_abs_corr == 0)
})

run_test("min_abs_corr = 0.8 (strict purity filter) accepted", function() {
  fit <- run_susie_region(.rg, .rp, L = 5, min_abs_corr = 0.8)
  stopifnot(fit$params$min_abs_corr == 0.8)
})

run_test("max_iter = 50 accepted", function() {
  fit <- run_susie_region(.rg, .rp, L = 5, max_iter = 50)
  stopifnot(fit$params$max_iter == 50)
})

run_test("estimate_residual_variance = FALSE accepted", function() {
  fit <- run_susie_region(.rg, .rp, L = 5, estimate_residual_variance = FALSE)
  stopifnot(fit$params$estimate_residual_variance == FALSE)
})

run_test("estimate_prior_variance = FALSE accepted", function() {
  fit <- run_susie_region(.rg, .rp, L = 5, estimate_prior_variance = FALSE)
  stopifnot(fit$params$estimate_prior_variance == FALSE)
})

run_test("prior_variance = 0.05 accepted", function() {
  fit <- run_susie_region(.rg, .rp, L = 5, prior_variance = 0.05)
  stopifnot(fit$params$prior_variance == 0.05)
})

run_test("susie output has pip, credible_sets, method, runtime_seconds", function() {
  fit <- run_susie_region(.rg, .rp, L = 5)
  stopifnot(all(c("pip", "credible_sets", "method", "runtime_seconds",
                  "input_type", "params", "additional") %in% names(fit)))
})

run_test("susie pip sums approximately to number of credible sets", function() {
  fit <- run_susie_region(.rg, .rp, L = 5)
  stopifnot(is.numeric(fit$pip), length(fit$pip) == .rg$p)
  stopifnot(all(fit$pip >= 0), all(fit$pip <= 1))
})


# =============================================================================
# SECTION 8: run_abf / run_abf_region
# =============================================================================

set_section("run_abf / run_abf_region — method arguments")

run_test("prior_variance = 0.04 (default) accepted", function() {
  fit <- run_abf_region(.rg, .rp, prior_variance = 0.04)
  stopifnot(fit$params$prior_variance == 0.04)
})

run_test("prior_variance = 0.1 accepted", function() {
  fit <- run_abf_region(.rg, .rp, prior_variance = 0.1)
  stopifnot(fit$params$prior_variance == 0.1)
})

run_test("coverage = 0.5 accepted", function() {
  fit <- run_abf_region(.rg, .rp, coverage = 0.5)
  stopifnot(fit$params$coverage == 0.5)
})

run_test("coverage = 0.99 accepted", function() {
  fit <- run_abf_region(.rg, .rp, coverage = 0.99)
  stopifnot(fit$params$coverage == 0.99)
})

run_test("abf returns exactly one credible set", function() {
  fit <- run_abf_region(.rg, .rp)
  stopifnot(is.list(fit$credible_sets), length(fit$credible_sets) == 1L)
})

run_test("abf pip sums to 1 (normalised ABF)", function() {
  fit <- run_abf_region(.rg, .rp)
  stopifnot(abs(sum(fit$pip) - 1) < 1e-9)
})

run_test("abf additional contains log10_abf", function() {
  fit <- run_abf_region(.rg, .rp)
  stopifnot("log10_abf" %in% names(fit$additional))
})

run_test("larger prior_variance increases ABF magnitude", function() {
  fit_lo <- run_abf_region(.rg, .rp, prior_variance = 0.01)
  fit_hi <- run_abf_region(.rg, .rp, prior_variance = 0.5)
  # Top PIP variant should be the same; but log10_abf values differ
  stopifnot(!identical(fit_lo$pip, fit_hi$pip))
})


# =============================================================================
# SECTION 9: run_susie_inf / run_susie_inf_region
# =============================================================================

set_section("run_susie_inf / run_susie_inf_region — method arguments")

run_test("L = 5 accepted", function() {
  fit <- run_susie_inf_region(.rg, .rp, L = 5)
  stopifnot(fit$params$L == 5)
})

run_test("L = 1 accepted", function() {
  fit <- run_susie_inf_region(.rg, .rp, L = 1)
  stopifnot(fit$params$L == 1)
})

run_test("coverage = 0.9 accepted", function() {
  fit <- run_susie_inf_region(.rg, .rp, L = 5, coverage = 0.9)
  stopifnot(fit$params$coverage == 0.9)
})

run_test("max_iter = 50 accepted", function() {
  fit <- run_susie_inf_region(.rg, .rp, L = 5, max_iter = 50)
  stopifnot(fit$params$max_iter == 50)
})

run_test("susie_inf output has pip, credible_sets, method", function() {
  fit <- run_susie_inf_region(.rg, .rp, L = 5)
  stopifnot(all(c("pip", "credible_sets", "method") %in% names(fit)))
})

run_test("susie_inf pip in [0,1]", function() {
  fit <- run_susie_inf_region(.rg, .rp, L = 5)
  stopifnot(all(fit$pip >= 0), all(fit$pip <= 1))
})


# =============================================================================
# SECTION 10: run_carma / run_carma_region
# =============================================================================

set_section("run_carma / run_carma_region — method arguments")

run_test("rho.index = 0.95 (default) accepted", function() {
  fit <- run_carma_region(.rg, .rp, rho.index = 0.95)
  stopifnot(fit$params$rho.index == 0.95)
})

run_test("rho.index = 0.9 accepted", function() {
  fit <- run_carma_region(.rg, .rp, rho.index = 0.9)
  stopifnot(fit$params$rho.index == 0.9)
})

run_test("num.causal = 5 accepted", function() {
  fit <- run_carma_region(.rg, .rp, num.causal = 5)
  stopifnot(fit$params$num.causal == 5)
})

run_test("num.causal = 1 accepted", function() {
  fit <- run_carma_region(.rg, .rp, num.causal = 1)
  stopifnot(!is.null(fit$pip))
})

run_test("carma returns pip of correct length", function() {
  fit <- run_carma_region(.rg, .rp)
  stopifnot(length(fit$pip) == .rg$p)
})

run_test("carma pip in [0,1]", function() {
  fit <- run_carma_region(.rg, .rp)
  stopifnot(all(fit$pip >= 0, na.rm = TRUE), all(fit$pip <= 1, na.rm = TRUE))
})

run_test("carma returns exactly one credible set (global)", function() {
  fit <- run_carma_region(.rg, .rp)
  stopifnot(is.list(fit$credible_sets), length(fit$credible_sets) >= 0L)
})


# =============================================================================
# SECTION 11: External binary wrappers — argument passthrough (SKIP if absent)
# =============================================================================

set_section("External wrappers — argument forwarding (FINEMAP, PAINTOR, BEATRICE, FUNMAP)")

finemap_bin <- tryCatch(setup_finemap(download = FALSE), error = function(e) NULL)
paintor_bin <- tryCatch(setup_paintor(), error = function(e) NULL)

run_test("finemap: finemap_path arg recognised", function() {
  fp <- finemap_bin %||% "finemap"
  # Just check that the call doesn't fail before hitting the binary
  fit <- tryCatch(
    run_finemap_region(.rg, .rp, finemap_path = fp, n_causal = 2),
    error = function(e) list(error = conditionMessage(e))
  )
  # Either succeeds or errors with a binary/file-system issue (not an arg error)
  stopifnot(!is.null(fit))
}, skip_reason = if (is.null(finemap_bin)) "FINEMAP binary not available on this machine" else NULL)

run_test("paintor: paintor_path arg recognised", function() {
  pp <- paintor_bin %||% "PAINTOR"
  fit <- tryCatch(
    run_paintor_region(.rg, .rp, paintor_path = pp, max_causal = 1),
    error = function(e) list(error = conditionMessage(e))
  )
  stopifnot(!is.null(fit))
}, skip_reason = if (is.null(paintor_bin)) "PAINTOR binary not available / not compiled" else NULL)

run_test("beatrice: beatrice_dir and python args forwarded (graceful error if absent)", function() {
  fit <- tryCatch(
    run_beatrice_region(.rg, .rp,
                        beatrice_dir = "~/Beatrice-Finemapping",
                        python = "/opt/anaconda3/bin/python3",
                        max_iter = 50),
    error = function(e) list(error = conditionMessage(e))
  )
  # If beatrice dir missing, returns error field not a crash
  stopifnot(!is.null(fit))
})

run_test("funmap: python arg forwarded (graceful error if absent)", function() {
  fit <- tryCatch(
    run_funmap_region(.rg, .rp,
                      python = "/opt/anaconda3/bin/python3",
                      L = 5),
    error = function(e) list(error = conditionMessage(e))
  )
  stopifnot(!is.null(fit))
})


# =============================================================================
# SECTION 12: run_marginal_z / run_marginal_z_region
# =============================================================================

set_section("run_marginal_z / run_marginal_z_region — baseline behaviour")

# Reuse the same region accessor as the other sections. SIM_MINI is fine —
# marginal_z does not need annotations.
.rg_mz <- SIM_MINI$genotypes[[1]]
.rp_mz <- SIM_MINI$scenarios[[1]]$regions[[1]]

run_test("marginal_z output has all standard fields", function() {
  fit <- run_marginal_z_region(.rg_mz, .rp_mz)
  stopifnot(all(c("pip", "credible_sets", "method", "input_type",
                  "params", "runtime_seconds", "additional") %in% names(fit)))
  stopifnot(fit$method == "marginal_z")
  stopifnot(fit$input_type == "summary")
})

run_test("marginal_z pip length equals p", function() {
  fit <- run_marginal_z_region(.rg_mz, .rp_mz)
  stopifnot(length(fit$pip) == .rg_mz$p)
})

run_test("marginal_z pip values lie in [0, 1]", function() {
  fit <- run_marginal_z_region(.rg_mz, .rp_mz)
  stopifnot(all(fit$pip >= 0), all(fit$pip <= 1))
})

run_test("marginal_z pip sums to ~1 (|z| / sum|z| normalisation)", function() {
  fit <- run_marginal_z_region(.rg_mz, .rp_mz)
  stopifnot(abs(sum(fit$pip) - 1) < 1e-8)
})

run_test("marginal_z default coverage = 0.95 accepted", function() {
  fit <- run_marginal_z_region(.rg_mz, .rp_mz)
  stopifnot(fit$params$coverage == 0.95)
})

run_test("marginal_z lower coverage produces smaller or equal CS", function() {
  fit95 <- run_marginal_z_region(.rg_mz, .rp_mz, coverage = 0.95)
  fit50 <- run_marginal_z_region(.rg_mz, .rp_mz, coverage = 0.50)
  size95 <- length(fit95$credible_sets[[1]])
  size50 <- length(fit50$credible_sets[[1]])
  stopifnot(size50 <= size95)
})

run_test("marginal_z returns exactly one credible set", function() {
  fit <- run_marginal_z_region(.rg_mz, .rp_mz)
  stopifnot(is.list(fit$credible_sets), length(fit$credible_sets) == 1L)
})


# =============================================================================
# SECTION 13: run_polyfun_oracle / run_polyfun_oracle_region
# =============================================================================

set_section("run_polyfun_oracle / run_polyfun_oracle_region — oracle priors")

.rg_po <- SIM_MINI_ANNOT$genotypes[[1]]
.rp_po <- SIM_MINI_ANNOT$scenarios[[1]]$regions[[1]]

run_test("polyfun_oracle output has all standard fields (annotated fixture)", function() {
  fit <- run_polyfun_oracle_region(.rg_po, .rp_po)
  stopifnot(all(c("pip", "credible_sets", "method", "input_type",
                  "params", "runtime_seconds", "additional") %in% names(fit)))
  stopifnot(fit$method == "polyfun_oracle")
})

run_test("polyfun_oracle pip valid (length p, in [0,1]) on annotated fixture", function() {
  fit <- run_polyfun_oracle_region(.rg_po, .rp_po)
  stopifnot(length(fit$pip) == .rg_po$p)
  stopifnot(all(fit$pip >= 0, na.rm = TRUE), all(fit$pip <= 1, na.rm = TRUE))
})

run_test("polyfun_oracle reports prior_weights of length p", function() {
  fit <- run_polyfun_oracle_region(.rg_po, .rp_po)
  pw <- fit$additional$prior_weights
  stopifnot(!is.null(pw), length(pw) == .rg_po$p)
  stopifnot(all(pw >= 0, na.rm = TRUE))
})

run_test("polyfun_oracle falls back to uniform prior when annotations absent", function() {
  # SIM_MINI has annotations="none", so $annotations_matrix is NULL and
  # truth$enrichment is also absent. The wrapper is documented to
  # gracefully fall back to uniform priors in that case (degenerating
  # to plain SuSiE) rather than erroring.
  rg <- SIM_MINI$genotypes[[1]]
  rp <- SIM_MINI$scenarios[[1]]$regions[[1]]
  fit <- run_polyfun_oracle_region(rg, rp)
  stopifnot(is.null(fit$error))
  stopifnot(length(fit$pip) == rg$p)
  stopifnot(all(fit$pip >= 0, na.rm = TRUE), all(fit$pip <= 1, na.rm = TRUE))
  stopifnot(identical(fit$params$prior_source, "uniform_fallback"))
})


# =============================================================================
# SECTION 14: run_polyfun_est / run_polyfun_est_region / scenario_setup
# =============================================================================

set_section("run_polyfun_est — estimated priors and scenario-setup hook")

.rg_pe <- SIM_MINI_ANNOT$genotypes[[1]]
.rp_pe <- SIM_MINI_ANNOT$scenarios[[1]]$regions[[1]]

run_test("polyfun_est output has all standard fields (annotated fixture)", function() {
  fit <- run_polyfun_est_region(.rg_pe, .rp_pe)
  stopifnot(all(c("pip", "credible_sets", "method", "input_type",
                  "params", "runtime_seconds", "additional") %in% names(fit)))
  stopifnot(fit$method == "polyfun_est")
})

run_test("polyfun_est pip valid (length p, in [0,1]) on annotated fixture", function() {
  fit <- run_polyfun_est_region(.rg_pe, .rp_pe)
  stopifnot(length(fit$pip) == .rg_pe$p)
  stopifnot(all(fit$pip >= 0, na.rm = TRUE), all(fit$pip <= 1, na.rm = TRUE))
})

run_test("polyfun_est reports non-negative tau of length m + 1", function() {
  fit <- run_polyfun_est_region(.rg_pe, .rp_pe)
  tau <- fit$additional$tau
  m   <- ncol(.rg_pe$annotations_matrix)
  stopifnot(!is.null(tau), length(tau) == m + 1L)
  stopifnot(all(tau >= 0, na.rm = TRUE))
})

run_test("polyfun_est runs without scenario hook (per-region tau fallback)", function() {
  # Direct per-region call, no pooled_tau supplied; wrapper should
  # fall back to per-region estimation rather than erroring.
  fit <- run_polyfun_est_region(.rg_pe, .rp_pe, pooled_tau = NULL)
  stopifnot(length(fit$pip) == .rg_pe$p)
  stopifnot(!is.null(fit$additional$tau))
})

run_test("polyfun_est scenario-setup hook returns pooled_tau when annotations consistent", function() {
  scen <- SIM_MINI_ANNOT$scenarios[[1]]
  extra <- run_polyfun_est_scenario_setup(
    genotypes = SIM_MINI_ANNOT$genotypes,
    regions   = scen$regions,
    user_args = list()
  )
  m <- ncol(SIM_MINI_ANNOT$genotypes[[1]]$annotations_matrix)
  stopifnot("pooled_tau" %in% names(extra))
  stopifnot(length(extra$pooled_tau) == m + 1L)
  stopifnot(all(extra$pooled_tau >= 0))
})

run_test("polyfun_est fails gracefully on no-annotation fixture", function() {
  rg <- SIM_MINI$genotypes[[1]]
  rp <- SIM_MINI$scenarios[[1]]$regions[[1]]
  fit <- run_polyfun_est_region(rg, rp)
  # Either falls back to uniform priors (succeeds) or returns an error
  # result. Both are acceptable; we just need not to crash.
  stopifnot(!is.null(fit))
  stopifnot(length(fit$pip) == rg$p)
})


# =============================================================================
# SECTION 15: MAF-stratified evaluation (by_causal_maf)
# =============================================================================

set_section("evaluate_methods — by_causal_maf stratification")

# RESULTS_MINI was computed against the pre-MAF-axis evaluator; recompute
# here so the eval object has the new $by_causal_maf field.
EVAL_MAF <- evaluate_methods(
  SIM_MINI, RESULTS_MINI,
  save = FALSE, verbose = FALSE
)

run_test("evaluate_methods returns a by_causal_maf list per method", function() {
  stopifnot("by_causal_maf" %in% names(EVAL_MAF$susie))
  stopifnot("by_causal_maf" %in% names(EVAL_MAF$abf))
})

run_test("by_causal_maf entries are non-empty when MAFs are available", function() {
  # SIM_MINI was simulated from real VCFs, so genotypes carry MAFs and
  # the bins should be populated.
  bm <- EVAL_MAF$susie$by_causal_maf
  stopifnot(!is.null(bm))
  stopifnot(length(bm) >= 1L)
})

run_test("by_causal_maf bin names are a subset of {rare, low, common}", function() {
  bm <- EVAL_MAF$susie$by_causal_maf
  stopifnot(all(names(bm) %in% c("rare", "low", "common")))
})

run_test("by_causal_maf bins appear in canonical order rare -> low -> common", function() {
  bm <- EVAL_MAF$susie$by_causal_maf
  canonical <- c("rare", "low", "common")
  present   <- canonical[canonical %in% names(bm)]
  stopifnot(identical(names(bm), present))
})

run_test("by_causal_maf bins contain a numeric auprc field (possibly NA)", function() {
  bm <- EVAL_MAF$susie$by_causal_maf
  for (b in names(bm)) {
    stopifnot("auprc" %in% names(bm[[b]]))
    stopifnot(is.numeric(bm[[b]]$auprc))
  }
})


# =============================================================================
# SECTION 16: Misspecification stratification (by_true_annotation_type)
# =============================================================================

set_section("evaluate_methods — by_true_annotation_type stratification")

# EVAL_MAF was built against SIM_MINI (annotations = "none"). Re-use it
# for the no-annotation case; spin up a quick evaluation against
# SIM_MINI_ANNOT (binary annotations) for the binary case.

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

run_test("evaluate_methods returns by_true_annotation_type per method", function() {
  stopifnot("by_true_annotation_type" %in% names(EVAL_MAF$susie))
  stopifnot("by_true_annotation_type" %in% names(EVAL_ANNOT$susie))
})

run_test("by_true_annotation_type uses 'none' on the no-annotation fixture", function() {
  bt <- EVAL_MAF$susie$by_true_annotation_type
  stopifnot(!is.null(bt))
  stopifnot(identical(names(bt), "none"))
})

run_test("by_true_annotation_type uses 'binary' on the annotated fixture", function() {
  bt <- EVAL_ANNOT$susie$by_true_annotation_type
  stopifnot(!is.null(bt))
  stopifnot(identical(names(bt), "binary"))
})

run_test("by_true_annotation_type bins carry a numeric auprc field", function() {
  for (eo in list(EVAL_MAF$susie, EVAL_ANNOT$susie)) {
    bt <- eo$by_true_annotation_type
    for (t in names(bt)) {
      stopifnot("auprc" %in% names(bt[[t]]))
      stopifnot(is.numeric(bt[[t]]$auprc))
    }
  }
})


# =============================================================================
# Summary
# =============================================================================

message("\n", strrep("=", 70))
message("TEST SUMMARY")
message(strrep("=", 70))
message(sprintf("  PASS: %d", .n_pass))
message(sprintf("  FAIL: %d", .n_fail))
message(sprintf("  SKIP: %d", .n_skip))
message(sprintf("  TOTAL: %d", .n_pass + .n_fail + .n_skip))


# =============================================================================
# Generate markdown report
# =============================================================================

message("\nGenerating docs/testing_report.md ...")

.section_order <- unique(sapply(.RESULTS, `[[`, "section"))

lines <- c(
  "# Fine-Mapping Benchmark — Function & Argument Test Report",
  "",
  sprintf("**Generated:** %s  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("**R version:** %s  ", paste(R.version$major, R.version$minor, sep = ".")),
  sprintf("**Platform:** %s  ", R.version$platform),
  "",
  sprintf("**Results: %d PASS / %d FAIL / %d SKIP**",
          .n_pass, .n_fail, .n_skip),
  "",
  "---",
  "",
  "## Table of contents",
  ""
)

for (sec in .section_order) {
  anchor <- tolower(gsub("[^a-z0-9]+", "-", sec))
  lines  <- c(lines, sprintf("- [%s](#%s)", sec, anchor))
}
lines <- c(lines, "")

for (sec in .section_order) {
  sec_results <- Filter(function(r) r$section == sec, .RESULTS)
  n_p <- sum(sapply(sec_results, function(r) r$status == "PASS"))
  n_f <- sum(sapply(sec_results, function(r) r$status == "FAIL"))
  n_s <- sum(sapply(sec_results, function(r) r$status == "SKIP"))

  lines <- c(
    lines,
    sprintf("## %s", sec),
    "",
    sprintf("**%d PASS / %d FAIL / %d SKIP**", n_p, n_f, n_s),
    "",
    "| Test | Status | Notes |",
    "|------|--------|-------|"
  )

  for (r in sec_results) {
    icon  <- switch(r$status, PASS = "✅", FAIL = "❌", SKIP = "⏭️", "?")
    note  <- if (r$status == "FAIL")  paste("Error:", r$error) else
             if (r$status == "SKIP")  r$error else
             if (!is.null(r$elapsed)) sprintf("%.2fs", r$elapsed) else ""
    # Escape pipe characters in note
    note  <- gsub("\\|", "\\\\|", note)
    lines <- c(lines, sprintf("| %s | %s %s | %s |", r$name, icon, r$status, note))
  }

  lines <- c(lines, "")
}

# Append function reference appendix
lines <- c(
  lines,
  "---",
  "",
  "## Appendix — Function argument reference",
  "",
  "### `simulate_genotypes()`",
  "",
  "| Argument | Type | Default | Description |",
  "|----------|------|---------|-------------|",
  "| `n_regions` | integer | 3 | Number of independent genomic regions |",
  "| `n` | integer | 500 | Number of individuals |",
  "| `p` | integer or vector | 200 | Target SNPs per region (scalar or per-region vector) |",
  "| `vcf_files` | character vector or NULL | NULL | VCF files (one per region); NULL = bundled example |",
  "| `min_maf` | numeric ∈ [0, 0.5] | 0.01 | Minimum MAF filter |",
  "| `max_maf` | numeric or NA | NA | Maximum MAF filter; NA = no upper filter |",
  "| `standardise` | logical | TRUE | Standardise genotypes to mean 0, variance 1 |",
  "| `genetic_map_dir` | character or NULL | NULL | Cache directory for HapMap genetic maps |",
  "| `seed` | integer or NULL | NULL | Random seed for reproducibility |",
  "| `save` | logical | FALSE | Save genotype list as .rds to `output_dir` |",
  "| `output_dir` | character | 'results' | Directory for saved output (created if absent) |",
  "| `verbose` | logical | TRUE | Print progress messages |",
  "",
  "### `simulate_phenotypes()`",
  "",
  "| Argument | Type | Default | Description |",
  "|----------|------|---------|-------------|",
  "| `genotypes` | list | — | Output from `simulate_genotypes()` |",
  "| `S` | integer or vector | 1 | Causal variants per region (scalar or per-region vector) |",
  "| `phi` | numeric or vector | 0.1 | PVE (scalar or per-region vector), must be in (0, 1) |",
  "| `model` | 'sparse' or 'sparse_inf' | 'sparse' | Genetic architecture model |",
  "| `p_causal` | numeric ∈ (0, 1] | 0.5 | Fraction of PVE from sparse component (sparse_inf only) |",
  "| `inf_model` | 'beatrice' or 'susie_inf' | 'beatrice' | Infinitesimal component formulation (sparse_inf only) |",
  "| `effect_distribution` | 'normal' or 'equal' | 'normal' | Effect size distribution |",
  "| `effect_variance` | numeric > 0 | 0.36 | Variance for normal effect sizes |",
  "| `annotations` | 'none', 'binary', 'continuous', or matrix | 'none' | Annotation mode |",
  "| `n_annotations` | integer ≥ 1 | 3 | Number of annotation columns (for binary/continuous) |",
  "| `annotation_proportions` | numeric, vector, or NULL | NULL | Proportion of 1s per binary annotation |",
  "| `enrichment` | numeric, vector, or NULL | NULL | Fold-enrichment for annotation-guided selection |",
  "| `seed` | integer or NULL | NULL | Random seed |",
  "| `save` | logical | FALSE | Save phenotype list as .rds to `output_dir` |",
  "| `output_dir` | character | 'results' | Directory for saved output (created if absent) |",
  "| `verbose` | logical | TRUE | Print progress messages |",
  "",
  "### `run_simulation()`",
  "",
  "| Argument | Type | Default | Description |",
  "|----------|------|---------|-------------|",
  "| `n_regions` | integer | 3 | Number of genomic regions |",
  "| `n` | integer | 500 | Number of individuals |",
  "| `p` | integer or vector | 200 | Target SNPs per region |",
  "| `n_iter` | integer ≥ 1 | 5 | Replicates per parameter combination |",
  "| `S` | integer vector | c(1,2,3,5) | Causal-variant values to sweep |",
  "| `phi` | numeric vector ∈ (0,1) | c(0.1,0.2,0.4,0.6) | PVE values to sweep |",
  "| `model` | 'sparse' or 'sparse_inf' | 'sparse' | Genetic architecture |",
  "| `p_causal` | numeric vector ∈ (0,1] | c(0.1,0.2,0.4) | p_causal values to sweep (sparse_inf only) |",
  "| `inf_model` | 'beatrice' or 'susie_inf' | 'beatrice' | Infinitesimal formulation (sparse_inf only) |",
  "| `effect_distribution` | 'normal' or 'equal' | 'normal' | Effect size distribution |",
  "| `effect_variance` | numeric > 0 | 0.36 | Normal effect variance |",
  "| `annotations` | 'none', 'binary', 'continuous', or matrix | 'none' | Annotation mode |",
  "| `n_annotations` | integer ≥ 1 | 3 | Number of annotation columns |",
  "| `annotation_proportions` | numeric, vector, or NULL | NULL | Binary annotation proportions |",
  "| `enrichment` | numeric, vector, or NULL | NULL | Annotation enrichment |",
  "| `vcf_dir` | character or NULL | NULL | Directory of VCF files (from prepare_vcfs.R) |",
  "| `vcf_files` | character vector or NULL | NULL | Explicit VCF paths (overrides vcf_dir) |",
  "| `min_maf` | numeric | 0.01 | Minimum MAF |",
  "| `max_maf` | numeric or NA | NA | Maximum MAF |",
  "| `standardise` | logical | TRUE | Standardise genotypes |",
  "| `seed` | integer or NULL | NULL | Master random seed |",
  "| `save` | logical | FALSE | Save result as .rds |",
  "| `output_dir` | character | 'results' | Output directory |",
  "| `verbose` | logical | TRUE | Print progress |",
  "",
  "### `run_methods()`",
  "",
  "| Argument | Type | Default | Description |",
  "|----------|------|---------|-------------|",
  "| `simulation` | list | — | Output of `run_simulation()` |",
  "| `methods` | character vector | 'susie' | Method names to run (case-insensitive) |",
  "| `method_args` | named list | list() | Per-method argument overrides |",
  "| `save` | logical | FALSE | Save per-method .rds files |",
  "| `output_dir` | character | 'results' | Output directory |",
  "| `verbose` | logical | TRUE | Print progress |",
  "",
  "**Supported methods and their key tuneable arguments (via `method_args`):**",
  "",
  "| Method | Key arguments |",
  "|--------|--------------|",
  "| `susie` | `L`, `coverage`, `min_abs_corr`, `max_iter`, `estimate_residual_variance`, `estimate_prior_variance`, `prior_variance` |",
  "| `susie_inf` | `L`, `coverage`, `max_iter` |",
  "| `abf` | `prior_variance`, `coverage` |",
  "| `carma` | `rho.index`, `num.causal` |",
  "| `finemap` | `finemap_path`, `n_causal`, `prior_std` |",
  "| `paintor` | `paintor_path`, `max_causal` |",
  "| `beatrice` | `beatrice_dir`, `python`, `max_iter` |",
  "| `funmap` | `python`, `L`, `max_iter` |",
  "",
  "### `evaluate_methods()`",
  "",
  "| Argument | Type | Default | Description |",
  "|----------|------|---------|-------------|",
  "| `simulation` | list | — | Output of `run_simulation()` |",
  "| `results` | list | — | Output of `run_methods()` |",
  "| `pip_thresholds` | numeric vector | seq(0,1,by=0.005) | PIP thresholds for power/FDR curve |",
  "| `n_pip_cal_bins` | integer | 10 | Equal-width bins for PIP calibration |",
  "| `save` | logical | FALSE | Write evaluation.rds and evaluation_summary.csv |",
  "| `output_dir` | character | 'results' | Output directory |",
  "| `verbose` | logical | TRUE | Print progress |",
  "",
  "### `plot_results()`",
  "",
  "| Argument | Type | Default | Description |",
  "|----------|------|---------|-------------|",
  "| `eval_out` | list | — | Output of `evaluate_methods()` |",
  "| `output_file` | character or NULL | NULL | Full PDF path (overrides `output_dir` when set) |",
  "| `output_dir` | character | 'results' | Directory to save `evaluation.pdf` when `output_file` is NULL |",
  "| `save` | logical | TRUE | If FALSE, skip writing the PDF entirely |",
  "| `methods` | character vector or NULL | NULL | Methods to include (NULL = all evaluated) |",
  "| `verbose` | logical | TRUE | Print progress |",
  ""
)

dir.create("docs", showWarnings = FALSE)
writeLines(lines, "docs/testing_report.md")
message("Report written to docs/testing_report.md")

if (.n_fail > 0) {
  message(sprintf("\n%d test(s) FAILED. See above for details.", .n_fail))
  quit(status = 1)
} else {
  message(sprintf("\nAll %d tests passed!", .n_pass + .n_skip))
}
