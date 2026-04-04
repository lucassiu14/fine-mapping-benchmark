# =============================================================================
# run_simulation.R
#
# Orchestration wrapper for benchmarking simulations.
#
# Simulates genotypes once, then iterates over a grid of parameter
# combinations (S, phi, and optionally p_causal) to produce multiple
# phenotype replicates for each setting.
# =============================================================================


#' Run a full benchmarking simulation
#'
#' Simulates genotypes once, then generates phenotypes across a grid of
#' parameter combinations (number of causal variants S, proportion of
#' variance explained phi, and optionally proportion of genetic variance
#' from sparse component p_causal). For each parameter combination,
#' \code{n_iter} independent replicates are generated.
#'
#' @param n_regions Integer. Number of independent genomic regions.
#'   Default: 3.
#' @param n Integer. Number of individuals. Default: 500.
#' @param p Integer or integer vector. Number of SNPs per region.
#'   Default: 200.
#' @param n_iter Integer. Number of independent replicates per parameter
#'   combination. Default: 5.
#' @param S Integer or integer vector. Number of causal variants to sweep
#'   over. Default: c(1, 2, 3, 5).
#' @param phi Numeric or numeric vector. PVE values to sweep over.
#'   Default: c(0.1, 0.2, 0.4, 0.6).
#' @param model Character. \code{"sparse"} or \code{"sparse_inf"}.
#'   Default: "sparse".
#' @param p_causal Numeric or numeric vector. Proportion of genetic variance
#'   from the sparse component. Only used when \code{model = "sparse_inf"}.
#'   If a vector, all values are swept. Default: c(0.1, 0.2, 0.4).
#' @param inf_model Character. Infinitesimal formulation: \code{"beatrice"}
#'   or \code{"susie_inf"}. Only used when \code{model = "sparse_inf"}.
#'   Default: "beatrice".
#' @param effect_distribution Character. \code{"normal"} or \code{"equal"}.
#'   Default: "normal".
#' @param effect_variance Numeric. Variance for normal effect sizes.
#'   Default: 0.36.
#' @param annotations Character or matrix. \code{"none"}, \code{"binary"},
#'   \code{"continuous"}, or a user-supplied matrix. Default: "none".
#' @param n_annotations Integer. Number of annotations (if applicable).
#'   Default: 3.
#' @param annotation_proportions Numeric or NULL. Proportions for binary
#'   annotations. Default: NULL (random).
#' @param enrichment Numeric or NULL. Fold-enrichment for annotations.
#'   Default: NULL (random).
#' @param vcf_files Character vector or NULL. VCF files for genotype
#'   simulation. Default: NULL (use bundled example).
#' @param min_maf Numeric. Minimum MAF filter. Default: 0.01.
#' @param max_maf Numeric or NA. Maximum MAF filter. Default: NA.
#' @param standardise Logical. Standardise genotypes. Default: TRUE.
#' @param seed Integer or NULL. Master random seed. Default: NULL.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{genotypes}{Output from \code{simulate_genotypes()}: a list of
#'       regions, each containing X, X_raw, maf, variant_ids, etc. The LD
#'       matrix is also pre-computed and stored here (shared across all
#'       scenarios).}
#'     \item{scenarios}{A list of simulation scenarios. Each scenario is a
#'       list containing:
#'       \describe{
#'         \item{scenario_id}{Unique integer ID.}
#'         \item{S}{Number of causal variants.}
#'         \item{phi}{PVE.}
#'         \item{p_causal}{Proportion of genetic variance from sparse
#'           component (NULL for sparse model).}
#'         \item{iter}{Replicate number (1 to n_iter).}
#'         \item{model}{Genetic architecture model.}
#'         \item{regions}{A list (one per region) each containing:
#'           z, beta_hat, se, y, annotations_matrix, and truth.}
#'       }
#'     }
#'     \item{params}{A list recording all simulation parameters for
#'       reproducibility.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Sparse model with defaults
#' result <- run_simulation(n_regions = 2, n = 200, p = 100, seed = 42)
#'
#' # How many scenarios?
#' length(result$scenarios)  # 4 * 4 * 5 = 80
#'
#' # Access first scenario
#' result$scenarios[[1]]$S
#' result$scenarios[[1]]$phi
#' result$scenarios[[1]]$regions[[1]]$z
#'
#' # Sparse + infinitesimal
#' result_inf <- run_simulation(
#'   n_regions = 2, n = 200, p = 100,
#'   model = "sparse_inf", seed = 42
#' )
#' length(result_inf$scenarios)  # 4 * 4 * 3 * 5 = 240
#' }
#'
#' @export
run_simulation <- function(n_regions = 3,
                           n = 500,
                           p = 200,
                           n_iter = 5,
                           S = c(1, 2, 3, 5),
                           phi = c(0.1, 0.2, 0.4, 0.6),
                           model = "sparse",
                           p_causal = c(0.1, 0.2, 0.4),
                           inf_model = "beatrice",
                           effect_distribution = "normal",
                           effect_variance = 0.36,
                           annotations = "none",
                           n_annotations = 3,
                           annotation_proportions = NULL,
                           enrichment = NULL,
                           vcf_files = NULL,
                           min_maf = 0.01,
                           max_maf = NA,
                           standardise = TRUE,
                           seed = NULL,
                           verbose = TRUE) {

  # --- Input validation -------------------------------------------------------

  model <- match.arg(model, choices = c("sparse", "sparse_inf"))
  inf_model <- match.arg(inf_model, choices = c("beatrice", "susie_inf"))
  effect_distribution <- match.arg(effect_distribution, choices = c("normal", "equal"))

  stopifnot(
    "n_iter must be a positive integer" =
      is.numeric(n_iter) && length(n_iter) == 1 &&
      n_iter == floor(n_iter) && n_iter >= 1
  )
  n_iter <- as.integer(n_iter)

  stopifnot(
    "S must be a vector of positive integers" =
      is.numeric(S) && all(S == floor(S)) && all(S >= 1)
  )
  S <- as.integer(S)

  stopifnot(
    "phi must be a numeric vector with values in (0, 1)" =
      is.numeric(phi) && all(phi > 0) && all(phi < 1)
  )

  if (model == "sparse_inf") {
    stopifnot(
      "p_causal must be a numeric vector with values in (0, 1]" =
        is.numeric(p_causal) && all(p_causal > 0) && all(p_causal <= 1)
    )
  }

  # --- Set master seed --------------------------------------------------------

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Build parameter grid ---------------------------------------------------

  if (model == "sparse") {
    param_grid <- expand.grid(
      S = S,
      phi = phi,
      iter = seq_len(n_iter),
      stringsAsFactors = FALSE
    )
    param_grid$p_causal <- NA_real_
  } else {
    param_grid <- expand.grid(
      S = S,
      phi = phi,
      p_causal = p_causal,
      iter = seq_len(n_iter),
      stringsAsFactors = FALSE
    )
  }

  n_scenarios <- nrow(param_grid)

  if (verbose) {
    if (model == "sparse") {
      message(sprintf(
        "Simulation grid: %d S values x %d phi values x %d iterations = %d scenarios",
        length(S), length(phi), n_iter, n_scenarios
      ))
    } else {
      message(sprintf(
        "Simulation grid: %d S x %d phi x %d p_causal x %d iter = %d scenarios",
        length(S), length(phi), length(p_causal), n_iter, n_scenarios
      ))
    }
    message(sprintf("Each scenario has %d region(s).", n_regions))
    message(sprintf("Total phenotype simulations: %d", n_scenarios * n_regions))
  }

  # --- Step 1: Simulate genotypes (once) --------------------------------------

  if (verbose) message("\n=== Simulating genotypes ===")

  genotypes <- simulate_genotypes(
    n_regions = n_regions,
    n = n,
    p = p,
    vcf_files = vcf_files,
    min_maf = min_maf,
    max_maf = max_maf,
    standardise = standardise,
    seed = NULL,  # seed already set above
    verbose = verbose
  )

  # Pre-compute LD matrices (shared across all scenarios)
  if (verbose) message("\nPre-computing LD matrices...")
  for (i in seq_len(n_regions)) {
    genotypes[[i]]$LD <- cor(genotypes[[i]]$X)
  }

  # Check that S values are feasible for all regions
  min_p <- min(vapply(genotypes, function(r) r$p, integer(1)))
  if (max(S) > min_p) {
    stop(
      sprintf(
        "Max S = %d exceeds the number of SNPs in the smallest region (%d). ",
        max(S), min_p
      ),
      "Reduce S or increase p.",
      call. = FALSE
    )
  }

  # --- Step 2: Generate phenotypes for each scenario --------------------------

  if (verbose) message("\n=== Simulating phenotypes ===")

  scenarios <- vector("list", n_scenarios)

  for (sc in seq_len(n_scenarios)) {

    S_sc <- param_grid$S[sc]
    phi_sc <- param_grid$phi[sc]
    iter_sc <- param_grid$iter[sc]
    p_causal_sc <- if (model == "sparse_inf") param_grid$p_causal[sc] else NULL

    if (verbose) {
      if (model == "sparse") {
        message(sprintf(
          "Scenario %d/%d: S=%d, phi=%.2f, iter=%d",
          sc, n_scenarios, S_sc, phi_sc, iter_sc
        ))
      } else {
        message(sprintf(
          "Scenario %d/%d: S=%d, phi=%.2f, p_causal=%.2f, iter=%d",
          sc, n_scenarios, S_sc, phi_sc, p_causal_sc, iter_sc
        ))
      }
    }

    # Run simulate_phenotypes on a copy of genotypes
    sim_result <- simulate_phenotypes(
      genotypes = genotypes,
      S = S_sc,
      phi = phi_sc,
      model = model,
      p_causal = if (!is.null(p_causal_sc)) p_causal_sc else 0.5,
      inf_model = inf_model,
      effect_distribution = effect_distribution,
      effect_variance = effect_variance,
      annotations = annotations,
      n_annotations = n_annotations,
      annotation_proportions = annotation_proportions,
      enrichment = enrichment,
      seed = NULL,  # let the RNG continue sequentially
      verbose = FALSE
    )

    # Extract only the phenotype-related fields (not the genotypes again)
    region_results <- vector("list", n_regions)
    for (i in seq_len(n_regions)) {
      region_results[[i]] <- list(
        y = sim_result[[i]]$y,
        z = sim_result[[i]]$z,
        beta_hat = sim_result[[i]]$beta_hat,
        se = sim_result[[i]]$se,
        annotations_matrix = sim_result[[i]]$annotations_matrix,
        truth = sim_result[[i]]$truth
      )
    }

    scenarios[[sc]] <- list(
      scenario_id = sc,
      S = S_sc,
      phi = phi_sc,
      p_causal = p_causal_sc,
      iter = iter_sc,
      model = model,
      regions = region_results
    )
  }

  # --- Store parameters for reproducibility -----------------------------------

  params <- list(
    n_regions = n_regions,
    n = n,
    p = p,
    n_iter = n_iter,
    S_values = S,
    phi_values = phi,
    model = model,
    p_causal_values = if (model == "sparse_inf") p_causal else NULL,
    inf_model = if (model == "sparse_inf") inf_model else NULL,
    effect_distribution = effect_distribution,
    effect_variance = effect_variance,
    annotation_type = if (is.character(annotations)) annotations else "user_supplied",
    n_annotations = n_annotations,
    annotation_proportions = annotation_proportions,
    enrichment = enrichment,
    min_maf = min_maf,
    seed = seed,
    n_scenarios = n_scenarios
  )

  # --- Return -----------------------------------------------------------------

  if (verbose) {
    message(sprintf(
      "\nDone. %d scenarios generated across %d region(s).",
      n_scenarios, n_regions
    ))
  }

  list(
    genotypes = genotypes,
    scenarios = scenarios,
    params = params
  )
}
