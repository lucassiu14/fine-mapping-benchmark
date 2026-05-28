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
#' @param vcf_dir Character or NULL. Path to a directory of VCF files prepared
#'   by \code{inst/scripts/prepare_vcfs.R} (e.g. \code{"data/vcf"}). If provided,
#'   \code{n_regions} files are sampled at random from this directory for each
#'   run (reproducibly if \code{seed} is set). Each file provides a distinct
#'   genomic region with real LD structure from 1000 Genomes Phase 3. Ignored
#'   when \code{vcf_files} is also supplied. Default: NULL.
#' @param vcf_files Character vector or NULL. VCF files for genotype
#'   simulation. Overrides \code{vcf_dir}. Default: NULL (use bundled example).
#' @param min_maf Numeric. Minimum MAF filter. Default: 0.01.
#' @param max_maf Numeric or NA. Maximum MAF filter. Default: NA.
#' @param standardise Logical. Standardise genotypes. Default: TRUE.
#' @param seed Integer or NULL. Master random seed. Default: NULL.
#' @param save Logical. If TRUE, save the result to \code{output_dir} as an
#'   \code{.rds} file. The filename encodes the key simulation parameters
#'   (model, n, p, S values, n_iter) and the seed (or "noseed" if NULL).
#'   Default: FALSE.
#' @param output_dir Character. Directory in which to save the result when
#'   \code{save = TRUE}. Created automatically if it does not exist.
#'   Default: \code{"results"}.
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
#'
#' # Save result to disk
#' result <- run_simulation(
#'   n_regions = 2, n = 200, p = 100, seed = 42,
#'   save = TRUE, output_dir = "results"
#' )
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
                           vcf_dir = NULL,
                           vcf_files = NULL,
                           genetic_map_dir = "data/genetic_maps",
                           min_maf = 0.01,
                           max_maf = NA,
                           standardise = TRUE,
                           seed = NULL,
                           save = FALSE,
                           output_dir = "results",
                           verbose = TRUE,
                           n_ref = NULL) {

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

  # --- Resolve vcf_dir → vcf_files --------------------------------------------

  if (!is.null(vcf_dir) && is.null(vcf_files)) {
    if (!dir.exists(vcf_dir)) {
      stop(
        "vcf_dir does not exist: ", vcf_dir, "\n",
        "Run inst/scripts/prepare_vcfs.R first to download the reference VCF files.",
        call. = FALSE
      )
    }
    available <- list.files(vcf_dir, pattern = "\\.vcf\\.gz$", full.names = TRUE)
    # Exclude index files (*.vcf.gz.tbi)
    available <- available[!grepl("\\.tbi$", available)]
    if (length(available) == 0) {
      stop(
        "No .vcf.gz files found in vcf_dir: ", vcf_dir, "\n",
        "Run inst/scripts/prepare_vcfs.R (development checkout) or ",
        "system.file(\"scripts/prepare_vcfs.R\", package = \"fmbenchmark\") first.",
        call. = FALSE
      )
    }
    if (length(available) < n_regions) {
      stop(
        sprintf(
          "vcf_dir contains %d VCF file(s) but n_regions = %d. ",
          length(available), n_regions
        ),
        "Add more regions to inst/extdata/regions.csv and re-run ",
        "inst/scripts/prepare_vcfs.R, or reduce n_regions.",
        call. = FALSE
      )
    }
    vcf_files <- sort(sample(available, n_regions))
    if (verbose) {
      message(sprintf(
        "Sampled %d region(s) from %d available in %s",
        n_regions, length(available), vcf_dir
      ))
    }
  }

  # --- Memory guard for LD matrix storage -------------------------------------
  # Each LD matrix is p_i x p_i doubles (8 bytes/element). When n_ref is set
  # we keep BOTH the in-sample LD (LD_true) and the ref-panel-derived LD per
  # region, doubling the LD-side memory footprint.
  ld_factor <- if (is.null(n_ref)) 1 else 2
  ld_bytes  <- 8 * ld_factor * sum(as.numeric(p)^2)
  ld_gb     <- ld_bytes / 1e9
  if (ld_gb > 4) {
    warning(sprintf(
      "Combined LD matrices will use ~%.1f GB%s. Reduce p, reduce n_regions, or run on a high-memory node.",
      ld_gb,
      if (!is.null(n_ref)) " (n_ref doubles the LD footprint)" else ""
    ), call. = FALSE)
  }

  # --- Step 1: Simulate genotypes (once) --------------------------------------

  if (verbose) message("\n=== Simulating genotypes ===")

  if (!is.null(genetic_map_dir) && !dir.exists(genetic_map_dir)) {
    dir.create(genetic_map_dir, recursive = TRUE)
  }

  genotypes <- simulate_genotypes(
    n_regions = n_regions,
    n = n,
    p = p,
    vcf_files = vcf_files,
    min_maf = min_maf,
    max_maf = max_maf,
    standardise = standardise,
    genetic_map_dir = genetic_map_dir,
    seed = NULL,  # seed already set above
    verbose = verbose,
    n_ref = n_ref
  )

  # Pre-compute LD matrices (shared across all scenarios).
  #
  # - LD_true is always the in-sample correlation matrix: cor(X). Used for
  #   diagnostics (e.g. mean((LD - LD_true)^2)) and as the "perfect" LD a
  #   method would receive if it had the GWAS genotypes themselves.
  # - LD is what methods actually receive. When n_ref is set we use the
  #   reference panel: cor(X_ref). Otherwise LD == LD_true exactly,
  #   preserving the pre-n_ref behaviour.
  if (verbose) message("\nPre-computing LD matrices...")
  for (i in seq_len(n_regions)) {
    genotypes[[i]]$LD_true <- cor(genotypes[[i]]$X)
    if (!is.null(genotypes[[i]]$X_ref)) {
      genotypes[[i]]$LD <- cor(genotypes[[i]]$X_ref)
    } else {
      genotypes[[i]]$LD <- genotypes[[i]]$LD_true
    }
  }

  # Pre-generate annotation matrices (once per region, shared across all scenarios)
  # Only applies to synthetic annotations; user-supplied matrices are already fixed.
  annotation_type_internal <- if (is.matrix(annotations)) "user_supplied" else annotations
  if (annotation_type_internal %in% c("binary", "continuous")) {
    if (verbose) message("\nPre-generating annotation matrices...")
    # Broadcast scalar annotation_proportions to the required length so that
    # simulate_annotations_for_region receives a proper vector (not a scalar
    # that would produce NA when indexed beyond position 1).
    ap_internal <- if (!is.null(annotation_proportions) && length(annotation_proportions) == 1L)
      rep(annotation_proportions, n_annotations) else annotation_proportions
    for (i in seq_len(n_regions)) {
      annot <- simulate_annotations_for_region(
        p = genotypes[[i]]$p,
        annotation_type = annotation_type_internal,
        n_annotations = n_annotations,
        annotation_proportions = ap_internal,
        user_annotation_matrix = NULL
      )
      genotypes[[i]]$annotations_matrix    <- annot$matrix
      genotypes[[i]]$annotation_proportions <- annot$proportions
    }
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

  # Realised per-region p may differ from the requested p when the underlying
  # VCF contains fewer variants than requested (after MAF filter / subsetting).
  # Record both so downstream analysis is unambiguous.
  p_actual <- vapply(genotypes, function(g) g$p, integer(1))

  params <- list(
    n_regions = n_regions,
    n = n,
    p = p,
    p_actual = p_actual,
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
    vcf_files_used = vcf_files,
    seed = seed,
    n_scenarios = n_scenarios,
    n_ref = n_ref
  )

  # --- Assemble result --------------------------------------------------------

  result <- list(
    genotypes = genotypes,
    scenarios = scenarios,
    params = params
  )

  # --- Save to disk (optional) ------------------------------------------------

  if (save) {
    stopifnot(
      "output_dir must be a single character string" =
        is.character(output_dir) && length(output_dir) == 1
    )
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    seed_tag  <- if (!is.null(seed)) paste0("seed", seed) else "noseed"
    S_tag     <- paste(S, collapse = "-")
    p_tag     <- if (length(unique(p)) == 1) as.character(p[1]) else paste(p, collapse = "-")
    fname     <- sprintf(
      "simulation_%s_n%d_p%s_S%s_iter%d_%s.rds",
      model, n, p_tag, S_tag, n_iter, seed_tag
    )
    fpath <- file.path(output_dir, fname)
    saveRDS(result, file = fpath)

    if (verbose) {
      message(sprintf("Result saved to: %s", fpath))
    }
  }

  # --- Return -----------------------------------------------------------------

  if (verbose) {
    message(sprintf(
      "\nDone. %d scenarios generated across %d region(s).",
      n_scenarios, n_regions
    ))
  }

  result
}
