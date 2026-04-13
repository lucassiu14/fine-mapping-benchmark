# =============================================================================
# simulate_gwfm_data.R
#
# Genome-wide fine-mapping simulation framework.
#
# Unlike the locus-based simulate_genotypes() / simulate_phenotypes() pipeline,
# this function simulates a SINGLE shared phenotype across all genomic regions,
# consistent with how a real GWAS works:
#
#   y = X_1 β_1 + X_2 β_2 + ... + X_K β_K + ε
#
# Causal variants are assigned genome-wide via Bernoulli(π) rather than
# fixing exactly S causal variants per region. This means some regions will
# have zero causal variants (null regions), which is realistic and important
# for evaluating false-positive control.
#
# Output format mirrors run_simulation() so that run_gwfm_methods() and
# evaluate_methods() can consume it without modification, with the exception
# that scenario metadata uses `pi` and `h2` instead of `S` and `phi`, and
# each region's truth includes `S_realized` (which may be 0).
# =============================================================================


#' Simulate genome-wide fine-mapping data
#'
#' Simulates genotypes across a representative set of genomic regions, then
#' generates a single shared phenotype (consistent with a real GWAS) for each
#' combination of genome-wide polygenicity (\code{pi}) and heritability
#' (\code{h2}). Supports sparse and sparse-infinitesimal genetic architectures
#' and functional annotations.
#'
#' @param n Integer. Number of individuals to simulate.
#' @param n_iter Integer. Number of independent phenotype replicates per
#'   \code{(pi, h2)} combination. Genotypes are simulated once and reused.
#'   Default: 5.
#' @param pi Numeric vector. Genome-wide polygenicity values to sweep over
#'   (probability that any given variant is causal). Default: \code{c(1e-4, 1e-3)}.
#' @param h2 Numeric vector. Total SNP heritability values to sweep over.
#'   Must be in (0, 1). Default: \code{c(0.1, 0.3)}.
#' @param model Character. Genetic architecture model. \code{"sparse"} (point
#'   mass at zero + normal slab) or \code{"sparse_inf"} (sparse component plus
#'   genome-wide infinitesimal background). Default: \code{"sparse"}.
#' @param p_causal Numeric. Proportion of total heritability attributable to
#'   the sparse (causal) component. Only used when \code{model = "sparse_inf"}.
#'   Must be in (0, 1]. Default: 0.5.
#' @param inf_model Character. Infinitesimal formulation. \code{"beatrice"}
#'   places infinitesimal effects on non-causal variants only; \code{"susie_inf"}
#'   places them on all variants. Only used when \code{model = "sparse_inf"}.
#'   Default: \code{"beatrice"}.
#' @param effect_distribution Character. Distribution for causal effect sizes.
#'   \code{"normal"} draws from N(0, \code{effect_variance}).
#'   \code{"equal"} assigns equal magnitudes (signed randomly).
#'   Default: \code{"normal"}.
#' @param effect_variance Numeric. Variance of the normal effect size
#'   distribution (before genome-wide scaling to \code{h2}). Only used when
#'   \code{effect_distribution = "normal"}. Default: 0.36.
#' @param annotations Character or matrix. \code{"none"}, \code{"binary"},
#'   \code{"continuous"}, or a user-supplied p_total x m matrix. Annotations
#'   modify per-variant causal probability via enrichment weights.
#'   Default: \code{"none"}.
#' @param n_annotations Integer. Number of annotation columns. Only used when
#'   \code{annotations} is \code{"binary"} or \code{"continuous"}.
#'   Default: 3.
#' @param annotation_proportions Numeric, scalar, vector, or NULL. Proportion
#'   of variants with value 1 per binary annotation column. If NULL, drawn
#'   from Uniform(0.01, 0.30). Only used when \code{annotations = "binary"}.
#'   Default: NULL.
#' @param enrichment Numeric, scalar, vector, or NULL. Fold-enrichment of each
#'   annotation for causal variant selection. If NULL, drawn from
#'   Uniform(2, 10). Values must be > 0. Default: NULL.
#' @param regions Character or data frame. \code{"representative"} loads the
#'   128 bundled genome-wide LD blocks from \code{data/gwfm_regions.csv}.
#'   A chromosome name (e.g. \code{"chr1"} or \code{"1"}) loads only blocks
#'   from that chromosome (useful for small-scale testing). A data frame with
#'   columns \code{region_id}, \code{chrom}, \code{start}, \code{end} uses
#'   user-supplied coordinates (e.g. a full LDetect partition loaded from
#'   \code{data/gwfm_regions_ldetect_EUR.csv}). Default: \code{"representative"}.
#' @param coverage Numeric in \code{(0, 1]}. Fraction of the loaded region set
#'   to use. Regions are subsampled \emph{proportionally within each chromosome}
#'   so that genome-wide spread is preserved at all coverage levels. For each
#'   chromosome, \code{round(coverage * n_chr_regions)} regions are sampled
#'   without replacement; chromosomes where this rounds to zero are omitted.
#'   \code{coverage = 1} (default) uses all available regions. Lower values
#'   reduce compute and download requirements proportionally. When using the
#'   bundled 128-region set: \code{coverage = 0.25} gives ~32 regions;
#'   \code{coverage = 0.1} gives ~13 regions. When using a full LDetect
#'   partition (~1,703 blocks), \code{coverage = 1} approximates the entire
#'   autosome. Default: \code{1}.
#' @param p Integer or integer vector. Target number of SNPs per region.
#'   Default: 200.
#' @param min_maf Numeric. Minimum minor allele frequency filter. Default: 0.01.
#' @param vcf_dir Character. Directory containing per-region VCF files
#'   (produced by \code{scripts/prepare_gwfm_vcfs.R}). Each file should be
#'   named \code{<region_id>.vcf.gz}. Default: \code{"data/gwfm_vcf"}.
#' @param genetic_map_dir Character or NULL. Directory for caching HapMap
#'   GRCh37 genetic maps downloaded by sim1000G. Default:
#'   \code{"data/genetic_maps"}.
#' @param seed Integer or NULL. Master random seed. Default: NULL.
#' @param save Logical. If TRUE, save the returned list as an \code{.rds} file
#'   inside \code{output_dir}. Default: FALSE.
#' @param output_dir Character. Directory for saved output. Created
#'   automatically if absent. Default: \code{"results"}.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list with three elements, mirroring the output of
#'   \code{\link{run_simulation}}:
#'   \describe{
#'     \item{genotypes}{List of length \code{n_regions}. Each element contains
#'       \code{X}, \code{X_raw}, \code{n}, \code{p}, \code{maf},
#'       \code{variant_ids}, \code{region_id}, \code{vcf_source}, \code{LD},
#'       and (if applicable) \code{annotations_matrix}.}
#'     \item{scenarios}{List of length \code{n_pi x n_h2 x n_iter} (x
#'       \code{n_p_causal} for \code{sparse_inf}). Each scenario contains:
#'       \describe{
#'         \item{scenario_id}{Integer.}
#'         \item{pi}{Genome-wide polygenicity used.}
#'         \item{h2}{Total SNP heritability used.}
#'         \item{p_causal}{Sparse fraction (NULL for sparse model).}
#'         \item{iter}{Replicate index.}
#'         \item{model}{Architecture model used.}
#'         \item{S_total}{Total number of causal variants realised genome-wide.}
#'         \item{regions}{List of length \code{n_regions}, each containing:
#'           \code{z}, \code{beta_hat}, \code{se}, \code{y} (shared phenotype),
#'           and \code{truth} (list with \code{causal_indices},
#'           \code{causal_effects}, \code{beta_true}, \code{pve_region},
#'           \code{S_realized}, \code{pi}, \code{h2}).}
#'       }
#'     }
#'     \item{params}{List recording all simulation parameters.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Small-scale test: chromosome 1 only, 2 iterations
#' sim <- simulate_gwfm_data(
#'   n = 500, n_iter = 2,
#'   pi = 1e-3, h2 = 0.3,
#'   regions = "1",
#'   vcf_dir = "data/gwfm_vcf",
#'   seed = 42
#' )
#'
#' # Full representative set, sparse model
#' sim <- simulate_gwfm_data(
#'   n = 5000, n_iter = 5,
#'   pi = c(1e-4, 1e-3), h2 = c(0.1, 0.3, 0.5),
#'   vcf_dir = "data/gwfm_vcf",
#'   seed = 1
#' )
#'
#' # Sparse + infinitesimal with binary annotations
#' sim <- simulate_gwfm_data(
#'   n = 2000, n_iter = 3,
#'   pi = 1e-3, h2 = 0.3,
#'   model = "sparse_inf", p_causal = 0.7,
#'   annotations = "binary", n_annotations = 3,
#'   regions = "1",
#'   vcf_dir = "data/gwfm_vcf",
#'   seed = 7
#' )
#' }
#'
#' @export
simulate_gwfm_data <- function(n,
                                n_iter              = 5,
                                pi                  = c(1e-4, 1e-3),
                                h2                  = c(0.1, 0.3),
                                model               = "sparse",
                                p_causal            = 0.5,
                                inf_model           = "beatrice",
                                effect_distribution = "normal",
                                effect_variance     = 0.36,
                                annotations         = "none",
                                n_annotations       = 3,
                                annotation_proportions = NULL,
                                enrichment          = NULL,
                                regions             = "representative",
                                coverage            = 1,
                                p                   = 200,
                                min_maf             = 0.01,
                                vcf_dir             = "data/gwfm_vcf",
                                genetic_map_dir     = "data/genetic_maps",
                                seed                = NULL,
                                save                = FALSE,
                                output_dir          = "results",
                                verbose             = TRUE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  if (!requireNamespace("sim1000G", quietly = TRUE)) {
    stop("Package 'sim1000G' is required. Install with: install.packages('sim1000G')",
         call. = FALSE)
  }

  stopifnot(
    "n must be a positive integer"     = is.numeric(n) && length(n) == 1 && n >= 1,
    "n_iter must be a positive integer" = is.numeric(n_iter) && length(n_iter) == 1 &&
      n_iter >= 1 && n_iter == floor(n_iter)
  )
  n      <- as.integer(n)
  n_iter <- as.integer(n_iter)

  stopifnot(
    "pi must be a numeric vector with values in (0, 1)" =
      is.numeric(pi) && all(pi > 0) && all(pi < 1),
    "h2 must be a numeric vector with values in (0, 1)" =
      is.numeric(h2) && all(h2 > 0) && all(h2 < 1)
  )

  model               <- match.arg(model, choices = c("sparse", "sparse_inf"))
  inf_model           <- match.arg(inf_model, choices = c("beatrice", "susie_inf"))
  effect_distribution <- match.arg(effect_distribution, choices = c("normal", "equal"))

  if (model == "sparse_inf") {
    stopifnot(
      "p_causal must be a single number in (0, 1]" =
        is.numeric(p_causal) && length(p_causal) == 1 &&
        p_causal > 0 && p_causal <= 1
    )
  }

  stopifnot(
    "effect_variance must be a positive number" =
      is.numeric(effect_variance) && length(effect_variance) == 1 &&
      effect_variance > 0,
    "min_maf must be a single number in [0, 0.5]" =
      is.numeric(min_maf) && length(min_maf) == 1 &&
      min_maf >= 0 && min_maf <= 0.5
  )

  # Annotation type
  user_annotation_matrix <- NULL
  if (is.matrix(annotations)) {
    user_annotation_matrix <- annotations
    annotation_type <- "user_supplied"
  } else if (is.character(annotations)) {
    annotation_type <- match.arg(annotations, choices = c("none", "binary", "continuous"))
  } else {
    stop("annotations must be 'none', 'binary', 'continuous', or a matrix.", call. = FALSE)
  }

  if (annotation_type %in% c("binary", "continuous")) {
    stopifnot(
      "n_annotations must be a positive integer" =
        is.numeric(n_annotations) && length(n_annotations) == 1 &&
        n_annotations >= 1 && n_annotations == floor(n_annotations)
    )
    n_annotations <- as.integer(n_annotations)
  }

  stopifnot(
    "coverage must be a single number in (0, 1]" =
      is.numeric(coverage) && length(coverage) == 1 &&
      coverage > 0 && coverage <= 1
  )

  if (annotation_type == "binary" && !is.null(annotation_proportions)) {
    annotation_proportions <- gwfm_validate_annotation_param(
      annotation_proportions, n_annotations, "annotation_proportions"
    )
    stopifnot(
      "annotation_proportions must be in (0, 1)" =
        all(annotation_proportions > 0 & annotation_proportions < 1)
    )
  }

  if (annotation_type != "none" && !is.null(enrichment)) {
    enrichment <- gwfm_validate_annotation_param(enrichment, n_annotations, "enrichment")
    stopifnot("enrichment values must be > 0" = all(enrichment > 0))
  }

  # ---------------------------------------------------------------------------
  # Load region specifications
  # ---------------------------------------------------------------------------

  region_df <- gwfm_load_regions(regions)

  # Subsample regions if coverage < 1
  if (coverage < 1) {
    region_df <- gwfm_subsample_regions(region_df, coverage)
    if (verbose) {
      message(sprintf(
        "coverage = %.2f: using %d of %d available regions (stratified by chromosome).",
        coverage, nrow(region_df),
        nrow(gwfm_load_regions(regions))
      ))
    }
  }

  n_regions <- nrow(region_df)

  # Map region IDs to VCF files
  vcf_files <- gwfm_resolve_vcf_files(region_df, vcf_dir)

  # Resolve p per region
  if (length(p) == 1) {
    p_vec <- rep(as.integer(p), n_regions)
  } else if (length(p) == n_regions) {
    p_vec <- as.integer(p)
  } else {
    stop(
      "p must be a scalar or a vector of length n_regions (", n_regions, ").",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Set master seed
  # ---------------------------------------------------------------------------

  if (!is.null(seed)) set.seed(seed)

  # ---------------------------------------------------------------------------
  # Build parameter grid
  # ---------------------------------------------------------------------------

  if (model == "sparse") {
    param_grid <- expand.grid(
      pi     = pi,
      h2     = h2,
      iter   = seq_len(n_iter),
      stringsAsFactors = FALSE
    )
    param_grid$p_causal <- NA_real_
  } else {
    param_grid <- expand.grid(
      pi      = pi,
      h2      = h2,
      p_causal = p_causal,
      iter    = seq_len(n_iter),
      stringsAsFactors = FALSE
    )
  }

  n_scenarios <- nrow(param_grid)

  if (verbose) {
    if (model == "sparse") {
      message(sprintf(
        "GWFM simulation grid: %d pi x %d h2 x %d iter = %d scenarios",
        length(pi), length(h2), n_iter, n_scenarios
      ))
    } else {
      message(sprintf(
        "GWFM simulation grid: %d pi x %d h2 x %d p_causal x %d iter = %d scenarios",
        length(pi), length(h2), length(p_causal), n_iter, n_scenarios
      ))
    }
    message(sprintf("Regions: %d (across chromosomes %s)",
                    n_regions,
                    paste(sort(unique(region_df$chrom)), collapse = ", ")))
  }

  # ---------------------------------------------------------------------------
  # Simulate genotypes once (shared across all scenarios)
  # ---------------------------------------------------------------------------

  if (verbose) message("\n--- Simulating genotypes ---")

  genotypes <- vector("list", n_regions)

  for (i in seq_len(n_regions)) {
    reg <- region_df[i, ]
    if (verbose) {
      message(sprintf(
        "  Region %d/%d: %s (chr%s:%d-%d) ~%d SNPs",
        i, n_regions, reg$region_id, reg$chrom, reg$start, reg$end, p_vec[i]
      ))
    }

    genotypes[[i]] <- simulate_single_region(
      vcf_file        = vcf_files[i],
      n               = n,
      p               = p_vec[i],
      min_maf         = min_maf,
      max_maf         = NA,
      standardise     = TRUE,
      genetic_map_dir = genetic_map_dir,
      region_id       = i,
      verbose         = verbose
    )

    # Attach region metadata
    genotypes[[i]]$region_label <- reg$region_id
    genotypes[[i]]$chrom        <- reg$chrom
    genotypes[[i]]$start        <- reg$start
    genotypes[[i]]$end          <- reg$end
  }

  # Pre-compute LD matrices (shared across all scenarios)
  if (verbose) message("\n  Pre-computing LD matrices...")
  for (i in seq_len(n_regions)) {
    genotypes[[i]]$LD <- cor(genotypes[[i]]$X)
  }

  # Genome-wide variant count
  p_total <- sum(vapply(genotypes, function(g) g$p, integer(1)))

  if (verbose) {
    message(sprintf("  %d regions, %d total variants.", n_regions, p_total))
  }

  # ---------------------------------------------------------------------------
  # Simulate per-region annotations (shared across all scenarios)
  # ---------------------------------------------------------------------------

  if (annotation_type != "none" && annotation_type != "user_supplied") {
    if (verbose) message("  Simulating annotation matrices...")

    for (i in seq_len(n_regions)) {
      p_i <- genotypes[[i]]$p
      annot_result <- gwfm_simulate_annotations(
        p                      = p_i,
        annotation_type        = annotation_type,
        n_annotations          = n_annotations,
        annotation_proportions = annotation_proportions,
        user_annotation_matrix = NULL
      )
      genotypes[[i]]$annotations_matrix      <- annot_result$matrix
      genotypes[[i]]$annotation_proportions  <- annot_result$proportions
    }
  } else if (annotation_type == "user_supplied") {
    # User supplied a single p_total x m matrix — split by region
    if (nrow(user_annotation_matrix) != p_total) {
      stop(sprintf(
        "User-supplied annotation matrix has %d rows but total variants = %d.",
        nrow(user_annotation_matrix), p_total
      ), call. = FALSE)
    }
    row_offset <- 0L
    for (i in seq_len(n_regions)) {
      p_i <- genotypes[[i]]$p
      genotypes[[i]]$annotations_matrix <- user_annotation_matrix[
        seq(row_offset + 1L, row_offset + p_i), , drop = FALSE
      ]
      row_offset <- row_offset + p_i
    }
  }

  # ---------------------------------------------------------------------------
  # Phenotype simulation loop
  # ---------------------------------------------------------------------------

  if (verbose) message("\n--- Simulating phenotypes ---")

  scenarios <- vector("list", n_scenarios)

  for (sc in seq_len(n_scenarios)) {
    pi_sc      <- param_grid$pi[sc]
    h2_sc      <- param_grid$h2[sc]
    p_causal_sc <- if (model == "sparse_inf") param_grid$p_causal[sc] else NA_real_
    iter_sc    <- param_grid$iter[sc]

    if (verbose) {
      if (model == "sparse") {
        message(sprintf(
          "  Scenario %d/%d: pi=%.1e, h2=%.2f, iter=%d",
          sc, n_scenarios, pi_sc, h2_sc, iter_sc
        ))
      } else {
        message(sprintf(
          "  Scenario %d/%d: pi=%.1e, h2=%.2f, p_causal=%.2f, iter=%d",
          sc, n_scenarios, pi_sc, h2_sc, p_causal_sc, iter_sc
        ))
      }
    }

    # --- Assign causal variants genome-wide -----------------------------------

    causal_assignment <- gwfm_assign_causal_variants(
      genotypes              = genotypes,
      pi                     = pi_sc,
      annotation_type        = annotation_type,
      enrichment             = enrichment,
      n_annotations          = n_annotations
    )

    S_total <- sum(vapply(causal_assignment, function(ca) length(ca$causal_indices), integer(1)))

    if (verbose) {
      message(sprintf("    Causal variants assigned: %d (expected: %.1f)",
                      S_total, pi_sc * p_total))
    }

    # --- Simulate joint phenotype ---------------------------------------------

    pheno_result <- gwfm_simulate_joint_phenotype(
      genotypes           = genotypes,
      causal_assignment   = causal_assignment,
      h2                  = h2_sc,
      model               = model,
      p_causal            = if (model == "sparse_inf") p_causal_sc else NULL,
      inf_model           = inf_model,
      effect_distribution = effect_distribution,
      effect_variance     = effect_variance,
      p_total             = p_total
    )

    y_shared   <- pheno_result$y
    beta_scaled <- pheno_result$beta_scaled  # list of length n_regions

    if (verbose) {
      message(sprintf("    Realised h2: %.4f", pheno_result$h2_realised))
    }

    # --- Compute per-region summary stats from shared y ----------------------

    region_results <- vector("list", n_regions)

    for (i in seq_len(n_regions)) {
      X_i          <- genotypes[[i]]$X
      beta_i       <- beta_scaled[[i]]
      causal_idx_i <- causal_assignment[[i]]$causal_indices

      # Marginal summary stats from the SHARED phenotype
      sumstats_i <- compute_summary_statistics(X = X_i, y = y_shared)

      # Region-level realized PVE: Var(X_i beta_i) / Var(y)
      g_i <- as.numeric(X_i %*% beta_i)
      pve_region <- var(g_i) / var(y_shared)

      region_results[[i]] <- list(
        z        = sumstats_i$z,
        beta_hat = sumstats_i$beta_hat,
        se       = sumstats_i$se,
        y        = y_shared,   # same y for all regions within a scenario
        truth = list(
          causal_indices = causal_idx_i,
          causal_effects = beta_i[causal_idx_i],
          beta_true      = beta_i,
          pve_region     = pve_region,
          S_realized     = length(causal_idx_i),
          pi             = pi_sc,
          h2             = h2_sc,
          model          = model,
          p_causal       = if (model == "sparse_inf") p_causal_sc else NULL,
          enrichment     = causal_assignment[[i]]$enrichment_used,
          annotation_type = annotation_type
        )
      )
    }

    scenarios[[sc]] <- list(
      scenario_id = sc,
      pi          = pi_sc,
      h2          = h2_sc,
      p_causal    = if (model == "sparse_inf") p_causal_sc else NULL,
      iter        = iter_sc,
      model       = model,
      S_total     = S_total,
      h2_realised = pheno_result$h2_realised,
      regions     = region_results
    )
  }

  # ---------------------------------------------------------------------------
  # Package result
  # ---------------------------------------------------------------------------

  params <- list(
    n                      = n,
    n_iter                 = n_iter,
    pi                     = pi,
    h2                     = h2,
    model                  = model,
    p_causal               = if (model == "sparse_inf") p_causal else NULL,
    inf_model              = if (model == "sparse_inf") inf_model else NULL,
    effect_distribution    = effect_distribution,
    effect_variance        = effect_variance,
    annotation_type        = annotation_type,
    n_annotations          = if (annotation_type %in% c("binary", "continuous")) n_annotations else NULL,
    annotation_proportions = annotation_proportions,
    enrichment             = enrichment,
    coverage               = coverage,
    n_regions              = n_regions,
    p_total                = p_total,
    regions                = region_df,
    min_maf                = min_maf,
    vcf_dir                = vcf_dir,
    seed                   = seed,
    simulation_type        = "gwfm"
  )

  result <- list(
    genotypes = genotypes,
    scenarios = scenarios,
    params    = params
  )

  # ---------------------------------------------------------------------------
  # Save (optional)
  # ---------------------------------------------------------------------------

  if (save) {
    stopifnot(
      "output_dir must be a single character string" =
        is.character(output_dir) && length(output_dir) == 1L
    )
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    seed_tag <- if (!is.null(seed)) paste0("seed", seed) else "noseed"
    pi_tag   <- paste(formatC(pi, format = "e", digits = 0), collapse = "-")
    h2_tag   <- paste(h2, collapse = "-")
    fname    <- sprintf("gwfm_%s_pi%s_h2%s_iter%d_%s.rds",
                        model, pi_tag, h2_tag, n_iter, seed_tag)
    fpath    <- file.path(output_dir, fname)
    saveRDS(result, file = fpath)
    if (verbose) message(sprintf("\nGWFM simulation saved to: %s", fpath))
  }

  if (verbose) {
    message(sprintf(
      "\nDone. %d scenarios, %d regions, %d total variants.",
      n_scenarios, n_regions, p_total
    ))
  }

  result
}


# =============================================================================
# Internal: stratified subsampling of regions by chromosome
# =============================================================================

# For each chromosome, samples round(coverage * n_chr) regions without
# replacement. Chromosomes where this rounds to zero are omitted entirely.
# The resulting region set preserves genome-wide spread at all coverage levels.

gwfm_subsample_regions <- function(region_df, coverage) {

  chrs <- sort(unique(region_df$chrom))
  sampled_list <- vector("list", length(chrs))

  for (k in seq_along(chrs)) {
    chr      <- chrs[k]
    chr_df   <- region_df[region_df$chrom == chr, ]
    n_chr    <- nrow(chr_df)
    n_keep   <- round(coverage * n_chr)

    if (n_keep == 0L) next

    n_keep <- min(n_keep, n_chr)
    sel    <- sort(sample(n_chr, n_keep))
    sampled_list[[k]] <- chr_df[sel, ]
  }

  result <- do.call(rbind, sampled_list[!vapply(sampled_list, is.null, logical(1))])

  if (is.null(result) || nrow(result) == 0) {
    stop(
      "coverage = ", coverage, " resulted in zero regions being selected. ",
      "Increase coverage or use a region set with more regions per chromosome.",
      call. = FALSE
    )
  }

  rownames(result) <- NULL
  result
}


# =============================================================================
# Internal: load region specification
# =============================================================================

gwfm_load_regions <- function(regions) {

  if (is.data.frame(regions)) {
    required <- c("region_id", "chrom", "start", "end")
    missing  <- setdiff(required, names(regions))
    if (length(missing) > 0) {
      stop(
        "regions data frame must have columns: ",
        paste(required, collapse = ", "),
        ". Missing: ", paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    return(regions)
  }

  if (!is.character(regions) || length(regions) != 1) {
    stop(
      "regions must be 'representative', a chromosome name (e.g. '1' or 'chr1'), ",
      "or a data frame with columns region_id, chrom, start, end.",
      call. = FALSE
    )
  }

  # Load the bundled gwfm_regions.csv
  bundled <- system.file("data/gwfm_regions.csv", package = "finemapbenchmark")
  if (!file.exists(bundled)) {
    # Fallback for development (running from project root)
    bundled <- "data/gwfm_regions.csv"
  }
  if (!file.exists(bundled)) {
    stop(
      "Cannot find data/gwfm_regions.csv. ",
      "Ensure the file exists in the project data/ directory.",
      call. = FALSE
    )
  }

  all_regions <- read.csv(bundled, stringsAsFactors = FALSE)

  if (regions == "representative") {
    return(all_regions)
  }

  # Chromosome filter: accept "chr1" or "1"
  chrom_query <- sub("^chr", "", regions)
  if (!grepl("^[0-9]+$", chrom_query)) {
    stop(
      "Unrecognised regions specification: '", regions, "'. ",
      "Use 'representative', a chromosome (e.g. '1' or 'chr1'), ",
      "or a data frame.",
      call. = FALSE
    )
  }

  subset_df <- all_regions[as.character(all_regions$chrom) == chrom_query, ]
  if (nrow(subset_df) == 0) {
    stop(
      "No regions found for chromosome ", chrom_query,
      " in data/gwfm_regions.csv.",
      call. = FALSE
    )
  }

  subset_df
}


# =============================================================================
# Internal: resolve VCF file paths for each region
# =============================================================================

gwfm_resolve_vcf_files <- function(region_df, vcf_dir) {

  if (!dir.exists(vcf_dir)) {
    stop(
      "vcf_dir does not exist: ", vcf_dir, "\n",
      "Run scripts/prepare_gwfm_vcfs.R first to download the VCF files.",
      call. = FALSE
    )
  }

  vcf_paths <- file.path(vcf_dir, paste0(region_df$region_id, ".vcf.gz"))
  missing   <- !file.exists(vcf_paths)

  if (any(missing)) {
    stop(
      "VCF file(s) missing for regions:\n",
      paste("  ", region_df$region_id[missing], vcf_paths[missing], sep = "  "),
      "\nRun scripts/prepare_gwfm_vcfs.R to download them.",
      call. = FALSE
    )
  }

  vcf_paths
}


# =============================================================================
# Internal: simulate annotation matrix for one region
# =============================================================================

gwfm_simulate_annotations <- function(p,
                                       annotation_type,
                                       n_annotations,
                                       annotation_proportions,
                                       user_annotation_matrix) {

  if (annotation_type == "none") {
    return(list(matrix = NULL, proportions = NULL))
  }

  if (annotation_type == "user_supplied") {
    return(list(matrix = user_annotation_matrix, proportions = NULL))
  }

  if (annotation_type == "binary") {
    props <- if (is.null(annotation_proportions)) {
      runif(n_annotations, min = 0.01, max = 0.30)
    } else {
      annotation_proportions
    }
    A <- matrix(0, nrow = p, ncol = n_annotations)
    for (k in seq_len(n_annotations)) {
      A[, k] <- rbinom(p, size = 1, prob = props[k])
    }
    colnames(A) <- paste0("annot_", seq_len(n_annotations))
    return(list(matrix = A, proportions = props))
  }

  if (annotation_type == "continuous") {
    A <- matrix(rnorm(p * n_annotations), nrow = p, ncol = n_annotations)
    colnames(A) <- paste0("annot_", seq_len(n_annotations))
    return(list(matrix = A, proportions = NULL))
  }
}


# =============================================================================
# Internal: assign causal variants genome-wide
# =============================================================================

# For each variant j, the causal probability is:
#   pi_j = pi * w_j / mean(w_j)
# where w_j = exp(sum_k A_jk * log(enrichment_k)) [1 if no annotations].
# Each variant is then independently Bernoulli(pi_j).

gwfm_assign_causal_variants <- function(genotypes,
                                         pi,
                                         annotation_type,
                                         enrichment,
                                         n_annotations) {

  n_regions <- length(genotypes)

  # Determine enrichment values (once, shared across regions)
  enrichment_used <- NULL
  if (annotation_type %in% c("binary", "continuous", "user_supplied") &&
      !is.null(genotypes[[1]]$annotations_matrix)) {

    m <- ncol(genotypes[[1]]$annotations_matrix)

    enrichment_used <- if (is.null(enrichment)) {
      runif(m, min = 2, max = 10)
    } else {
      enrichment
    }
  }

  # Per-region causal assignment
  causal_assignment <- vector("list", n_regions)

  for (i in seq_len(n_regions)) {
    p_i <- genotypes[[i]]$p
    A_i <- genotypes[[i]]$annotations_matrix

    if (!is.null(A_i) && !is.null(enrichment_used)) {
      # Compute per-variant weights
      log_enrich  <- log(enrichment_used)
      log_weights <- as.numeric(A_i %*% log_enrich)
      weights     <- exp(log_weights - max(log_weights))
      # Rescale pi so mean(pi_j) = pi
      pi_vec <- pi * weights / mean(weights)
      # Clamp to (0, 1) for safety
      pi_vec <- pmin(pmax(pi_vec, 1e-10), 1 - 1e-10)
    } else {
      pi_vec <- rep(pi, p_i)
    }

    is_causal       <- as.logical(rbinom(p_i, size = 1, prob = pi_vec))
    causal_indices  <- which(is_causal)

    causal_assignment[[i]] <- list(
      causal_indices  = causal_indices,
      enrichment_used = enrichment_used
    )
  }

  causal_assignment
}


# =============================================================================
# Internal: simulate joint phenotype from shared y = sum_i X_i beta_i + eps
# =============================================================================

gwfm_simulate_joint_phenotype <- function(genotypes,
                                           causal_assignment,
                                           h2,
                                           model,
                                           p_causal,
                                           inf_model,
                                           effect_distribution,
                                           effect_variance,
                                           p_total) {

  n_regions <- length(genotypes)
  n         <- genotypes[[1]]$n

  # --- Draw raw effect sizes --------------------------------------------------

  beta_raw <- vector("list", n_regions)

  for (i in seq_len(n_regions)) {
    p_i            <- genotypes[[i]]$p
    causal_idx_i   <- causal_assignment[[i]]$causal_indices
    beta_i         <- rep(0.0, p_i)
    S_i            <- length(causal_idx_i)

    if (S_i > 0) {
      if (effect_distribution == "normal") {
        beta_i[causal_idx_i] <- rnorm(S_i, mean = 0, sd = sqrt(effect_variance))
      } else {
        # "equal": unit magnitudes with random signs
        beta_i[causal_idx_i] <- sample(c(-1, 1), S_i, replace = TRUE)
      }
    }

    beta_raw[[i]] <- beta_i
  }

  # --- Compute total genetic signal -------------------------------------------

  g_sparse_list <- lapply(seq_len(n_regions), function(i) {
    as.numeric(genotypes[[i]]$X %*% beta_raw[[i]])
  })

  g_sparse_total <- Reduce("+", g_sparse_list)
  var_g_sparse   <- var(g_sparse_total)

  # --- Handle degenerate case (no causal variants drawn) ----------------------

  if (var_g_sparse < .Machine$double.eps) {
    warning(
      "No causal variants were drawn (or all effects are zero). ",
      "The phenotype will be pure noise. Try increasing pi or effect_variance.",
      call. = FALSE
    )
    beta_scaled  <- lapply(seq_len(n_regions), function(i) rep(0.0, genotypes[[i]]$p))
    y            <- rnorm(n)
    h2_realised  <- 0.0

    return(list(
      y           = y,
      beta_scaled = beta_scaled,
      h2_realised = h2_realised
    ))
  }

  # --- Scale sparse component to explain its target fraction of h2 ------------

  if (model == "sparse") {
    sparse_target_h2 <- h2
    inf_target_h2    <- 0.0
  } else {
    sparse_target_h2 <- p_causal * h2
    inf_target_h2    <- (1 - p_causal) * h2
  }

  sparse_scale  <- sqrt(sparse_target_h2 / var_g_sparse)
  g_sparse_sc   <- g_sparse_total * sparse_scale

  # Rescale beta accordingly
  beta_scaled <- lapply(beta_raw, function(b) b * sparse_scale)

  # --- Infinitesimal component (sparse_inf only) ------------------------------

  g_inf_sc <- rep(0.0, n)

  if (model == "sparse_inf") {
    g_inf_total <- gwfm_infinitesimal_component(
      genotypes      = genotypes,
      causal_assignment = causal_assignment,
      inf_model      = inf_model,
      p_total        = p_total
    )
    var_g_inf <- var(g_inf_total)

    if (var_g_inf < .Machine$double.eps) {
      warning("Infinitesimal component has zero variance. Skipping inf component.", call. = FALSE)
    } else {
      inf_scale <- sqrt(inf_target_h2 / var_g_inf)
      g_inf_sc  <- g_inf_total * inf_scale
    }
  }

  # --- Total genetic signal and residual noise --------------------------------

  g_total <- g_sparse_sc + g_inf_sc

  # Residual noise to achieve target h2 in expectation
  # Var(y) ≈ Var(g_total) + sigma2 = h2 + (1-h2) = 1 (for standardised y)
  sigma2 <- 1 - h2
  eps    <- rnorm(n, mean = 0, sd = sqrt(sigma2))
  y      <- g_total + eps

  # Realised h2
  h2_realised <- var(g_total) / var(y)

  list(
    y           = y,
    beta_scaled = beta_scaled,
    h2_realised = h2_realised
  )
}


# =============================================================================
# Internal: infinitesimal genetic component for sparse_inf model
# =============================================================================

gwfm_infinitesimal_component <- function(genotypes,
                                          causal_assignment,
                                          inf_model,
                                          p_total) {

  n_regions <- length(genotypes)
  n         <- genotypes[[1]]$n
  g_inf     <- rep(0.0, n)

  if (inf_model == "beatrice") {
    # Infinitesimal effects from NON-causal variants only, across all regions
    for (i in seq_len(n_regions)) {
      p_i           <- genotypes[[i]]$p
      causal_idx_i  <- causal_assignment[[i]]$causal_indices
      noncausal_idx <- setdiff(seq_len(p_i), causal_idx_i)
      m_nc          <- length(noncausal_idx)

      if (m_nc == 0) next

      X_nc   <- genotypes[[i]]$X[, noncausal_idx, drop = FALSE]
      alpha  <- rnorm(m_nc, mean = 0, sd = 1 / sqrt(p_total))
      g_inf  <- g_inf + as.numeric(X_nc %*% alpha)
    }

  } else if (inf_model == "susie_inf") {
    # Infinitesimal effects from ALL variants, across all regions
    for (i in seq_len(n_regions)) {
      p_i   <- genotypes[[i]]$p
      alpha <- rnorm(p_i, mean = 0, sd = 1 / sqrt(p_total))
      g_inf <- g_inf + as.numeric(genotypes[[i]]$X %*% alpha)
    }
  }

  g_inf
}


# =============================================================================
# Internal: validate annotation parameter (scalar or vector of length m)
# =============================================================================

gwfm_validate_annotation_param <- function(x, n_annotations, name) {
  if (!is.numeric(x)) {
    stop(sprintf("%s must be numeric.", name), call. = FALSE)
  }
  if (length(x) == 1) {
    x <- rep(x, n_annotations)
  } else if (length(x) != n_annotations) {
    stop(
      sprintf(
        "%s must be a scalar or vector of length n_annotations (%d). Got length %d.",
        name, n_annotations, length(x)
      ),
      call. = FALSE
    )
  }
  as.numeric(x)
}
