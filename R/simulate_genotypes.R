# =============================================================================
# simulate_genotypes.R
#
# Simulate genotype matrices from 1000 Genomes haplotypes using sim1000G.
#
# Produces one genotype matrix per region, each with realistic LD structure.
# =============================================================================

#' Simulate genotype matrices for multiple independent genomic regions
#'
#' Uses the sim1000G package to simulate genotypes from 1000 Genomes
#' haplotypes. Each region produces an independent genotype matrix with
#' realistic linkage disequilibrium structure.
#'
#' @param n_regions Integer. Number of independent genomic regions to simulate.
#'   Default: 3.
#' @param n Integer. Number of unrelated individuals to simulate. Default: 500.
#' @param p Integer or integer vector. Number of SNPs per region.
#'   If a single integer, the same target is used for all regions.
#'   If a vector, must have length equal to \code{n_regions}, and \code{p[i]}
#'   is the target number of SNPs for region \code{i}. Default: 200.
#'   Note: when using the bundled sim1000G example VCF (i.e. \code{vcf_files = NULL}),
#'   the maximum usable value is approximately 500. Values above this are
#'   automatically capped with a warning. To use larger p, supply your own
#'   VCF files via \code{vcf_files}.
#' @param vcf_files Character vector of paths to VCF files (one per region),
#'   or NULL to use the example VCF bundled with sim1000G. If a single path
#'   is provided and \code{n_regions > 1}, that VCF is reused for all regions
#'   (with independent draws of individuals). Default: NULL.
#' @param min_maf Numeric. Minimum minor allele frequency for filtering
#'   variants from the VCF. Default: 0.01.
#' @param max_maf Numeric or NA. Maximum minor allele frequency. NA means
#'   no upper filter. Default: NA.
#' @param standardise Logical. If TRUE, standardise each genotype column to
#'   mean 0, variance 1 (as in BEATRICE). If FALSE, return raw 0/1/2
#'   genotype coding. Default: TRUE.
#' @param seed Integer or NULL. Random seed for reproducibility. If NULL,
#'   no seed is set. Default: NULL.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list of length \code{n_regions}. Each element is a list with:
#'   \describe{
#'     \item{X}{Genotype matrix (n x p_actual). Standardised if requested.}
#'     \item{X_raw}{Raw 0/1/2 genotype matrix (always included).}
#'     \item{n}{Number of individuals (rows).}
#'     \item{p}{Actual number of SNPs (columns). May differ slightly from
#'       the requested \code{p} depending on available variants in the VCF.}
#'     \item{maf}{Minor allele frequencies of the simulated SNPs.}
#'     \item{variant_ids}{Variant identifiers from the VCF.}
#'     \item{region_id}{Integer index of this region.}
#'     \item{vcf_source}{Path to the VCF file used for this region.}
#'   }
#'
#' @details
#' The sim1000G package simulates genotypes by recombining haplotypes from
#' 1000 Genomes reference panels. This preserves realistic LD patterns.
#'
#' When using the bundled example VCF for multiple regions: each region is
#' simulated as an independent draw of individuals from the same haplotype
#' panel. The LD structure within each region is realistic, and regions are
#' independent of each other (as they would be on different chromosomes).
#'
#' The \code{p} argument controls the \code{maxNumberOfVariants} parameter
#' passed to \code{sim1000G::readVCF}. The actual number of SNPs returned
#' may be slightly less than \code{p} if fewer variants pass the MAF filter.
#'
#' @examples
#' \dontrun{
#' # Default: 3 regions, 500 individuals, 200 SNPs each
#' geno <- simulate_genotypes()
#'
#' # Custom: 5 regions, 1000 individuals, 500 SNPs each
#' geno <- simulate_genotypes(n_regions = 5, n = 1000, p = 500)
#'
#' # Different number of SNPs per region
#' geno <- simulate_genotypes(n_regions = 3, p = c(200, 500, 1000))
#'
#' # User-supplied VCF files
#' geno <- simulate_genotypes(
#'   n_regions = 2,
#'   vcf_files = c("region1.vcf.gz", "region2.vcf.gz"),
#'   p = c(300, 500)
#' )
#' }
#'
#' @export
simulate_genotypes <- function(n_regions = 3,
                               n = 500,
                               p = 200,
                               vcf_files = NULL,
                               min_maf = 0.01,
                               max_maf = NA,
                               standardise = TRUE,
                               seed = NULL,
                               verbose = TRUE) {

  # --- Input validation -------------------------------------------------------

  if (!requireNamespace("sim1000G", quietly = TRUE)) {
    stop(
      "Package 'sim1000G' is required but not installed.\n",
      "Install it with: install.packages('sim1000G')",
      call. = FALSE
    )
  }

  # n_regions
  stopifnot(
    "n_regions must be a positive integer" =
      is.numeric(n_regions) && length(n_regions) == 1 &&
      n_regions == floor(n_regions) && n_regions >= 1
  )
  n_regions <- as.integer(n_regions)

  # n
  stopifnot(
    "n must be a positive integer" =
      is.numeric(n) && length(n) == 1 &&
      n == floor(n) && n >= 1
  )
  n <- as.integer(n)

  # p: scalar or vector of length n_regions
  stopifnot(
    "p must be a positive integer or vector of positive integers" =
      is.numeric(p) && all(p == floor(p)) && all(p >= 1)
  )
  if (length(p) == 1) {
    p <- rep(as.integer(p), n_regions)
  } else if (length(p) != n_regions) {
    stop(
      "If p is a vector, its length must equal n_regions (", n_regions, "). ",
      "Got length ", length(p), ".",
      call. = FALSE
    )
  }
  p <- as.integer(p)

  # Cap p at 500 when using the bundled VCF (which has ~567 variants)
  if (is.null(vcf_files) && any(p > 500)) {
    warning(
      "The bundled sim1000G example VCF contains ~567 variants. ",
      "Capping p at 500 for regions where p > 500. ",
      "To use larger p, supply your own VCF files via the vcf_files argument.",
      call. = FALSE
    )
    p[p > 500] <- 500L
  }

  # vcf_files
  if (!is.null(vcf_files)) {
    if (length(vcf_files) == 1) {
      # Reuse single VCF for all regions
      if (!file.exists(vcf_files)) {
        stop("VCF file not found: ", vcf_files, call. = FALSE)
      }
      vcf_files <- rep(vcf_files, n_regions)
    } else if (length(vcf_files) != n_regions) {
      stop(
        "If vcf_files is provided, its length must be 1 or equal to n_regions (",
        n_regions, "). Got length ", length(vcf_files), ".",
        call. = FALSE
      )
    } else {
      # Check all files exist
      missing <- !file.exists(vcf_files)
      if (any(missing)) {
        stop(
          "VCF file(s) not found:\n",
          paste("  ", vcf_files[missing], collapse = "\n"),
          call. = FALSE
        )
      }
    }
  }

  # MAF
  stopifnot(
    "min_maf must be a single number in [0, 0.5]" =
      is.numeric(min_maf) && length(min_maf) == 1 &&
      min_maf >= 0 && min_maf <= 0.5
  )

  # --- Set seed ---------------------------------------------------------------

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Resolve VCF file paths -------------------------------------------------

  if (is.null(vcf_files)) {
    # Use the bundled example VCF from sim1000G
    examples_dir <- system.file("examples", package = "sim1000G")
    default_vcf <- file.path(examples_dir, "region.vcf.gz")
    if (!file.exists(default_vcf)) {
      stop(
        "Could not find the bundled sim1000G example VCF file.\n",
        "Ensure sim1000G is properly installed.",
        call. = FALSE
      )
    }
    vcf_files <- rep(default_vcf, n_regions)
    if (verbose) {
      message(
        "Using bundled sim1000G example VCF for all ", n_regions, " region(s)."
      )
    }
  }

  # --- Simulate each region ---------------------------------------------------

  regions <- vector("list", n_regions)

  for (i in seq_len(n_regions)) {
    if (verbose) {
      message(
        sprintf(
          "Region %d/%d: simulating %d individuals x ~%d SNPs from %s",
          i, n_regions, n, p[i], basename(vcf_files[i])
        )
      )
    }

    regions[[i]] <- simulate_single_region(
      vcf_file = vcf_files[i],
      n = n,
      p = p[i],
      min_maf = min_maf,
      max_maf = max_maf,
      standardise = standardise,
      region_id = i,
      verbose = verbose
    )
  }

  # --- Return -----------------------------------------------------------------

  if (verbose) {
    total_snps <- sum(vapply(regions, function(r) r$p, integer(1)))
    message(
      sprintf(
        "Done. %d region(s), %d individuals, %d total SNPs.",
        n_regions, n, total_snps
      )
    )
  }

  regions
}


# =============================================================================
# Internal: simulate genotypes for a single region
# =============================================================================

simulate_single_region <- function(vcf_file,
                                   n,
                                   p,
                                   min_maf,
                                   max_maf,
                                   standardise,
                                   region_id,
                                   verbose) {

  # --- Read VCF ---------------------------------------------------------------
  # We read more variants than requested, because the MAF filter and
  # monomorphic-in-simulation filtering may reduce the count.
  # We request up to 2x the target p, then trim to exactly p after simulation.

  buffer_factor <- 1.5
  max_variants_to_read <- as.integer(ceiling(p * buffer_factor))

  # sim1000G uses cat() for progress messages; capture when not verbose.
  # We wrap in a helper that captures stdout (cat) and suppresses messages.
  read_vcf_call <- function() {
    sim1000G::readVCF(
      vcf_file,
      maxNumberOfVariants = max_variants_to_read,
      min_maf = min_maf,
      max_maf = if (is.na(max_maf)) NA else max_maf
    )
  }

  if (!verbose) {
    invisible(utils::capture.output({
      vcf_obj <- read_vcf_call()
    }))
  } else {
    vcf_obj <- read_vcf_call()
  }

  n_variants_available <- length(vcf_obj$maf)

  if (n_variants_available == 0) {
    stop(
      "No variants passed the MAF filter for region ", region_id,
      " (file: ", vcf_file, ").",
      call. = FALSE
    )
  }

  if (n_variants_available < p) {
    warning(
      sprintf(
        "Region %d: requested %d SNPs but only %d variants available after MAF filtering. Using all %d.",
        region_id, p, n_variants_available, n_variants_available
      ),
      call. = FALSE
    )
    target_p <- n_variants_available
  } else {
    target_p <- p
  }

  # --- Subset VCF to target number of SNPs ------------------------------------
  # If we have more variants than needed, randomly sample target_p of them
  # to give the user the requested number. This also introduces variation
  # between regions when using the same VCF.

  if (n_variants_available > target_p) {
    selected_indices <- sort(sample(n_variants_available, target_p))
    vcf_obj <- sim1000G::subsetVCF(vcf_obj, var_index = selected_indices)
  }

  # --- Initialise simulation and generate individuals -------------------------

  # totalNumberOfIndividuals must be >= n; give some headroom
  total_capacity <- as.integer(n * 1.2) + 10

  # sim1000G uses cat() and global state; capture output to keep things tidy
  if (!verbose) {
    invisible(utils::capture.output({
      sim1000G::generateUniformGeneticMap()
      sim1000G::startSimulation(vcf_obj, totalNumberOfIndividuals = total_capacity)
      ids <- sim1000G::generateUnrelatedIndividuals(n)
      genotype_raw <- sim1000G::retrieveGenotypes(ids)
    }))
  } else {
    sim1000G::generateUniformGeneticMap()
    sim1000G::startSimulation(vcf_obj, totalNumberOfIndividuals = total_capacity)
    ids <- sim1000G::generateUnrelatedIndividuals(n)
    genotype_raw <- sim1000G::retrieveGenotypes(ids)
  }

  # genotype_raw is n x p_actual (0/1/2 coding)
  # Rows = individuals, columns = variants

  # --- Post-processing: remove monomorphic SNPs in simulated data -------------

  col_means <- colMeans(genotype_raw)
  sim_maf <- pmin(col_means / 2, 1 - col_means / 2)
  polymorphic <- sim_maf > 0
  if (sum(polymorphic) < ncol(genotype_raw)) {
    n_removed <- sum(!polymorphic)
    if (verbose) {
      message(
        sprintf(
          "  Removed %d monomorphic SNP(s) in simulated data.", n_removed
        )
      )
    }
    genotype_raw <- genotype_raw[, polymorphic, drop = FALSE]
    sim_maf <- sim_maf[polymorphic]
  }

  p_actual <- ncol(genotype_raw)

  # --- Standardise if requested -----------------------------------------------

  if (standardise) {
    X <- scale_genotypes(genotype_raw)
  } else {
    X <- genotype_raw
  }

  # --- Collect variant metadata -----------------------------------------------

  # Try to get variant IDs from the VCF object
  variant_ids <- if (!is.null(vcf_obj$varid)) {
    # Subset to polymorphic if we removed any
    if (length(vcf_obj$varid) == length(polymorphic)) {
      vcf_obj$varid[polymorphic]
    } else {
      paste0("SNP_", seq_len(p_actual))
    }
  } else {
    paste0("SNP_", seq_len(p_actual))
  }

  # --- Return -----------------------------------------------------------------

  list(
    X = X,
    X_raw = genotype_raw,
    n = nrow(X),
    p = p_actual,
    maf = sim_maf,
    variant_ids = variant_ids,
    region_id = region_id,
    vcf_source = vcf_file
  )
}


# =============================================================================
# Internal: standardise genotype matrix to mean 0, variance 1 per column
# =============================================================================

#' Standardise a genotype matrix
#'
#' Centres each column to mean 0 and scales to variance 1.
#' Columns with zero variance (monomorphic) are set to all zeros.
#'
#' @param X_raw Numeric matrix of raw genotypes (0/1/2 coding).
#' @return Numeric matrix of the same dimensions, standardised.
#' @keywords internal
scale_genotypes <- function(X_raw) {
  col_means <- colMeans(X_raw)
  col_sds <- apply(X_raw, 2, sd)

  # Protect against zero-variance columns (should already be removed,
  # but be defensive)
  zero_var <- col_sds == 0
  col_sds[zero_var] <- 1  # avoid division by zero; column will be set to 0

  X <- sweep(X_raw, 2, col_means, "-")
  X <- sweep(X, 2, col_sds, "/")

  # Zero out any columns that had zero variance
  if (any(zero_var)) {
    X[, zero_var] <- 0
  }

  X
}
