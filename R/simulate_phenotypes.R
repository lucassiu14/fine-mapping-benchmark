# =============================================================================
# simulate_phenotypes.R
#
# Simulate phenotypes, compute summary statistics, and record ground truth
# for fine-mapping benchmarking.
#
# Takes the output of simulate_genotypes() and adds:
#   - Phenotype vector y
#   - Summary statistics (z-scores, beta-hat, se)
#   - LD matrix
#   - Ground truth (causal indices, effect sizes, etc.)
#   - Annotations (if applicable)
# =============================================================================


#' Simulate phenotypes and summary statistics for fine-mapping benchmarking
#'
#' Takes the output of \code{\link{simulate_genotypes}} and, for each region,
#' selects causal variants, generates a phenotype vector, computes marginal
#' summary statistics, and records the ground truth.
#'
#' @param genotypes List. Output from \code{\link{simulate_genotypes}}.
#' @param S Integer or integer vector. Number of causal variants per region.
#'   If scalar, the same value is used for all regions. If a vector, must
#'   have length equal to the number of regions. Default: 1.
#' @param phi Numeric or numeric vector. Proportion of variance explained
#'   (PVE) by the genetic component. If scalar, same for all regions.
#'   Must be in (0, 1). Default: 0.1.
#' @param model Character. Genetic architecture model. Either \code{"sparse"}
#'   (standard sparse model, y = Xb + e) or \code{"sparse_inf"} (sparse +
#'   infinitesimal effects). Default: "sparse".
#' @param p_causal Numeric. Proportion of genetic variance attributable to the
#'   sparse (causal) component. Only used when \code{model = "sparse_inf"}.
#'   Must be in (0, 1]. Default: 0.5.
#' @param inf_model Character. Which infinitesimal formulation to use. Either
#'   \code{"beatrice"} (noncausal variants only) or \code{"susie_inf"} (all
#'   variants). Only used when \code{model = "sparse_inf"}. Default: "beatrice".
#' @param effect_distribution Character. Distribution for causal effect sizes.
#'   \code{"normal"} draws from N(0, effect_variance).
#'   \code{"equal"} partitions variance equally among causal variants.
#'   Default: "normal".
#' @param effect_variance Numeric. Variance of the normal effect size
#'   distribution. Only used when \code{effect_distribution = "normal"}.
#'   Default: 0.36 (i.e. sd = 0.6, following SuSiE).
#' @param annotations Character or matrix. Controls functional annotation
#'   simulation. \code{"none"} for no annotations (causal variants selected
#'   uniformly at random). \code{"binary"} for synthetic binary annotations.
#'   \code{"continuous"} for synthetic continuous annotations (drawn from
#'   N(0,1)). Or a user-supplied matrix (p x m) of annotation values.
#'   Default: "none".
#' @param n_annotations Integer. Number of annotation categories. Only used
#'   when \code{annotations} is \code{"binary"} or \code{"continuous"}.
#'   Default: 3.
#' @param annotation_proportions Numeric, scalar, vector, or NULL. Controls
#'   the proportion of variants with value 1 in each binary annotation.
#'   If NULL, proportions are drawn randomly from Uniform(0.01, 0.30) for
#'   each annotation. If scalar, that proportion is used for all annotations.
#'   If a vector, must have length \code{n_annotations}. Only used when
#'   \code{annotations = "binary"}. Default: NULL.
#' @param enrichment Numeric, scalar, vector, or NULL. Fold-enrichment of
#'   each annotation for causal variant selection. If NULL, enrichments are
#'   drawn from Uniform(2, 10) for each annotation. If scalar, that
#'   enrichment is used for all annotations. If a vector, must have length
#'   \code{n_annotations}. Values must be > 0 (values < 1 indicate
#'   depletion). Default: NULL.
#' @param seed Integer or NULL. Random seed. Default: NULL.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return The input \code{genotypes} list, with additional fields appended
#'   to each region:
#'   \describe{
#'     \item{y}{Phenotype vector (n x 1).}
#'     \item{z}{Marginal z-scores (p x 1).}
#'     \item{beta_hat}{Marginal effect size estimates (p x 1).}
#'     \item{se}{Standard errors of marginal effects (p x 1).}
#'     \item{LD}{LD (correlation) matrix (p x p).}
#'     \item{annotations_matrix}{Annotation matrix (p x m), or NULL if none.}
#'     \item{truth}{List containing ground truth:
#'       \describe{
#'         \item{causal_indices}{Integer vector of causal variant indices.}
#'         \item{causal_effects}{Numeric vector of true effect sizes.}
#'         \item{beta_true}{Full p-length vector of true effects (0 for non-causal).}
#'         \item{pve}{Realised PVE.}
#'         \item{S}{Number of causal variants.}
#'         \item{phi}{Target PVE.}
#'         \item{model}{Genetic architecture model used.}
#'         \item{effect_distribution}{Effect size distribution used.}
#'         \item{annotation_type}{Type of annotations used.}
#'         \item{enrichment}{Enrichment values used (or NULL).}
#'         \item{annotation_proportions}{Proportions used for binary annotations (or NULL).}
#'       }
#'     }
#'   }
#'
#' @examples
#' \dontrun{
#' geno <- simulate_genotypes(n_regions = 2, n = 200, p = 100, seed = 1)
#'
#' # Sparse model, no annotations
#' sim <- simulate_phenotypes(geno, S = 3, phi = 0.2, seed = 1)
#'
#' # Sparse + infinitesimal, with binary annotations
#' sim <- simulate_phenotypes(
#'   geno, S = 2, phi = 0.3,
#'   model = "sparse_inf", p_causal = 0.5,
#'   annotations = "binary", n_annotations = 3,
#'   seed = 1
#' )
#' }
#'
#' @export
simulate_phenotypes <- function(genotypes,
                                S = 1,
                                phi = 0.1,
                                model = "sparse",
                                p_causal = 0.5,
                                inf_model = "beatrice",
                                effect_distribution = "normal",
                                effect_variance = 0.36,
                                annotations = "none",
                                n_annotations = 3,
                                annotation_proportions = NULL,
                                enrichment = NULL,
                                seed = NULL,
                                verbose = TRUE) {

  # --- Input validation -------------------------------------------------------

  n_regions <- length(genotypes)

  # S: scalar or vector
  S <- validate_per_region_param(S, n_regions, "S", integer_valued = TRUE)

  # phi: scalar or vector
  phi <- validate_per_region_param(phi, n_regions, "phi")
  stopifnot(
    "phi must be in (0, 1)" = all(phi > 0 & phi < 1)
  )

  # model
  model <- match.arg(model, choices = c("sparse", "sparse_inf"))

  # p_causal
  if (model == "sparse_inf") {
    stopifnot(
      "p_causal must be a single number in (0, 1]" =
        is.numeric(p_causal) && length(p_causal) == 1 &&
        p_causal > 0 && p_causal <= 1
    )
  }

  # inf_model
  inf_model <- match.arg(inf_model, choices = c("beatrice", "susie_inf"))

  # effect_distribution
  effect_distribution <- match.arg(effect_distribution, choices = c("normal", "equal"))

  # effect_variance
  stopifnot(
    "effect_variance must be a positive number" =
      is.numeric(effect_variance) && length(effect_variance) == 1 &&
      effect_variance > 0
  )

  # annotations
  user_annotation_matrix <- NULL
  if (is.matrix(annotations)) {
    user_annotation_matrix <- annotations
    annotation_type <- "user_supplied"
  } else if (is.character(annotations)) {
    annotation_type <- match.arg(annotations, choices = c("none", "binary", "continuous"))
  } else {
    stop("annotations must be 'none', 'binary', 'continuous', or a matrix.", call. = FALSE)
  }

  # n_annotations
  if (annotation_type %in% c("binary", "continuous")) {
    stopifnot(
      "n_annotations must be a positive integer" =
        is.numeric(n_annotations) && length(n_annotations) == 1 &&
        n_annotations == floor(n_annotations) && n_annotations >= 1
    )
    n_annotations <- as.integer(n_annotations)
  }

  # annotation_proportions (binary only)
  if (annotation_type == "binary" && !is.null(annotation_proportions)) {
    annotation_proportions <- validate_annotation_param(
      annotation_proportions, n_annotations, "annotation_proportions"
    )
    stopifnot(
      "annotation_proportions must be in (0, 1)" =
        all(annotation_proportions > 0 & annotation_proportions < 1)
    )
  }

  # enrichment
  if (annotation_type != "none" && !is.null(enrichment)) {
    enrichment <- validate_annotation_param(
      enrichment, n_annotations, "enrichment"
    )
    stopifnot(
      "enrichment must be positive" = all(enrichment > 0)
    )
  }

  # Check S <= p for each region
  for (i in seq_len(n_regions)) {
    if (S[i] > genotypes[[i]]$p) {
      stop(
        sprintf(
          "Region %d: S = %d but only %d SNPs available. S must be <= p.",
          i, S[i], genotypes[[i]]$p
        ),
        call. = FALSE
      )
    }
  }

  # --- Set seed ---------------------------------------------------------------

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Simulate each region ---------------------------------------------------

  for (i in seq_len(n_regions)) {
    if (verbose) {
      message(sprintf(
        "Region %d/%d: S=%d, phi=%.3f, model=%s",
        i, n_regions, S[i], phi[i], model
      ))
    }

    p_i <- genotypes[[i]]$p
    n_i <- genotypes[[i]]$n
    X_i <- genotypes[[i]]$X  # standardised genotype matrix

    # --- Simulate annotations for this region ---------------------------------

    annot_result <- simulate_annotations_for_region(
      p = p_i,
      annotation_type = annotation_type,
      n_annotations = if (annotation_type %in% c("binary", "continuous")) n_annotations else 0,
      annotation_proportions = annotation_proportions,
      user_annotation_matrix = user_annotation_matrix
    )

    A_i <- annot_result$matrix           # p x m annotation matrix, or NULL
    props_i <- annot_result$proportions   # realised proportions, or NULL

    # --- Select causal variants -----------------------------------------------

    causal_result <- select_causal_variants(
      p = p_i,
      S = S[i],
      annotation_matrix = A_i,
      enrichment = enrichment,
      n_annotations = if (annotation_type %in% c("binary", "continuous")) n_annotations else 0,
      annotation_type = annotation_type
    )

    causal_indices <- causal_result$causal_indices
    enrichment_used <- causal_result$enrichment

    # --- Draw effect sizes ----------------------------------------------------

    beta_true <- rep(0, p_i)

    if (effect_distribution == "normal") {
      beta_true[causal_indices] <- rnorm(S[i], mean = 0, sd = sqrt(effect_variance))
    } else if (effect_distribution == "equal") {
      # Effect sizes set later during variance calibration
      beta_true[causal_indices] <- 1  # placeholder, will be rescaled
    }

    # --- Generate phenotype ---------------------------------------------------

    if (model == "sparse") {
      pheno_result <- generate_phenotype_sparse(
        X = X_i,
        beta_true = beta_true,
        phi = phi[i],
        effect_distribution = effect_distribution,
        causal_indices = causal_indices
      )
    } else if (model == "sparse_inf") {
      pheno_result <- generate_phenotype_sparse_inf(
        X = X_i,
        beta_true = beta_true,
        phi = phi[i],
        p_causal = p_causal,
        inf_model = inf_model,
        effect_distribution = effect_distribution,
        causal_indices = causal_indices
      )
    }

    y_i <- pheno_result$y
    beta_true <- pheno_result$beta_true  # may have been rescaled for "equal"

    # --- Compute summary statistics -------------------------------------------

    sumstats <- compute_summary_statistics(X = X_i, y = y_i)

    # --- Store results --------------------------------------------------------

    genotypes[[i]]$y <- y_i
    genotypes[[i]]$z <- sumstats$z
    genotypes[[i]]$beta_hat <- sumstats$beta_hat
    genotypes[[i]]$se <- sumstats$se
    genotypes[[i]]$LD <- sumstats$LD
    genotypes[[i]]$annotations_matrix <- A_i

    genotypes[[i]]$truth <- list(
      causal_indices = causal_indices,
      causal_effects = beta_true[causal_indices],
      beta_true = beta_true,
      pve = pheno_result$pve_realised,
      S = S[i],
      phi = phi[i],
      model = model,
      p_causal = if (model == "sparse_inf") p_causal else NULL,
      inf_model = if (model == "sparse_inf") inf_model else NULL,
      effect_distribution = effect_distribution,
      effect_variance = effect_variance,
      annotation_type = annotation_type,
      enrichment = enrichment_used,
      annotation_proportions = props_i
    )

    if (verbose) {
      message(sprintf(
        "  Causal indices: {%s}, realised PVE: %.4f",
        paste(causal_indices, collapse = ", "),
        pheno_result$pve_realised
      ))
    }
  }

  if (verbose) {
    message("Phenotype simulation complete.")
  }

  genotypes
}


# =============================================================================
# Internal: validate a parameter that can be scalar or per-region vector
# =============================================================================

validate_per_region_param <- function(x, n_regions, name, integer_valued = FALSE) {
  if (!is.numeric(x)) {
    stop(sprintf("%s must be numeric.", name), call. = FALSE)
  }
  if (integer_valued && any(x != floor(x))) {
    stop(sprintf("%s must be integer-valued.", name), call. = FALSE)
  }
  if (length(x) == 1) {
    x <- rep(x, n_regions)
  } else if (length(x) != n_regions) {
    stop(
      sprintf(
        "If %s is a vector, it must have length n_regions (%d). Got length %d.",
        name, n_regions, length(x)
      ),
      call. = FALSE
    )
  }
  if (integer_valued) as.integer(x) else as.numeric(x)
}


# =============================================================================
# Internal: validate annotation parameter (scalar or vector of length m)
# =============================================================================

validate_annotation_param <- function(x, n_annotations, name) {
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


# =============================================================================
# Internal: simulate annotations for a single region
# =============================================================================

simulate_annotations_for_region <- function(p,
                                            annotation_type,
                                            n_annotations,
                                            annotation_proportions,
                                            user_annotation_matrix) {

  if (annotation_type == "none") {
    return(list(matrix = NULL, proportions = NULL))
  }

  if (annotation_type == "user_supplied") {
    if (nrow(user_annotation_matrix) != p) {
      stop(
        sprintf(
          "User-supplied annotation matrix has %d rows but region has %d SNPs.",
          nrow(user_annotation_matrix), p
        ),
        call. = FALSE
      )
    }
    return(list(matrix = user_annotation_matrix, proportions = NULL))
  }

  if (annotation_type == "binary") {
    # Determine proportions
    if (is.null(annotation_proportions)) {
      props <- runif(n_annotations, min = 0.01, max = 0.30)
    } else {
      props <- annotation_proportions
    }

    # Generate binary matrix
    A <- matrix(0, nrow = p, ncol = n_annotations)
    for (k in seq_len(n_annotations)) {
      A[, k] <- rbinom(p, size = 1, prob = props[k])
    }
    colnames(A) <- paste0("annot_", seq_len(n_annotations))

    return(list(matrix = A, proportions = props))
  }

  if (annotation_type == "continuous") {
    # Draw from N(0, 1) following Funmap
    A <- matrix(rnorm(p * n_annotations), nrow = p, ncol = n_annotations)
    colnames(A) <- paste0("annot_", seq_len(n_annotations))

    return(list(matrix = A, proportions = NULL))
  }
}


# =============================================================================
# Internal: select causal variants
# =============================================================================

select_causal_variants <- function(p,
                                   S,
                                   annotation_matrix,
                                   enrichment,
                                   n_annotations,
                                   annotation_type) {

  # --- No annotations: uniform random selection -------------------------------

  if (annotation_type == "none" || is.null(annotation_matrix)) {
    causal_indices <- sort(sample(p, S))
    return(list(causal_indices = causal_indices, enrichment = NULL))
  }

  # --- With annotations: weighted selection -----------------------------------

  m <- ncol(annotation_matrix)

  # Determine enrichment values
  if (is.null(enrichment)) {
    enrichment_vals <- runif(m, min = 2, max = 10)
  } else {
    enrichment_vals <- enrichment
  }

  # Compute unnormalised weights: w_j = exp(sum_k A_jk * log(enrichment_k))
  log_enrichment <- log(enrichment_vals)
  log_weights <- as.numeric(annotation_matrix %*% log_enrichment)
  weights <- exp(log_weights - max(log_weights))  # subtract max for numerical stability

  # Normalise to probabilities
  probs <- weights / sum(weights)

  # Sample S causal indices without replacement
  causal_indices <- sort(sample(p, S, replace = FALSE, prob = probs))

  return(list(causal_indices = causal_indices, enrichment = enrichment_vals))
}


# =============================================================================
# Internal: generate phenotype under sparse model
# =============================================================================

generate_phenotype_sparse <- function(X, beta_true, phi, effect_distribution,
                                      causal_indices) {

  n <- nrow(X)

  # Compute genetic signal
  g_sparse <- as.numeric(X %*% beta_true)
  var_g <- var(g_sparse)

  # Handle "equal" effect distribution: rescale so each causal variant
  # contributes equally to genetic variance
  if (effect_distribution == "equal") {
    S <- length(causal_indices)
    # Set effects so that Var(Xb) gives a reasonable signal,
    # then calibrate via sigma^2
    # Start with unit effects, compute variance, then proceed
    if (var_g == 0) {
      # All effects are 1 (placeholder), but Var(Xb) = 0 is unlikely
      # with standardised genotypes. If it happens, redraw.
      warning("Var(Xb) = 0 with equal effects. Trying random signs.", call. = FALSE)
      beta_true[causal_indices] <- sample(c(-1, 1), S, replace = TRUE)
      g_sparse <- as.numeric(X %*% beta_true)
      var_g <- var(g_sparse)
    }
  }

  # Check for degenerate case
  if (var_g < .Machine$double.eps) {
    stop(
      "Var(Xb) is effectively zero. The drawn effect sizes produced no signal. ",
      "This can happen with very small S or very small effect_variance. ",
      "Try increasing effect_variance or S.",
      call. = FALSE
    )
  }

  # Calibrate residual variance to achieve target PVE
  # phi = Var(Xb) / (Var(Xb) + sigma^2)
  # => sigma^2 = Var(Xb) * (1 - phi) / phi
  sigma2 <- var_g * (1 - phi) / phi

  # Draw phenotype
  e <- rnorm(n, mean = 0, sd = sqrt(sigma2))
  y <- g_sparse + e

  # Compute realised PVE
  pve_realised <- var_g / (var_g + sigma2)

  list(
    y = y,
    beta_true = beta_true,
    pve_realised = pve_realised
  )
}


# =============================================================================
# Internal: generate phenotype under sparse + infinitesimal model
# =============================================================================

generate_phenotype_sparse_inf <- function(X, beta_true, phi, p_causal,
                                          inf_model, effect_distribution,
                                          causal_indices) {

  n <- nrow(X)
  p <- ncol(X)
  S <- length(causal_indices)

  # --- Sparse component -------------------------------------------------------

  g_sparse <- as.numeric(X %*% beta_true)
  var_g_sparse <- var(g_sparse)

  if (effect_distribution == "equal" && var_g_sparse < .Machine$double.eps) {
    warning("Var(Xb) = 0 with equal effects. Trying random signs.", call. = FALSE)
    beta_true[causal_indices] <- sample(c(-1, 1), S, replace = TRUE)
    g_sparse <- as.numeric(X %*% beta_true)
    var_g_sparse <- var(g_sparse)
  }

  if (var_g_sparse < .Machine$double.eps) {
    stop(
      "Var(Xb) is effectively zero for the sparse component. ",
      "Try increasing effect_variance or S.",
      call. = FALSE
    )
  }

  # --- Infinitesimal component ------------------------------------------------

  if (inf_model == "beatrice") {
    # BEATRICE: infinitesimal effects from noncausal variants only
    noncausal_indices <- setdiff(seq_len(p), causal_indices)
    m_nc <- length(noncausal_indices)
    X_nc <- X[, noncausal_indices, drop = FALSE]

    # g_NC ~ N(0, (1/(m-d)) * X_NC * X_NC^T)
    # Equivalent to: alpha_nc ~ N(0, 1/m_nc * I), g_NC = X_NC * alpha_nc
    alpha_nc <- rnorm(m_nc, mean = 0, sd = 1 / sqrt(m_nc))
    g_inf <- as.numeric(X_nc %*% alpha_nc)

  } else if (inf_model == "susie_inf") {
    # SuSiE-inf: infinitesimal effects from all variants
    # alpha_j ~ N(0, tau^2) for all j
    # We draw alpha with unit variance, then rescale during normalisation
    alpha_all <- rnorm(p, mean = 0, sd = 1 / sqrt(p))
    g_inf <- as.numeric(X %*% alpha_all)
  }

  var_g_inf <- var(g_inf)

  if (var_g_inf < .Machine$double.eps) {
    stop(
      "Var(g_inf) is effectively zero for the infinitesimal component. ",
      "This is unexpected with standardised genotypes.",
      call. = FALSE
    )
  }

  # --- Variance normalisation (BEATRICE equation 21) --------------------------
  # Total genetic variance = phi (PVE)
  # Sparse component explains p_causal * phi of total phenotypic variance
  # Infinitesimal component explains (1 - p_causal) * phi
  # Residual explains 1 - phi

  # Scale sparse component
  g_sparse_scaled <- g_sparse * sqrt(p_causal * phi / var_g_sparse)

  # Scale infinitesimal component
  g_inf_scaled <- g_inf * sqrt((1 - p_causal) * phi / var_g_inf)

  # Residual noise
  sigma2 <- 1 - phi
  e <- rnorm(n, mean = 0, sd = sqrt(sigma2))

  # Phenotype
  y <- g_sparse_scaled + g_inf_scaled + e

  # The true beta_true needs to be rescaled consistently
  scale_factor <- sqrt(p_causal * phi / var_g_sparse)
  beta_true <- beta_true * scale_factor

  # Compute realised PVE
  var_total <- var(y)
  pve_realised <- (var(g_sparse_scaled) + var(g_inf_scaled)) / var_total

  list(
    y = y,
    beta_true = beta_true,
    pve_realised = pve_realised
  )
}


# =============================================================================
# Internal: compute marginal summary statistics and LD matrix
# =============================================================================

compute_summary_statistics <- function(X, y) {

  n <- nrow(X)
  p <- ncol(X)

  beta_hat <- numeric(p)
  se <- numeric(p)

  # Marginal OLS regression of y on each column of X (with implicit intercept)
  # Since X is standardised (mean 0) and we centre y, the intercept is 0
  # and beta_hat_j = (x_j^T y) / (x_j^T x_j)
  y_centered <- y - mean(y)

  for (j in seq_len(p)) {
    xj <- X[, j]
    xtx <- sum(xj^2)
    xty <- sum(xj * y_centered)
    beta_hat[j] <- xty / xtx

    # Residuals and standard error (df = n - 2 for intercept + slope)
    resid <- y_centered - xj * beta_hat[j]
    sigma2_j <- sum(resid^2) / (n - 2)
    se[j] <- sqrt(sigma2_j / xtx)
  }

  # z-scores
  z <- beta_hat / se

  # LD matrix: correlation matrix of standardised genotypes
  # For standardised X (mean 0, var 1), this equals (1/(n-1)) * X^T X
  # We use cor() to be safe regardless of exact standardisation
  LD <- cor(X)

  list(
    beta_hat = beta_hat,
    se = se,
    z = z,
    LD = LD
  )
}
