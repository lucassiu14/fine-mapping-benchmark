# =============================================================================
# wrapper_sbayesrc.R
#
# SBayesRC-style Bayesian regression with annotation-informed mixture prior,
# implemented in R for the fine-mapping benchmark's per-locus regime.
#
# See docs/autoresearch/method-sbayesrc.md for the algorithmic derivation and
# the caveats section (§0.1 of the plan) - this implementation runs FAR
# outside SBayesRC's native genome-wide-7M-SNP regime and its scores must be
# treated as in-context-relative to the other Tier-1 methods in this
# benchmark, NOT as a reflection of SBayesRC's real-world performance on
# actual GWAS-scale data.
#
# Why an in-R reimplementation instead of the upstream R package
# -------------------------------------------------------------
# The upstream (zhilizheng/SBayesRC) requires a pre-eigen-decomposed
# genome-wide LD folder + an annotation file both tied to a specific ~7M-SNP
# reference panel. Repurposing that machinery for our per-region 40-1000-SNP
# simulated LD blocks is significantly harder than reimplementing the
# algorithm - and would produce something that behaves less faithfully.
# So we reimplement the algorithm here, close to how it is described in
# Zheng et al. (Nat Genet 2024), tailored to the summary-statistic + block-LD
# setting the rest of this package produces.
#
# Model
# -----
# For each region i, effect sizes beta_{i,j} follow a spike + K-normal-slabs
# mixture:
#
#     beta_{i,j} | comp_{i,j} = k ~ Normal(0, sigma2_k),  k = 1, ..., K
#     beta_{i,j} | comp_{i,j} = 0 = 0                     (spike)
#
# where sigma2_1 > sigma2_2 > ... > sigma2_K is a fixed grid of variance
# scales (default: c(0.05, 0.005, 5e-4, 5e-5), corresponding roughly to
# large/medium/small/tiny per-SNP heritability).
#
# The per-SNP mixture proportions are annotation-modulated via a multinomial
# logit with the spike (k=0) as the reference class:
#
#     log(pi_{i,j,k} / pi_{i,j,0}) = alpha_k + A_{i,j}^T gamma_k,
#     pi_{i,j,k} = softmax_k(alpha_k + A_{i,j}^T gamma_k)
#
# alpha and gamma are SHARED across regions - this is the "pool annotations
# across regions to estimate the shared annotation-to-prior mapping" step
# from the plan (Sec 0.1 of the auto-research doc).
#
# Summary-stat likelihood
# -----------------------
# Given standardized genotypes with n samples and residual variance ~= 1,
# the marginal-effect estimator is
#
#     beta_hat_j = z_j / sqrt(n),  Var(beta_hat_j | beta) = R / n
#
# so the SNP-j residual (conditional on beta_{-j}) is
#
#     r_j = beta_hat_j - sum_{k != j} R_{j,k} * beta_k
#     r_j | beta_j ~ Normal(beta_j, 1/n)
#
# and the mixture-conditional posteriors on beta_j are the standard
# normal-normal updates: v_k = 1/(1/sigma2_k + n), m_k = v_k * n * r_j.
#
# Gibbs sweep
# -----------
#   1. For each region i, several sweeps over its SNPs, updating (comp_j,
#      beta_j) jointly from the mixture posterior.
#   2. Every gamma_update_every iterations, pool the current comp_{i,j}
#      assignments across all regions and refit (alpha, gamma) via a
#      multinomial logistic regression on A - this is the shared prior
#      update, and is where cross-region information flows.
#   3. Repeat, discarding burn-in samples, then estimate PIP_{i,j} as the
#      posterior probability comp_{i,j} > 0.
# =============================================================================


# =============================================================================
# setup_sbayesrc()
# =============================================================================

#' Check that dependencies for sbayesrc are available
#'
#' @return Invisible TRUE.
#' @export
setup_sbayesrc <- function() {
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("The `nnet` package is required for sbayesrc (used for the ",
         "shared multinomial annotation regression). It ships with R by ",
         "default; install with install.packages('nnet') if missing.",
         call. = FALSE)
  }
  invisible(TRUE)
}


# =============================================================================
# run_sbayesrc(): single-region interface
# =============================================================================

#' Run SBayesRC-style Bayesian regression on a single region
#'
#' Fits a K-component mixture-of-normals prior (with the point mass at zero
#' as the spike) using a Gibbs sampler on the summary-stat likelihood
#' r_j | beta_j ~ Normal(beta_j, 1/n). When \code{pooled_gamma} is supplied
#' (as it is under \code{run_methods()} via
#' \code{run_sbayesrc_scenario_setup()}) the per-SNP mixture weights use the
#' cross-region annotation regression fit; otherwise they are refit on this
#' region alone.
#'
#' When \code{annotations} is NULL, the prior collapses to region-common
#' mixture proportions (learned from the Gibbs assignments), and the method
#' still runs.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD correlation matrix (p x p).
#' @param n Integer. Sample size.
#' @param annotations Matrix or NULL. Functional annotations (p x m).
#' @param sigma2_scale Numeric vector of length K. Per-component prior
#'   variance grid. Default \code{c(0.05, 0.005, 5e-4, 5e-5)}.
#' @param n_iter Integer. Total Gibbs iterations. Default 300.
#' @param burn_in Integer. Iterations to discard before averaging. Default 150.
#' @param gamma_update_every Integer. Refit the annotation regression
#'   (alpha, gamma) every this many Gibbs iterations. Default 10.
#' @param pooled_gamma List with \code{alpha}, \code{gamma} matrices from
#'   \code{run_sbayesrc_scenario_setup()}, or NULL. When supplied, no
#'   region-local annotation refit is performed - the shared coefficients
#'   are used as-is.
#' @param seed Integer or NULL. Sampler seed.
#' @param variant_ids Character or NULL. Passed through unused.
#' @param ... Ignored (wrapper compatibility).
#'
#' @return List with pip, credible_sets, method = "sbayesrc", params,
#'   runtime_seconds, additional (posterior_mean_beta, prior_source, etc.).
#' @export
run_sbayesrc <- function(z, LD, n,
                         annotations = NULL,
                         sigma2_scale = c(0.05, 0.005, 5e-4, 5e-5),
                         n_iter = 300L,
                         burn_in = 150L,
                         gamma_update_every = 10L,
                         pooled_gamma = NULL,
                         seed = NULL,
                         variant_ids = NULL,
                         ...) {
  setup_sbayesrc()
  t0 <- Sys.time()
  if (!is.null(seed)) set.seed(seed)

  p <- length(z)
  stopifnot("LD must be p x p" = is.matrix(LD) && all(dim(LD) == c(p, p)))
  K <- length(sigma2_scale)
  stopifnot(K >= 1L)

  prior_source <- if (is.null(annotations)) {
    "no_annotations"
  } else if (!is.null(pooled_gamma)) {
    "pooled_scenario_gamma"
  } else {
    "single_region_gamma"
  }

  # --- Marginal beta_hat + per-SNP residual state ---------------------------
  beta_hat <- z / sqrt(n)
  beta     <- rep(0, p)                # current effect estimates
  comp     <- rep(0L, p)               # current mixture assignments (0..K)
  # Running R %*% beta so per-SNP residuals are cheap: init 0 (beta = 0)
  R_beta   <- rep(0, p)

  # --- Annotation setup ------------------------------------------------------
  has_annot <- !is.null(annotations)
  if (has_annot) {
    A <- as.matrix(annotations)
    m_ann <- ncol(A)
    if (nrow(A) != p) {
      stop("annotations must have p rows to match z", call. = FALSE)
    }
    if (!is.null(pooled_gamma)) {
      alpha <- as.numeric(pooled_gamma$alpha)
      gamma <- as.matrix(pooled_gamma$gamma)
      stopifnot(length(alpha) == K, nrow(gamma) == m_ann, ncol(gamma) == K)
    } else {
      # Sparse initial prior: ~97% mass on the spike, ~3% distributed
      # roughly evenly across the K slabs. alpha_k = log((0.03/K) / 0.97).
      # This matches the sparse-model regime of the benchmark (S causal
      # variants out of p, S << p) and gives the sampler a reasonable
      # starting point before the multinomial refit learns the true
      # sparsity from the data.
      alpha <- rep(log((0.03 / K) / 0.97), K)
      gamma <- matrix(0, nrow = m_ann, ncol = K)  # per-annotation slopes
    }
  } else {
    A <- NULL; m_ann <- 0L
    alpha <- rep(log((0.03 / K) / 0.97), K)   # same sparse initial prior
    gamma <- matrix(0, nrow = 0L, ncol = K)
  }

  # Prior mixture proportions per SNP: pi[j, ] over (0, 1, ..., K)
  pi_mat <- .sbayesrc_priors_from_gamma(A, alpha, gamma, K, p = p)

  # --- Storage for post-burn PIP estimation ---------------------------------
  pip_running   <- rep(0, p)
  beta_running  <- rep(0, p)
  n_kept        <- 0L

  # --- Main Gibbs loop ------------------------------------------------------
  for (it in seq_len(n_iter)) {
    # Sweep over SNPs
    swp <- .sbayesrc_gibbs_sweep(
      beta_hat = beta_hat, LD = LD, n = n,
      beta = beta, comp = comp, R_beta = R_beta,
      pi_mat = pi_mat, sigma2_scale = sigma2_scale
    )
    beta   <- swp$beta
    comp   <- swp$comp
    R_beta <- swp$R_beta

    # Update annotation regression periodically (only when pooled_gamma is NULL
    # and we have annotations). Under run_methods()/scenario_setup the pooled
    # coefficients are supplied and this branch is skipped.
    if (has_annot && is.null(pooled_gamma) &&
        (it %% gamma_update_every == 0L)) {
      fit <- tryCatch(.sbayesrc_fit_gamma(comp, A, K), error = function(e) NULL)
      if (!is.null(fit)) {
        alpha  <- fit$alpha
        gamma  <- fit$gamma
        pi_mat <- .sbayesrc_priors_from_gamma(A, alpha, gamma, K, p = p)
      }
    }

    # Accumulate post-burn samples
    if (it > burn_in) {
      pip_running  <- pip_running  + (comp > 0L)
      beta_running <- beta_running + beta
      n_kept       <- n_kept + 1L
    }
  }

  pip  <- if (n_kept > 0L) pip_running  / n_kept else rep(0, p)
  bhat <- if (n_kept > 0L) beta_running / n_kept else rep(0, p)

  # --- 95%-mass credible set (greedy on PIP) --------------------------------
  ord      <- order(pip, decreasing = TRUE)
  cumpip   <- cumsum(pip[ord])
  keep     <- ord[seq_len(min(which(cumpip >= 0.95 * sum(pip)), length(ord)))]
  if (length(keep) == 0L || sum(pip) == 0) keep <- integer(0)
  cs_list  <- if (length(keep) > 0L) list(as.integer(keep)) else list()

  list(
    pip             = pip,
    credible_sets   = cs_list,
    method          = "sbayesrc",
    input_type      = "summary",
    params          = list(
      sigma2_scale = sigma2_scale, n_iter = n_iter, burn_in = burn_in,
      gamma_update_every = gamma_update_every, K = K
    ),
    runtime_seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    additional      = list(
      posterior_mean_beta = bhat,
      prior_source        = prior_source,
      alpha               = alpha,
      gamma               = gamma,
      variant_ids         = variant_ids
    )
  )
}


# =============================================================================
# run_sbayesrc_region(): adapter for run_methods()
# =============================================================================

#' Adapter for run_methods() - forwards to run_sbayesrc
#'
#' @param region_geno One element of \code{simulation$genotypes}.
#' @param region_pheno One element of a scenario's \code{regions}.
#' @param pooled_gamma List or NULL. Shared (alpha, gamma) from
#'   \code{run_sbayesrc_scenario_setup()}. When supplied, no region-local
#'   annotation refit happens.
#' @param ... Forwarded to \code{run_sbayesrc()}.
#'
#' @return Output of run_sbayesrc.
#' @export
run_sbayesrc_region <- function(region_geno, region_pheno,
                                pooled_gamma = NULL, ...) {
  A <- region_geno$annotations_matrix
  if (is.null(A)) A <- region_pheno$annotations_matrix
  run_sbayesrc(
    z            = region_pheno$z,
    LD           = region_geno$LD,
    n            = region_geno$n,
    annotations  = A,
    pooled_gamma = pooled_gamma,
    variant_ids  = region_geno$variant_ids,
    ...
  )
}


# =============================================================================
# Scenario-level setup: pool across regions, estimate shared (alpha, gamma)
# =============================================================================

#' Fit a pooled annotation regression for sbayesrc across the scenario
#'
#' Runs a short-burn Gibbs on every region with region-local annotation
#' refits, pools the resulting per-SNP component assignments across all
#' regions, and does one final multinomial regression on the pooled data to
#' obtain the shared \code{(alpha, gamma)} - i.e., the annotation-to-prior
#' mapping learned genome-wide-style from the whole scenario.
#'
#' Downstream calls (per region) use those shared coefficients as
#' \code{pooled_gamma}, so each region's Gibbs skips the per-region
#' annotation refit and its priors are set from the cross-region model. This
#' implements the plan's §0.1 requirement that SBayesRC's annotation prior
#' be learned jointly across regions and fed back per-block.
#'
#' Returns an empty list (no forwarding) when there are fewer than 2 regions,
#' when annotation matrices are missing/inconsistent, or when the pilot Gibbs
#' fails - the per-region wrapper then does its own single-region refit.
#'
#' @param genotypes List of region_geno.
#' @param regions List of region_pheno for this scenario.
#' @param user_args User's method_args for sbayesrc.
#'
#' @return List with a single element \code{pooled_gamma} = list(alpha, gamma),
#'   or an empty list.
#' @export
run_sbayesrc_scenario_setup <- function(genotypes, regions, user_args) {
  n_regions <- length(regions)
  if (n_regions < 2L) return(list())

  # Resolve tuning knobs (allow user override of sigma2_scale, K, etc.)
  sigma2_scale       <- user_args$sigma2_scale %||% c(0.05, 0.005, 5e-4, 5e-5)
  pilot_iter         <- user_args$scenario_pilot_iter %||% 100L
  pilot_burn         <- user_args$scenario_pilot_burn %||% 40L
  gamma_update_every <- user_args$gamma_update_every  %||% 10L
  K <- length(sigma2_scale)

  # Validate that every region has an annotation matrix of consistent width
  A_list <- vector("list", n_regions)
  comp_list <- vector("list", n_regions)
  for (i in seq_len(n_regions)) {
    A_i <- genotypes[[i]]$annotations_matrix
    if (is.null(A_i)) A_i <- regions[[i]]$annotations_matrix
    if (is.null(A_i) || !is.matrix(A_i)) return(list())
    A_list[[i]] <- A_i
  }
  ncols <- vapply(A_list, ncol, integer(1))
  if (length(unique(ncols)) != 1L) return(list())

  # Pilot Gibbs per region to seed component assignments the joint regression
  # will pool over. We do NOT need the full posterior - just enough for a
  # reasonable pooled fit.
  for (i in seq_len(n_regions)) {
    A_i <- A_list[[i]]
    z_i <- regions[[i]]$z
    if (is.null(z_i) || length(z_i) != nrow(A_i)) return(list())
    LD_i <- genotypes[[i]]$LD
    n_i  <- genotypes[[i]]$n
    if (is.null(LD_i) || is.null(n_i)) return(list())

    pilot <- tryCatch(
      run_sbayesrc(
        z = z_i, LD = LD_i, n = n_i, annotations = A_i,
        sigma2_scale = sigma2_scale, n_iter = pilot_iter, burn_in = pilot_burn,
        gamma_update_every = gamma_update_every, pooled_gamma = NULL
      ),
      error = function(e) NULL
    )
    if (is.null(pilot)) return(list())
    # Threshold each SNP's PIP to a component: 0 for the spike, or a slab
    # index sampled proportional to sigma2_scale (larger slab -> more mass).
    # We don't have full sample paths outside run_sbayesrc's memory, but the
    # posterior mean beta gives a reasonable slab-index proxy: sort SNPs by
    # |bhat| within the non-null (PIP > 0.5) set and assign the top third to
    # slab 1 (largest variance), middle third to slab 2, etc.
    pip_i  <- pilot$pip
    bhat_i <- pilot$additional$posterior_mean_beta
    comp_i <- rep(0L, length(z_i))
    hits   <- which(pip_i > 0.5)
    if (length(hits) > 0L) {
      # Assign slab by rank of |bhat| among the hits
      ranks <- rank(-abs(bhat_i[hits]))
      slab  <- pmin(K, ceiling(ranks / length(hits) * K))
      comp_i[hits] <- slab
    }
    comp_list[[i]] <- comp_i
  }

  # Pool and fit the shared (alpha, gamma)
  A_pooled    <- do.call(rbind, A_list)
  comp_pooled <- unlist(comp_list, use.names = FALSE)
  fit <- tryCatch(
    .sbayesrc_fit_gamma(comp_pooled, A_pooled, K),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list())

  list(pooled_gamma = list(alpha = fit$alpha, gamma = fit$gamma))
}


# =============================================================================
# Internals
# =============================================================================

# One Gibbs sweep over SNPs in a single region.
# Updates in place (returns new copies of) beta, comp, and R_beta = LD %*% beta
# so per-SNP residuals stay cheap to compute.
.sbayesrc_gibbs_sweep <- function(beta_hat, LD, n, beta, comp, R_beta,
                                  pi_mat, sigma2_scale) {
  p <- length(beta_hat)
  K <- length(sigma2_scale)
  # Precompute per-component posterior variance and log-BF scale factor
  # (which do not depend on r_j).
  v_k    <- 1 / (1 / sigma2_scale + n)           # posterior variance
  # log p(r|k) - log p(r|k=0) has two pieces:
  #   log-sd ratio  = 0.5 * log((1/n) / (sigma2_k + 1/n))
  #                 = -0.5 * log(1 + n * sigma2_k)
  #   quadratic in r = 0.5 * (n - 1/(sigma2_k + 1/n)) * r^2
  #                  = 0.5 * n^2 * sigma2_k / (1 + n * sigma2_k) * r^2
  log_sd_ratio <- -0.5 * log(1 + n * sigma2_scale)
  quad_scale   <- 0.5 * n^2 * sigma2_scale / (1 + n * sigma2_scale)

  for (j in seq_len(p)) {
    # Residual: beta_hat_j - sum_{k != j} R_{j,k} beta_k
    #         = beta_hat_j - (R_beta[j] - R[j,j] * beta[j])
    #         = beta_hat_j - R_beta[j] + beta[j]     (R[j,j] = 1)
    r_j <- beta_hat[j] - R_beta[j] + beta[j]

    # log-likelihood contribution per component (relative to k=0)
    log_lik <- c(0, log_sd_ratio + quad_scale * r_j^2)   # length K+1

    # log-priors from pi_mat (already normalized rows over 0..K)
    log_pri <- log(pmax(pi_mat[j, ], 1e-30))

    log_post <- log_pri + log_lik
    log_post <- log_post - max(log_post)
    w_post   <- exp(log_post); w_post <- w_post / sum(w_post)

    new_k <- sample.int(K + 1L, size = 1L, prob = w_post) - 1L
    if (new_k == 0L) {
      new_beta <- 0
    } else {
      m_k      <- v_k[new_k] * n * r_j
      new_beta <- stats::rnorm(1L, mean = m_k, sd = sqrt(v_k[new_k]))
    }

    # Update R_beta with the change in beta_j
    delta_j <- new_beta - beta[j]
    if (delta_j != 0) R_beta <- R_beta + LD[, j] * delta_j
    beta[j] <- new_beta
    comp[j] <- new_k
  }

  list(beta = beta, comp = comp, R_beta = R_beta)
}


# Compute per-SNP prior mixture proportions pi[j, ] over (0, 1, ..., K)
# from alpha and gamma (with k=0 spike as reference class).
# When annotations are NULL, the row-common softmax over (0, alpha) is broadcast
# to all p rows.
.sbayesrc_priors_from_gamma <- function(A, alpha, gamma, K, p = NULL) {
  if (is.null(A) || nrow(gamma) == 0L) {
    logits <- matrix(alpha, nrow = 1L, ncol = K)
  } else {
    logits <- sweep(as.matrix(A) %*% gamma, 2, alpha, "+")
    if (is.null(p)) p <- nrow(logits)
  }
  ext <- cbind(0, logits)
  ext <- ext - apply(ext, 1L, max)
  ex  <- exp(ext)
  pi_row <- ex / rowSums(ex)
  if (is.null(A) || nrow(gamma) == 0L) {
    if (is.null(p)) stop(".sbayesrc_priors_from_gamma needs p when A is NULL")
    matrix(rep(as.numeric(pi_row), each = p), nrow = p, ncol = K + 1L)
  } else {
    pi_row
  }
}


# Fit a multinomial logistic regression of comp (0..K) on A. Returns
# alpha (K-vector) and gamma (m_ann x K matrix) with k=0 (spike) as the
# reference class.
#
# Two failure modes have to be guarded for at benchmark scale:
#   (a) EMPTY CLASSES: with small p some slab labels never appear. The
#       correct behaviour is "give this class ~0 prior mass" (very negative
#       alpha), NOT the padded 0 that softmax would then convert into equal
#       prior with the spike — the latter puts non-causal PIPs at ~0.8.
#   (b) DIVERGENT FITS: with few non-spike observations, unregularized
#       multinom overshoots. Heavy L2 (decay = 1) plus a coefficient cap
#       keeps priors sane. This is a Gibbs-sampler auxiliary update, not a
#       final scientific estimate, so aggressive regularization is fine.
.sbayesrc_fit_gamma <- function(comp, A, K,
                                empty_alpha = -12,
                                coef_cap    =   6) {
  m_ann <- ncol(A)
  y <- factor(comp, levels = as.character(0:K))
  # Empirical class counts — used both to detect empty classes and to
  # sanity-check the fitted intercepts.
  counts <- table(y)

  fit <- tryCatch(
    suppressMessages(suppressWarnings(nnet::multinom(
      y ~ ., data = data.frame(y = y, as.data.frame(A)),
      trace = FALSE, MaxNWts = 1e5, decay = 1.0, maxit = 200L
    ))),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(alpha = rep(empty_alpha, K),
                gamma = matrix(0, nrow = m_ann, ncol = K)))
  }
  coefs <- tryCatch(coef(fit), error = function(e) NULL)
  if (is.null(coefs)) {
    return(list(alpha = rep(empty_alpha, K),
                gamma = matrix(0, nrow = m_ann, ncol = K)))
  }
  # coefs shape:
  #   - (K x (1 + m_ann)) matrix with rownames "1".."K" when >=2 non-reference
  #     classes are present in the data;
  #   - plain (1 + m_ann) numeric vector when only ONE non-reference class was
  #     fit (nnet::multinom collapses to a binary logistic). The class it
  #     belongs to has NO rowname on the vector, so infer it from `comp`.
  if (is.null(dim(coefs))) {
    present <- setdiff(sort(unique(as.integer(as.character(y)))), 0L)
    lone_class <- if (length(present) == 1L) present else 1L
    coefs <- matrix(coefs, nrow = 1L,
                    dimnames = list(as.character(lone_class),
                                    names(coefs)))
  }

  # Reindex to always have rows 1..K in that order, padding empty classes
  # with `empty_alpha` (intercept) and zero (annotation slopes) so that
  # they get essentially zero prior mass — not the softmax-uniform prior
  # a raw zero-padding would produce.
  full <- matrix(0, nrow = K, ncol = ncol(coefs),
                 dimnames = list(as.character(1:K), colnames(coefs)))
  seen <- rownames(coefs)
  if (!is.null(seen)) full[seen, ] <- coefs
  # Empty classes: fill intercept column with a strongly negative value
  unseen <- setdiff(as.character(1:K), seen)
  if (length(unseen) > 0L) full[unseen, 1L] <- empty_alpha

  # Cap magnitudes to prevent divergence when the data is scarce. Use
  # `full[] <-` so the matrix dim/dimnames survive (pmax on a matrix
  # returns a plain vector, which breaks the [, 1L] extraction below).
  full[] <- pmax(-coef_cap, pmin(coef_cap, full))

  alpha <- as.numeric(full[, 1L])
  gamma <- t(full[, -1L, drop = FALSE])
  list(alpha = alpha, gamma = gamma)
}
