# =============================================================================
# wrapper_sbayesrc.R
#
# SBayesRC (Zheng et al., Nat Genet 56:767, 2024), implemented from the paper
# alone: the Methods ("Summary-data-based low-rank model", "SBayesRC") and
# Supplementary Notes 1, 3, 5, 6, 7, 8 and 10. No code from the SBayesRC
# software (GPL-3) is used; this package stays MIT (user decision,
# 2026-09-16). Where the paper is specific the implementation follows it
# exactly; every point the paper leaves open is listed below with the choice
# made here.
#
# Mapping to this benchmark
# -------------------------
# SBayesRC fits all SNPs of one GWAS jointly, with a block-diagonal LD matrix
# over quasi-independent LD blocks. Here each region of a scenario is one
# block. Note that the simulator gives every region its own phenotype, so a
# scenario is not a single GWAS; results describe SBayesRC applied to that
# design.
#
# The model (paper)
# -----------------
#   * Marginal effects on the standardised scale, b_j = s_j b*_j with
#     s_j = sqrt(1 / (N_j sigma_j^2 + b*_j^2)), assuming unit phenotypic
#     variance (Supp. Note 6).
#   * For each block, R = U Lambda U'. Keep the q leading eigenpairs whose
#     eigenvalues explain at least a proportion rho of the sum of the nonzero
#     eigenvalues, and fit the low-rank model w = Q beta + epsilon with
#     w = Lambda_q^-1/2 U_q' b, Q = Lambda_q^1/2 U_q' and
#     Var(epsilon) = I sigma_e^2 / N (Supp. Note 1).
#   * beta_j ~ sum_k pi_jk N(0, gamma_k sigma_g^2), five components with
#     gamma = [0, 0.001, 0.01, 0.1, 1]% (Methods).
#   * Stick-breaking membership, p_jk = Pr(delta_j >= k | delta_j >= k - 1)
#     for k = 2..5, with pi_j1 = 1 - p_j2, pi_j2 = (1 - p_j3) p_j2, ...,
#     pi_j5 = p_j5 p_j4 p_j3 p_j2 (Supp. Note 5).
#   * Probit link p_jk = Phi(mu_k + A_j' alpha_k); flat prior on mu_k;
#     alpha_kc ~ N(0, sigma_ak^2); sigma_ak^2 ~ scaled-inv-chi2(4, 1)
#     (Methods). Binary annotations are 0/1; quantitative annotations are
#     standardised to mean 0 and variance 1.
#   * sigma_g^2 = sum over blocks of w_hat' w_hat with w_hat = Q beta, computed
#     in every iteration (Supp. Note 7, Algorithm line 23).
#   * sigma_e^2 has a scaled-inv-chi2(nu_e, tau_e^2) prior and is sampled for
#     each block (Supp. Note 8, Algorithm line 25).
#
# The sampler (Supp. Note 8 and its Algorithm)
# --------------------------------------------
#   For each iteration: for each SNP in turn, sample delta_j from its full
#   conditional (beta_j integrated out) and then beta_j given delta_j, and
#   update the residual; then, for k = 2..5 and the SNPs with
#   delta_j >= k - 1, sample the latent l_jk from a normal truncated by
#   z_jk = 1(delta_j >= k) (Albert and Chib), then mu_k and each alpha_kc by
#   single-site Gibbs, then sigma_ak^2; then sigma_g^2; then sigma_e^2 per
#   block. The Supplementary Note's printed full conditionals for beta_j and
#   sigma_e^2 drop the factor N that its stated joint distribution carries
#   (and set Q_j'Q_j to one); here both are derived from that joint
#   distribution:
#     precision_jk = N Q_j'Q_j / sigma_e^2 + 1 / (gamma_k sigma_g^2),
#     mean_jk      = (N / sigma_e^2) r_j / precision_jk,
#     r_j          = Q_j' (w - sum_{j' != j} Q_j' beta_j'),
#     sigma_e^2    ~ scaled-inv-chi2(q + nu_e, (N e'e + nu_e tau_e^2)/(q + nu_e)).
#   The note prints the latent-variable coding both ways round; here
#   z_jk = 1 when the SNP reaches component k, with l_jk > 0, which is the
#   coding its truncated-normal full conditional uses.
#   MCMC: 3,000 iterations, the first 1,000 discarded (Methods). The PIP is
#   the posterior probability of a non-zero component.
#
# Tuning rho (Supp. Note 10)
# --------------------------
#   Pseudo summary statistics b_t = b + sqrt(1/n_t - 1/n) U Lambda^1/2 xi,
#   xi ~ N(0, I), using the eigenpairs kept at rho = 0.995, and
#   b_v = (b n - b_t n_t) / n_v. For rho in (0.995, 0.99, 0.95, 0.9), a
#   150-iteration run without annotations on b_t gives the posterior mean of
#   beta over its last 50 iterations, and the pseudo-validation correlation
#   R = beta' b_v / sqrt(m Var(beta)). rho moves from 0.995 when R > 0 and
#   |R / R_0.995| > 1.25; if the best R is at rho = 0.9 the user is asked to
#   extend the grid.
#
# Points the paper leaves open, and the choice made here
# -------------------------------------------------------
#   * Starting values ("Initialize model parameters"): beta = 0, delta = 1,
#     sigma_g^2 = 0.5, sigma_e^2 = 1, alpha = 0, sigma_ak^2 = 1, and mu_k set
#     so that the starting mixture proportions are
#     (0.990, 0.005, 0.003, 0.001, 0.001), the software's documented defaults.
#   * The sigma_e^2 prior (nu_e, tau_e^2) is not given; nu_e = 4 and
#     tau_e^2 = 1 are used, centring it on the unit phenotypic variance the
#     paper assumes.
#   * "Nonzero" eigenvalues are those above the usual numerical-rank
#     tolerance, max(lambda) * m * machine epsilon.
#   * The training share of the pseudo split is not given; 90% is used.
#   * When several rho values pass the 1.25 rule, the one with the largest R
#     is taken. When R at 0.995 is not positive, the rho with the largest
#     positive R is taken. "Prompt the user" is implemented as an error when
#     the chosen rho is the smallest in the grid.
#   * If every beta is zero in an iteration, sigma_g^2 keeps its previous
#     value (the mixture variances would otherwise all be zero).
#   * The flat prior on mu_k is taken as uniform on [-8, 8] on the probit
#     scale (probabilities down to about 1e-15). On the whole line it gives
#     an improper posterior whenever the SNPs conditioning component k are
#     completely separated (all or none reach k), which is common for the
#     few SNPs in the upper components of a small scenario, and the chain
#     then diverges. mu_k is drawn from the resulting truncated normal.
#   * A component k whose conditioning set (delta_j >= k - 1) is empty leaves
#     mu_k and alpha_k at their current values that iteration.
#   * SBayesRC reports no credible sets, so none are returned.
#
# This file provides:
#   - setup_sbayesrc()                 : dependency check
#   - sbayesrc()                       : the joint fit across regions (blocks)
#   - run_sbayesrc()                   : SBayesRC on a single region
#   - run_sbayesrc_region()            : adapter called by run_methods()
#   - run_sbayesrc_scenario_setup()    : runs sbayesrc() once per scenario
# =============================================================================


# =============================================================================
# setup_sbayesrc()
# =============================================================================

#' Check that the dependencies of sbayesrc are available
#'
#' The implementation uses base R only.
#'
#' @return Invisible TRUE.
#' @export
setup_sbayesrc <- function() {
  invisible(TRUE)
}


# =============================================================================
# sbayesrc(): the joint fit
# =============================================================================

#' SBayesRC across a set of LD blocks
#'
#' Fits SBayesRC (see the file header) with each region as one LD block.
#'
#' @param z_list List of z-score vectors, one per region.
#' @param ld_list List of LD correlation matrices, one per region.
#' @param n Integer. GWAS sample size.
#' @param annot_list List of annotation matrices (variants x annotations) with
#'   the same columns in every region, or NULL to fit without annotations.
#' @param beta_hat_list,se_list Optional lists of marginal effects and their
#'   standard errors. When NULL, \code{b = z / sqrt(n + z^2)}, which is the
#'   same scaling.
#' @param rho Numeric or NULL. Eigenvalue cut-off. NULL (default) tunes it by
#'   pseudo-validation over \code{tune_grid}.
#' @param tune_grid Numeric. Cut-offs to tune over. Default
#'   \code{c(0.995, 0.99, 0.95, 0.9)}.
#' @param tune_iter,tune_keep Integer. Iterations of each tuning run, and the
#'   final iterations averaged. Defaults 150 and 50.
#' @param train_prop Numeric. Training share of the pseudo split. Default 0.9.
#' @param n_iter,burn_in Integer. MCMC length and burn-in. Defaults 3000 and
#'   1000.
#' @param gamma Numeric length 5. Component variance scales (fractions of
#'   sigma_g^2). Default \code{c(0, 1e-5, 1e-4, 1e-3, 1e-2)}.
#' @param start_h2,start_pi Starting sigma_g^2 and mixture proportions.
#' @param nu_alpha,tau2_alpha Prior of the annotation-effect variances.
#' @param nu_e,tau2_e Prior of the residual variances.
#' @param seed Integer or NULL. Random seed.
#'
#' @return A list with \code{pip} and \code{beta} (per-region lists),
#'   \code{rho}, \code{tuning} (data frame or NULL), \code{q} (kept
#'   eigenpairs per region), \code{sigma2_g} (posterior mean), \code{mu} and
#'   \code{alpha} (posterior means).
#' @export
sbayesrc <- function(z_list, ld_list, n, annot_list = NULL,
                     beta_hat_list = NULL, se_list = NULL,
                     rho = NULL, tune_grid = c(0.995, 0.99, 0.95, 0.9),
                     tune_iter = 150L, tune_keep = 50L, train_prop = 0.9,
                     n_iter = 3000L, burn_in = 1000L,
                     gamma = c(0, 1e-5, 1e-4, 1e-3, 1e-2),
                     start_h2 = 0.5,
                     start_pi = c(0.990, 0.005, 0.003, 0.001, 0.001),
                     nu_alpha = 4, tau2_alpha = 1, nu_e = 4, tau2_e = 1,
                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  R <- length(z_list)
  stopifnot(R >= 1L, length(ld_list) == R, length(gamma) == 5L,
            length(start_pi) == 5L, burn_in < n_iter)
  sizes <- vapply(z_list, length, integer(1))

  # --- scaled marginal effects (Supp. Note 6) --------------------------------
  b_list <- lapply(seq_len(R), function(i) {
    if (!is.null(beta_hat_list) && !is.null(se_list)) {
      bh <- as.numeric(beta_hat_list[[i]]); s <- as.numeric(se_list[[i]])
      bh / sqrt(n * s^2 + bh^2)
    } else {
      z <- as.numeric(z_list[[i]])
      z / sqrt(n + z^2)
    }
  })

  # --- annotations: binary kept, quantitative standardised --------------------
  A <- NULL
  if (!is.null(annot_list)) {
    A <- do.call(rbind, lapply(annot_list, as.matrix))
    stopifnot(nrow(A) == sum(sizes))
    for (cc in seq_len(ncol(A))) {
      if (!all(A[, cc] %in% c(0, 1))) {
        s_cc <- stats::sd(A[, cc])
        if (!is.finite(s_cc) || s_cc == 0) {
          stop("annotation ", cc, " is constant", call. = FALSE)
        }
        A[, cc] <- (A[, cc] - mean(A[, cc])) / s_cc
      }
    }
  }

  eig <- lapply(ld_list, .sbrc_eigen)

  # --- tune rho by pseudo-validation (Supp. Note 10) --------------------------
  tuning <- NULL
  if (is.null(rho)) {
    rho_max <- max(tune_grid)
    n_t <- round(train_prop * n); n_v <- n - n_t
    b_t <- vector("list", R); b_v <- vector("list", R)
    for (i in seq_len(R)) {
      e  <- eig[[i]]
      qm <- .sbrc_q(e$values, rho_max)
      xi <- stats::rnorm(qm)
      noise <- e$vectors[, seq_len(qm), drop = FALSE] %*%
        (sqrt(e$values[seq_len(qm)]) * xi)
      b_t[[i]] <- b_list[[i]] + sqrt(1 / n_t - 1 / n) * as.numeric(noise)
      b_v[[i]] <- (b_list[[i]] * n - b_t[[i]] * n_t) / n_v
    }
    bv_all <- unlist(b_v, use.names = FALSE)
    Rcor <- vapply(tune_grid, function(r) {
      blocks <- .sbrc_blocks(eig, b_t, r)
      fit <- .sbrc_mcmc(blocks, A = NULL, n = n_t, n_iter = tune_iter,
                        burn_in = tune_iter - tune_keep, gamma = gamma,
                        start_h2 = start_h2, start_pi = start_pi,
                        nu_alpha = nu_alpha, tau2_alpha = tau2_alpha,
                        nu_e = nu_e, tau2_e = tau2_e)
      bm <- unlist(fit$beta, use.names = FALSE)
      sum(bm * bv_all) / sqrt(length(bm) * stats::var(bm))
    }, numeric(1))
    tuning <- data.frame(rho = tune_grid, R = Rcor)
    rho <- .sbrc_choose_rho(tune_grid, Rcor)
  }

  # --- main run ------------------------------------------------------------------
  blocks <- .sbrc_blocks(eig, b_list, rho)
  fit <- .sbrc_mcmc(blocks, A = A, n = n, n_iter = n_iter, burn_in = burn_in,
                    gamma = gamma, start_h2 = start_h2, start_pi = start_pi,
                    nu_alpha = nu_alpha, tau2_alpha = tau2_alpha,
                    nu_e = nu_e, tau2_e = tau2_e)
  c(fit, list(rho = rho, tuning = tuning,
              q = vapply(blocks, function(bk) nrow(bk$Q), integer(1))))
}


# =============================================================================
# run_sbayesrc(): one region
# =============================================================================

#' Run SBayesRC on a single region
#'
#' Treats the region as the only LD block. Under \code{run_methods()} the
#' fit is joint across the regions of a scenario instead (see
#' \code{run_sbayesrc_scenario_setup()}).
#'
#' @param z Numeric vector. Marginal z-scores.
#' @param LD Matrix. LD correlation matrix.
#' @param n Integer. Sample size.
#' @param annotations Matrix or NULL. Annotation matrix.
#' @param beta_hat,se Numeric or NULL. Marginal effects and standard errors.
#' @param variant_ids Character or NULL.
#' @param ... Passed to \code{\link{sbayesrc}}.
#'
#' @return List with pip, credible_sets (empty), method = "sbayesrc",
#'   params, runtime_seconds and additional.
#' @export
run_sbayesrc <- function(z, LD, n, annotations = NULL, beta_hat = NULL,
                         se = NULL, variant_ids = NULL, ...) {
  t0 <- Sys.time()
  p <- length(z)
  stopifnot("LD must be p x p" = is.matrix(LD) && all(dim(LD) == c(p, p)))
  extra <- list(...)
  extra <- extra[intersect(names(extra), names(formals(sbayesrc)))]
  fit <- do.call(sbayesrc, c(list(
    z_list = list(z), ld_list = list(LD), n = n,
    annot_list = if (is.null(annotations)) NULL else list(annotations),
    beta_hat_list = if (is.null(beta_hat)) NULL else list(beta_hat),
    se_list = if (is.null(se)) NULL else list(se)), extra))
  .sbrc_region_output(fit, 1L, runtime = as.numeric(difftime(Sys.time(), t0,
                                                             units = "secs")),
                      joint = FALSE, variant_ids = variant_ids)
}


# =============================================================================
# run_sbayesrc_region(): adapter for run_methods()
# =============================================================================

#' Region adapter for SBayesRC
#'
#' Returns the region's result from the joint fit computed by
#' \code{run_sbayesrc_scenario_setup()}; without that result it fits the
#' region on its own.
#'
#' @param region_geno,region_pheno Region data.
#' @param .sbayesrc_cache Named list of per-region results keyed by z
#'   fingerprint. Default NULL.
#' @param .sbayesrc_error Character or NULL. Error from the scenario setup.
#' @param ... Passed to \code{\link{run_sbayesrc}}.
#'
#' @return The per-region result list.
#' @export
run_sbayesrc_region <- function(region_geno, region_pheno,
                                .sbayesrc_cache = NULL,
                                .sbayesrc_error = NULL, ...) {
  if (!is.null(.sbayesrc_error)) stop(.sbayesrc_error, call. = FALSE)
  if (!is.null(.sbayesrc_cache)) {
    hit <- .sbayesrc_cache[[.fb_fingerprint(region_pheno$z)]]
    if (!is.null(hit)) return(hit)
  }
  A <- region_geno$annotations_matrix
  if (is.null(A)) A <- region_pheno$annotations_matrix
  run_sbayesrc(z = region_pheno$z, LD = region_geno$LD, n = region_geno$n,
               annotations = A, beta_hat = region_pheno$beta_hat,
               se = region_pheno$se, variant_ids = region_geno$variant_ids, ...)
}


# =============================================================================
# Scenario-level setup: the joint fit across the regions
# =============================================================================

#' Scenario-level SBayesRC fit
#'
#' Fits all regions of a scenario jointly, one LD block per region, and
#' returns the per-region results keyed by z fingerprint.
#'
#' @param genotypes List of per-region genotype data.
#' @param regions List of per-region phenotype data for one scenario.
#' @param user_args Named list of arguments for \code{\link{sbayesrc}}.
#'
#' @return A named list merged into each region's arguments.
#' @export
run_sbayesrc_scenario_setup <- function(genotypes, regions, user_args) {
  R <- length(regions)
  t0 <- Sys.time()
  A_list <- lapply(seq_len(R), function(i) {
    A <- genotypes[[i]]$annotations_matrix
    if (is.null(A)) regions[[i]]$annotations_matrix else A
  })
  has_A <- !vapply(A_list, is.null, logical(1))
  tryCatch({
    if (any(has_A) && !all(has_A)) {
      stop("some regions have no annotation matrix", call. = FALSE)
    }
    known <- names(formals(sbayesrc))
    extra <- user_args[intersect(names(user_args), known)]
    fit <- do.call(sbayesrc, c(list(
      z_list        = lapply(regions, `[[`, "z"),
      ld_list       = lapply(genotypes, `[[`, "LD"),
      n             = genotypes[[1]]$n,
      annot_list    = if (all(has_A)) A_list else NULL,
      beta_hat_list = lapply(regions, `[[`, "beta_hat"),
      se_list       = lapply(regions, `[[`, "se")), extra))
    runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cache <- lapply(seq_len(R), function(i)
      .sbrc_region_output(fit, i, runtime = runtime / R, joint = TRUE,
                          variant_ids = genotypes[[i]]$variant_ids))
    names(cache) <- vapply(regions, function(r) .fb_fingerprint(r$z),
                           character(1))
    list(.sbayesrc_cache = cache)
  }, error = function(e) list(.sbayesrc_error = paste("SBayesRC failed:",
                                                      conditionMessage(e))))
}


# =============================================================================
# Internals
# =============================================================================

.sbrc_region_output <- function(fit, i, runtime, joint, variant_ids) {
  list(
    pip             = fit$pip[[i]],
    credible_sets   = list(),
    method          = "sbayesrc",
    input_type      = "summary",
    params          = list(rho = fit$rho, q = fit$q[[i]], joint = joint),
    runtime_seconds = runtime,
    additional      = list(
      posterior_mean_beta = fit$beta[[i]],
      tuning              = fit$tuning,
      sigma2_g            = fit$sigma2_g,
      mu                  = fit$mu,
      alpha               = fit$alpha,
      variant_ids         = variant_ids
    )
  )
}

# Eigen-decomposition of one LD block with the nonzero eigenvalues marked.
.sbrc_eigen <- function(LD) {
  e <- eigen(as.matrix(LD), symmetric = TRUE)
  tol <- max(e$values) * nrow(LD) * .Machine$double.eps
  nz <- sum(e$values > tol)
  list(values = e$values[seq_len(nz)],
       vectors = e$vectors[, seq_len(nz), drop = FALSE])
}

# Smallest q whose leading eigenvalues explain at least rho of the total.
.sbrc_q <- function(values, rho) {
  cs <- cumsum(values) / sum(values)
  min(which(cs >= rho - 1e-12))
}

# Low-rank blocks at cut-off rho: Q (q x m), w (q), column norms d.
.sbrc_blocks <- function(eig, b_list, rho) {
  lapply(seq_along(eig), function(i) {
    e <- eig[[i]]
    q <- .sbrc_q(e$values, rho)
    Ut <- t(e$vectors[, seq_len(q), drop = FALSE])
    lam <- e$values[seq_len(q)]
    Q <- sqrt(lam) * Ut
    w <- as.numeric((Ut %*% b_list[[i]]) / sqrt(lam))
    list(Q = Q, w = w, d = colSums(Q * Q))
  })
}

# Choice of rho from the pseudo-validation correlations.
.sbrc_choose_rho <- function(grid, Rcor) {
  base <- max(grid)
  R0 <- Rcor[grid == base]
  ok <- is.finite(Rcor) & Rcor > 0
  if (!any(ok)) {
    stop("all pseudo-validation correlations are non-positive; check the ",
         "summary statistics", call. = FALSE)
  }
  chosen <- base
  if (is.finite(R0) && R0 > 0) {
    pass <- ok & abs(Rcor / R0) > 1.25
    if (any(pass)) chosen <- grid[pass][which.max(Rcor[pass])]
  } else {
    chosen <- grid[ok][which.max(Rcor[ok])]
  }
  if (chosen == min(grid)) {
    stop("the best eigen cut-off is the smallest in the tuning grid (",
         chosen, "); extend the grid or check the summary statistics",
         call. = FALSE)
  }
  chosen
}

# Stick-breaking log mixture proportions from the probit linear predictors.
# eta: SNPs x 4 matrix for k = 2..5. Returns SNPs x 5.
.sbrc_logpi <- function(eta) {
  lp <- stats::pnorm(eta, log.p = TRUE)
  lq <- stats::pnorm(eta, lower.tail = FALSE, log.p = TRUE)
  c2 <- lp[, 1L]
  c3 <- c2 + lp[, 2L]
  c4 <- c3 + lp[, 3L]
  cbind(lq[, 1L],
        c2 + lq[, 2L],
        c3 + lq[, 3L],
        c4 + lq[, 4L],
        c4 + lp[, 4L])
}

# Normal(eta, 1) truncated to (0, Inf) where upper is TRUE, (-Inf, 0) otherwise.
.sbrc_rtruncnorm <- function(eta, upper) {
  lu <- log(stats::runif(length(eta)))
  out <- numeric(length(eta))
  i1 <- which(upper); i0 <- which(!upper)
  if (length(i1)) {
    out[i1] <- eta[i1] - stats::qnorm(lu[i1] + stats::pnorm(eta[i1], log.p = TRUE),
                                      log.p = TRUE)
  }
  if (length(i0)) {
    out[i0] <- eta[i0] + stats::qnorm(lu[i0] + stats::pnorm(-eta[i0], log.p = TRUE),
                                      log.p = TRUE)
  }
  out
}

# One draw from N(m, s^2) truncated to (lo, hi), stable in the tails.
.sbrc_rtnorm1 <- function(m, s, lo, hi) {
  a <- (lo - m) / s; b <- (hi - m) / s
  if (a > 0) {
    la <- stats::pnorm(a, lower.tail = FALSE, log.p = TRUE)
    lb <- stats::pnorm(b, lower.tail = FALSE, log.p = TRUE)
    lv <- la + log1p(-stats::runif(1L) * (1 - exp(lb - la)))
    return(m + s * stats::qnorm(lv, lower.tail = FALSE, log.p = TRUE))
  }
  if (b < 0) {
    la <- stats::pnorm(a, log.p = TRUE)
    lb <- stats::pnorm(b, log.p = TRUE)
    lv <- lb + log1p(-stats::runif(1L) * (1 - exp(la - lb)))
    return(m + s * stats::qnorm(lv, log.p = TRUE))
  }
  pa <- stats::pnorm(a); pb <- stats::pnorm(b)
  m + s * stats::qnorm(pa + stats::runif(1L) * (pb - pa))
}

# One Gibbs sweep over the SNPs of a block.
.sbrc_sweep <- function(Q, d, eps, beta, delta, logpi, v, n, s2e) {
  nk <- n / s2e
  vk <- v[-1L]
  lv <- log(vk)
  for (j in seq_along(beta)) {
    qj  <- Q[, j]
    bj  <- beta[j]
    rhs <- sum(qj * eps) + d[j] * bj
    P   <- nk * d[j] + 1 / vk
    num <- nk * rhs
    ll  <- c(0, -0.5 * (lv + log(P)) + 0.5 * num * num / P) + logpi[j, ]
    pr  <- exp(ll - max(ll))
    k   <- 1L + sum(stats::runif(1L) * sum(pr) > cumsum(pr))
    bnew <- if (k == 1L) 0 else stats::rnorm(1L, num / P[k - 1L], sqrt(1 / P[k - 1L]))
    if (bnew != bj) eps <- eps + qj * (bj - bnew)
    beta[j]  <- bnew
    delta[j] <- k
  }
  list(eps = eps, beta = beta, delta = delta)
}

# The SBayesRC MCMC over a list of low-rank blocks.
.sbrc_mcmc <- function(blocks, A, n, n_iter, burn_in, gamma, start_h2,
                       start_pi, nu_alpha, tau2_alpha, nu_e, tau2_e,
                       mu_bound = 8) {
  B <- length(blocks)
  sizes <- vapply(blocks, function(bk) ncol(bk$Q), integer(1))
  m <- sum(sizes)
  starts <- cumsum(c(0L, sizes))
  C <- if (is.null(A)) 0L else ncol(A)

  # starting values
  pi0 <- start_pi / sum(start_pi)
  p0 <- c(1 - pi0[1],
          sum(pi0[3:5]) / sum(pi0[2:5]),
          sum(pi0[4:5]) / sum(pi0[3:5]),
          pi0[5] / sum(pi0[4:5]))
  mu    <- stats::qnorm(p0)
  alpha <- matrix(0, 4L, C)
  s2a   <- rep(tau2_alpha, 4L)
  s2g   <- start_h2
  s2e   <- rep(1, B)
  beta  <- lapply(sizes, numeric)
  delta <- lapply(sizes, function(s) rep(1L, s))
  eps   <- lapply(blocks, function(bk) bk$w)

  eta_all <- function() {
    base <- matrix(mu, m, 4L, byrow = TRUE)
    if (C > 0L) base <- base + A %*% t(alpha)
    base
  }
  logpi <- .sbrc_logpi(eta_all())

  pip_sum <- lapply(sizes, numeric)
  beta_sum <- lapply(sizes, numeric)
  s2g_sum <- 0; mu_sum <- numeric(4L); alpha_sum <- matrix(0, 4L, C)
  n_kept <- 0L

  for (it in seq_len(n_iter)) {
    v <- gamma * s2g
    # SNP effects and memberships, block by block
    for (b in seq_len(B)) {
      rows <- starts[b] + seq_len(sizes[b])
      sw <- .sbrc_sweep(blocks[[b]]$Q, blocks[[b]]$d, eps[[b]], beta[[b]],
                        delta[[b]], logpi[rows, , drop = FALSE], v, n, s2e[b])
      beta[[b]] <- sw$beta; delta[[b]] <- sw$delta
    }
    dall <- unlist(delta, use.names = FALSE)

    # annotation effects (probit, stick-breaking)
    for (k in 2:5) {
      idx <- which(dall >= k - 1L)
      if (length(idx) == 0L) next
      zk <- dall[idx] >= k
      Ak <- if (C > 0L) A[idx, , drop = FALSE] else NULL
      eta <- mu[k - 1L] + (if (C > 0L) as.numeric(Ak %*% alpha[k - 1L, ])
                           else numeric(length(idx)))
      l <- .sbrc_rtruncnorm(eta, zk)
      res <- l - eta
      # mu_k, flat prior
      res <- res + mu[k - 1L]
      mu[k - 1L] <- .sbrc_rtnorm1(mean(res), sqrt(1 / length(idx)),
                                  -mu_bound, mu_bound)
      res <- res - mu[k - 1L]
      if (C > 0L) {
        ss <- colSums(Ak * Ak)
        for (cc in seq_len(C)) {
          ac <- Ak[, cc]
          res <- res + ac * alpha[k - 1L, cc]
          Ckc <- ss[cc] + 1 / s2a[k - 1L]
          alpha[k - 1L, cc] <- stats::rnorm(1L, sum(ac * res) / Ckc, sqrt(1 / Ckc))
          res <- res - ac * alpha[k - 1L, cc]
        }
        s2a[k - 1L] <- (sum(alpha[k - 1L, ]^2) + nu_alpha * tau2_alpha) /
          stats::rchisq(1L, C + nu_alpha)
      }
    }
    logpi <- .sbrc_logpi(eta_all())

    # genetic variance and residual variances
    s2g_new <- 0
    for (b in seq_len(B)) {
      what <- as.numeric(blocks[[b]]$Q %*% beta[[b]])
      eps[[b]] <- blocks[[b]]$w - what
      s2g_new <- s2g_new + sum(what * what)
      q <- length(eps[[b]])
      s2e[b] <- (n * sum(eps[[b]]^2) + nu_e * tau2_e) / stats::rchisq(1L, q + nu_e)
    }
    if (s2g_new > 0) s2g <- s2g_new

    if (it > burn_in) {
      for (b in seq_len(B)) {
        pip_sum[[b]]  <- pip_sum[[b]] + (delta[[b]] >= 2L)
        beta_sum[[b]] <- beta_sum[[b]] + beta[[b]]
      }
      s2g_sum <- s2g_sum + s2g
      mu_sum <- mu_sum + mu
      alpha_sum <- alpha_sum + alpha
      n_kept <- n_kept + 1L
    }
  }

  list(pip      = lapply(pip_sum, function(x) x / n_kept),
       beta     = lapply(beta_sum, function(x) x / n_kept),
       sigma2_g = s2g_sum / n_kept,
       mu       = mu_sum / n_kept,
       alpha    = alpha_sum / n_kept)
}
