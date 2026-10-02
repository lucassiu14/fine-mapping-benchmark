# =============================================================================
# wrapper_polyfun_ldsc.R
#
# PolyFun + SuSiE (Weissbrod et al., Nat Genet 52:1355, 2020), following the
# paper's Online Methods ("PolyFun fine-mapping method") step by step.
#
# Where the paper does not state an implementation detail, the released
# PolyFun code is followed (github.com/omerwe/polyfun, MIT licence: polyfun.py,
# ldsc_polyfun/regressions.py, ldsc_polyfun/jackknife.py, ldsc_polyfun/irwls.py,
# compute_ldscores_from_ld.py, finemapper.py). Every such detail is listed
# below so that the provenance of each step is explicit.
#
# Mapping to this benchmark
# -------------------------
# PolyFun works on whole chromosomes. Here each region of a scenario is
# treated as its own chromosome, numbered 1..R in region order, so the
# paper's odd/even split and the exclusion of the target chromosome apply to
# regions. This needs at least four regions (two per parity, for the
# leave-one-chromosome-out choice of the ridge penalty).
#
# The five steps (paper)
# ----------------------
#   1. Estimate annotation coefficients tau and intercept b by L2-regularised
#      S-LDSC using only SNPs on even (resp. odd) chromosomes, minimising
#        sum_i (chi2_i - n sum_c tau_c l(i,c) - n b - 1)^2 + lambda sum_c tau_c^2,
#      with lambda chosen from a geometric grid of 100 values from 1e-8 to
#      100 by the average out-of-chromosome error within that set.
#   2. Per-SNP heritability for SNPs on odd (resp. even) chromosomes,
#      var[beta_i | a_i] = sum_c tau_c a_ic, using the step-1 coefficients.
#   3. Partition all SNPs into 20 bins of similar per-SNP heritability with
#      Ckmedian.1d.dp.
#   4. For target chromosome r, apply S-LDSC with non-negativity constraints
#      to the 20 bins using SNPs on chromosomes of the same parity as r,
#      excluding r. Floor the per-bin estimates at q * max (q = 1/100) and
#      rescale them to the same genome-wide total.
#   5. Prior causal probability of each SNP proportional to its bin's
#      per-SNP heritability, summing to one within the locus.
# Fine-mapping (paper, simulations): SuSiE with 10 causal SNPs per locus and
# the per-locus causal effect variance (scaled_prior_variance) set by the
# modified HESS estimator; all other SuSiE parameters at their defaults.
#
# Implementation details taken from the released code (paper silent)
# -------------------------------------------------------------------
#   * LD scores l(i,c) = sum_j r2_ij a_jc with the unbiased estimator
#     r2 = r^2 (n-1)/(n-2) - 1/(n-2) (compute_ldscores_from_ld.py), over each
#     region's full LD matrix. n is the LD sample size (the GWAS N for
#     in-sample LD, n_ref for a reference panel).
#   * The regression-weight LD score w_ld is the non-partitioned LD score
#     sum_j r2_ij over the regression SNPs.
#   * M per annotation is the column sum of the annotation over all SNPs.
#   * SNPs with chi2 >= max(0.001 N, 80) are dropped from both regressions.
#   * Regression weights: LDSC Hsq.weights at the aggregate h2 estimate
#     M_tot (mean(chi2) - 1) / mean(N ld_tot), clipped to [0, 1], with
#     ld and w_ld floored at 1; square-rooted and normalised to sum to one
#     (the old_weights path of LD_Score_Regression, IRWLS._weight).
#   * Free intercept as the last column. Step 1 centres the non-intercept
#     columns before weighting, scales every column to unit L2 norm, grids
#     lambda over mean(diag(X'X)) * 10^[-8, 2], leaves the intercept
#     unpenalised, and picks lambda by the sklearn r2_score of the
#     out-of-chromosome predictions (Jackknife_Ridge).
#   * Coefficients are divided back by the column norms and by mean(N).
#   * After Ckmedian, bins are rebalanced so that no bin except the first has
#     fewer than 10 SNPs (PolyFun.create_df_bins).
#   * Step 4 solves the weighted problem, intercept included, by exact
#     non-negative least squares (scipy.optimize.nnls; Lawson-Hanson here).
#   * Modified HESS (finemapper.py): SNPs with p below the 0.005 quantile of
#     the locus p-values; one of each pair with |r| > 0.99 removed by a
#     random maximal independent set; h2 = a' R^-1 a - m/n with a = z/sqrt(n);
#     averaged over 100 random sets; prior variance h2 / L; an error when
#     that is <= 0 or >= 1.
#   * SuSiE is called with bhat = z, shat = 1, R, n, L, the prior weights,
#     the estimated residual variance and max_iter = 100. susieR 0.15 no
#     longer has susie_suff_stat(); susie_rss() takes the same inputs.
#
# Where the paper and the released code differ, the paper is followed
# (user decision, 2026-09-16)
# ----------------------------------------------------------------------
#   (a) Step 2 applies the coefficients fitted on the OPPOSITE parity set.
#       As released, compute_snpvar_chr() takes est_loco_ridge for the set
#       that does not contain r, i.e. the fit that excludes that set, which
#       is the fit on r's own parity including r.
#   (b) lambda is chosen within each parity set. As released, one lambda is
#       chosen on all chromosomes (reestimate_lambda = False).
#
# Choices this setting forces (not in the paper)
# ----------------------------------------------
#   * A base annotation (all ones) is added. PolyFun's baseline-LF model has
#     no base annotation because its MAF bins already sum to one for every
#     SNP; the simulated annotations have no such columns.
#   * Without annotations there is no prior to estimate, and the
#     non-functionally-informed SuSiE of finemapper.py is run (uniform
#     prior, the same SuSiE settings including HESS).
#   * Ties in the per-SNP heritability are ordered by SNP index before
#     binning; the released code sorts with pandas' default (unstable) sort.
#
# This file provides:
#   - setup_polyfun_ldsc()                : dependency check
#   - polyfun_priors()                    : steps 1-5 across the regions
#   - run_polyfun_ldsc()                  : PolyFun + SuSiE on one region,
#                                           given its prior
#   - run_polyfun_ldsc_region()           : adapter called by run_methods()
#   - run_polyfun_ldsc_scenario_setup()   : runs polyfun_priors() once per
#                                           scenario
# =============================================================================


# =============================================================================
# setup_polyfun_ldsc()
# =============================================================================

#' Check that the dependencies of polyfun_ldsc are available
#'
#' PolyFun needs susieR for fine-mapping and Ckmeans.1d.dp for the binning
#' step (the same R package the released PolyFun calls).
#'
#' @param packages Character. Packages to check. Default both.
#'
#' @return Invisible TRUE if all are available.
#' @export
setup_polyfun_ldsc <- function(packages = c("susieR", "Ckmeans.1d.dp")) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(pkg, " is required for polyfun_ldsc but is not installed.\n",
           "Install it with install.packages('", pkg, "').", call. = FALSE)
    }
  }
  invisible(TRUE)
}


# =============================================================================
# polyfun_priors(): steps 1-5 across all regions of a scenario
# =============================================================================

#' Estimate PolyFun prior causal probabilities across regions
#'
#' Runs steps 1-5 of PolyFun (see the file header) with each region treated
#' as its own chromosome, numbered by its position in the input lists.
#'
#' @param z_list List of numeric z-score vectors, one per region.
#' @param ld_list List of LD correlation matrices, one per region.
#' @param annot_list List of annotation matrices (variants x annotations),
#'   one per region, all with the same columns. A base annotation is added.
#' @param n Integer. GWAS sample size.
#' @param n_ld Integer vector or NULL. Sample size behind each LD matrix,
#'   used in the unbiased r^2. NULL uses \code{n} for every region.
#' @param num_bins Integer. Number of bins in step 3. Default 20 (paper).
#' @param q Numeric. Ratio between the largest and the floored per-SNP
#'   heritability in step 4. Default 100 (paper: q = 1/100 of the maximum).
#' @param min_bin_size Integer. Bin rebalancing threshold. Default 10.
#' @param num_lambdas Integer. Size of the ridge grid. Default 100.
#'
#' @return A list with \code{prior} (list of per-region prior probabilities,
#'   each summing to one), \code{snpvar} (per-region per-SNP heritabilities
#'   after step 4), \code{snpvar_ridge} (after step 2), \code{bins}
#'   (per-region bin labels), \code{taus_ridge} (step-1 coefficients by
#'   parity set) and \code{lambda} (the chosen ridge penalties).
#' @export
polyfun_priors <- function(z_list, ld_list, annot_list, n, n_ld = NULL,
                           num_bins = 20L, q = 100, min_bin_size = 10L,
                           num_lambdas = 100L) {
  R <- length(z_list)
  stopifnot(length(ld_list) == R, length(annot_list) == R, R >= 1L)
  if (R < 4L) {
    stop("PolyFun needs at least four regions (two per odd/even set) to choose ",
         "the ridge penalty out of chromosome; got ", R, call. = FALSE)
  }
  setup_polyfun_ldsc()
  if (is.null(n_ld)) n_ld <- rep(n, R)
  n_ld <- rep_len(as.numeric(n_ld), R)
  n_annot <- unique(vapply(annot_list, ncol, integer(1)))
  if (length(n_annot) != 1L) {
    stop("every region must carry the same annotation columns", call. = FALSE)
  }

  # --- per-region LD scores (the base annotation is column 1) --------------
  A_list <- vector("list", R); L_list <- vector("list", R)
  wld_list <- vector("list", R); R2_list <- vector("list", R)
  for (i in seq_len(R)) {
    z <- as.numeric(z_list[[i]]); LD <- as.matrix(ld_list[[i]])
    A <- cbind(base = 1, as.matrix(annot_list[[i]]))
    if (nrow(A) != length(z) || any(dim(LD) != length(z))) {
      stop("region ", i, ": z, LD and annotations have inconsistent sizes",
           call. = FALSE)
    }
    R2 <- .pf_r2_unbiased(LD, n_ld[i])
    R2_list[[i]]  <- R2
    A_list[[i]]   <- A
    L_list[[i]]   <- R2 %*% A
    wld_list[[i]] <- rowSums(R2)
  }
  chr    <- rep(seq_len(R), vapply(z_list, length, integer(1)))
  chisq  <- unlist(lapply(z_list, function(z) as.numeric(z)^2), use.names = FALSE)
  L_all  <- do.call(rbind, L_list)
  A_all  <- do.call(rbind, A_list)
  wld    <- unlist(wld_list, use.names = FALSE)
  M      <- colSums(A_all)

  # --- steps 1-2: L2-regularised S-LDSC, opposite-parity per-SNP h2 --------
  s1 <- .pf_ridge_step(L_all, wld, chisq, chr, n, M, num_lambdas)
  snpvar_ridge <- numeric(length(chisq))
  for (r in seq_len(R)) {
    other <- if (r %% 2L == 0L) "odd" else "even"
    rows  <- which(chr == r)
    snpvar_ridge[rows] <- as.numeric(A_all[rows, , drop = FALSE] %*%
                                       s1$tau[[other]][seq_len(ncol(A_all))])
  }

  # --- step 3: 20 bins by Ckmedian.1d.dp ------------------------------------
  bin <- .pf_bins(snpvar_ridge, num_bins, min_bin_size)
  n_bins <- max(bin)
  B_all <- matrix(0, length(bin), n_bins)
  B_all[cbind(seq_along(bin), bin)] <- 1

  # --- step 4: non-negative S-LDSC on the bins, same parity minus target ----
  Lb_all <- do.call(rbind, lapply(seq_len(R), function(i)
    R2_list[[i]] %*% B_all[chr == i, , drop = FALSE]))
  snpvar <- .pf_bins_step(Lb_all, B_all, wld, chisq, chr, n, colSums(B_all))

  # constrain_range: floor at max / q, rescale to the same total
  h2_total <- sum(snpvar)
  if (!is.finite(h2_total) || h2_total <= 0) {
    stop("PolyFun estimated zero heritability in every bin", call. = FALSE)
  }
  snpvar[snpvar < max(snpvar) / q] <- max(snpvar) / q
  snpvar <- snpvar * h2_total / sum(snpvar)

  # --- step 5: prior proportional to per-SNP h2 within each locus -----------
  split_r <- function(v) lapply(seq_len(R), function(i) v[chr == i])
  sv <- split_r(snpvar)
  list(
    prior        = lapply(sv, function(v) v / sum(v)),
    snpvar       = sv,
    snpvar_ridge = split_r(snpvar_ridge),
    bins         = split_r(bin),
    taus_ridge   = s1$tau,
    lambda       = s1$lambda
  )
}


# =============================================================================
# run_polyfun_ldsc(): PolyFun + SuSiE on one region
# =============================================================================

#' Run PolyFun + SuSiE on a single region
#'
#' Fine-maps one region with SuSiE as PolyFun does (finemapper.py), using
#' the prior causal probabilities from \code{\link{polyfun_priors}}. With
#' \code{prior = NULL} this is PolyFun's non-functionally-informed mode.
#'
#' @param z Numeric vector. Marginal z-scores.
#' @param LD Matrix. LD correlation matrix.
#' @param n Integer. GWAS sample size.
#' @param prior Numeric vector or NULL. Prior causal probabilities for the
#'   region (from \code{polyfun_priors()}); NULL gives a uniform prior.
#' @param L Integer. Maximum number of causal SNPs. Default 10 (paper).
#' @param hess Logical. Set the causal effect variance by modified HESS.
#'   Default TRUE (paper).
#' @param hess_iter Integer. HESS draws to average. Default 100.
#' @param max_iter Integer. SuSiE IBSS iterations. Default 100.
#' @param variant_ids Character or NULL. Optional variant labels.
#' @param ... Ignored (for wrapper compatibility).
#'
#' @return List with pip, credible_sets, method = "polyfun_ldsc", params,
#'   runtime_seconds and additional (prior_weights, prior_var, h2_hess,
#'   prior_source).
#' @export
run_polyfun_ldsc <- function(z, LD, n, prior = NULL, L = 10, hess = TRUE,
                             hess_iter = 100L, max_iter = 100L,
                             variant_ids = NULL, ...) {
  setup_polyfun_ldsc("susieR")
  t0 <- Sys.time()
  p <- length(z)
  stopifnot("LD must be p x p" = is.matrix(LD) && all(dim(LD) == c(p, p)))
  if (!is.null(prior)) {
    stopifnot(length(prior) == p, all(is.finite(prior)), all(prior > 0))
    prior <- prior / sum(prior)
  }

  prior_var <- NULL; h2_hess <- NA_real_
  if (hess) {
    if (L == 1) stop("HESS cannot be used with a single causal SNP per locus",
                     call. = FALSE)
    h2_hess <- .pf_hess_h2(z, LD, n, n_samples = hess_iter)
    prior_var <- h2_hess / L
    if (prior_var <= 0) {
      stop("HESS estimates that the locus causally explains zero heritability",
           call. = FALSE)
    }
    if (prior_var >= 1) {
      stop("HESS-estimated prior-var >1. The HESS estimator cannot be used in ",
           "this locus.", call. = FALSE)
    }
  }

  # susieR prints a reminder that residual-variance estimation assumes
  # in-sample LD; PolyFun estimates it regardless, so the note is silenced.
  fit <- suppressMessages(susieR::susie_rss(
    bhat = as.numeric(z), shat = rep(1, p), R = LD, n = n, L = L,
    scaled_prior_variance      = if (is.null(prior_var)) 1e-4 else prior_var,
    estimate_prior_variance    = is.null(prior_var),
    residual_variance          = NULL,
    estimate_residual_variance = TRUE,
    max_iter      = max_iter,
    prior_weights = prior
  ))

  cs_list <- if (is.null(fit$sets$cs)) list() else lapply(fit$sets$cs, as.integer)
  list(
    pip             = as.numeric(susieR::susie_get_pip(fit)),
    credible_sets   = cs_list,
    method          = "polyfun_ldsc",
    input_type      = "summary",
    params          = list(L = L, hess = hess, hess_iter = hess_iter,
                           max_iter = max_iter),
    runtime_seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    additional      = list(
      prior_weights = if (is.null(prior)) rep(1 / p, p) else prior,
      prior_var     = prior_var,
      h2_hess       = h2_hess,
      prior_source  = if (is.null(prior)) "non_functional" else "polyfun",
      variant_ids   = variant_ids
    )
  )
}


# =============================================================================
# run_polyfun_ldsc_region(): the adapter called by run_methods()
# =============================================================================

#' Region adapter for PolyFun + SuSiE
#'
#' Looks up the region's prior in the scenario-level result of
#' \code{run_polyfun_ldsc_scenario_setup()}. A region with annotations but
#' no prior fails, because PolyFun's prior needs every region of the
#' scenario.
#'
#' @param region_geno List. Genotype-side region data.
#' @param region_pheno List. Phenotype-side region data (contains \code{z}).
#' @param .polyfun_cache Named list of priors keyed by z fingerprint, from
#'   the scenario setup. Default NULL.
#' @param .polyfun_error Character or NULL. Error from the scenario setup.
#' @param ... Passed to \code{\link{run_polyfun_ldsc}}.
#'
#' @return The output of \code{\link{run_polyfun_ldsc}}.
#' @export
run_polyfun_ldsc_region <- function(region_geno, region_pheno,
                                    .polyfun_cache = NULL,
                                    .polyfun_error = NULL, ...) {
  A <- region_geno$annotations_matrix
  if (is.null(A)) A <- region_pheno$annotations_matrix
  z <- region_pheno$z
  prior <- NULL
  if (!is.null(A)) {
    prior <- if (is.null(.polyfun_cache)) NULL else
      .polyfun_cache[[.fb_fingerprint(z)]]
    if (is.null(prior)) {
      stop(if (!is.null(.polyfun_error)) .polyfun_error else
             "PolyFun priors were not computed for this region; run polyfun_ldsc through run_methods() with every region of the scenario",
           call. = FALSE)
    }
  }
  run_polyfun_ldsc(z = z, LD = region_geno$LD, n = region_geno$n,
                   prior = prior, variant_ids = region_geno$variant_ids, ...)
}


# =============================================================================
# Scenario-level setup hook: steps 1-5 once per scenario
# =============================================================================

#' Scenario-level PolyFun priors
#'
#' Runs \code{\link{polyfun_priors}} over all regions of a scenario and
#' returns the priors keyed by z fingerprint. Returns an empty list when the
#' scenario has no annotations (PolyFun then runs non-functionally), and an
#' error message when the priors cannot be estimated.
#'
#' @param genotypes List of per-region genotype data.
#' @param regions List of per-region phenotype data for one scenario.
#' @param user_args Named list of method arguments (\code{num_bins},
#'   \code{q} are honoured).
#'
#' @return A named list merged into each region's arguments.
#' @export
run_polyfun_ldsc_scenario_setup <- function(genotypes, regions, user_args) {
  R <- length(regions)
  A_list <- lapply(seq_len(R), function(i) {
    A <- genotypes[[i]]$annotations_matrix
    if (is.null(A)) regions[[i]]$annotations_matrix else A
  })
  if (all(vapply(A_list, is.null, logical(1)))) return(list())
  res <- tryCatch({
    if (any(vapply(A_list, is.null, logical(1)))) {
      stop("some regions have no annotation matrix", call. = FALSE)
    }
    n_ld <- vapply(seq_len(R), function(i) {
      nr <- genotypes[[i]]$n_ref
      as.numeric(if (is.null(nr)) genotypes[[i]]$n else nr)
    }, numeric(1))
    pf <- polyfun_priors(
      z_list     = lapply(regions, `[[`, "z"),
      ld_list    = lapply(genotypes, `[[`, "LD"),
      annot_list = A_list,
      n          = genotypes[[1]]$n,
      n_ld       = n_ld,
      num_bins   = if (is.null(user_args$num_bins)) 20L else user_args$num_bins,
      q          = if (is.null(user_args$q)) 100 else user_args$q
    )
    cache <- stats::setNames(pf$prior,
                             vapply(regions, function(r) .fb_fingerprint(r$z),
                                    character(1)))
    list(.polyfun_cache = cache)
  }, error = function(e) list(.polyfun_error = paste("PolyFun prior estimation failed:",
                                                     conditionMessage(e))))
  res
}


# =============================================================================
# Internal helpers
# =============================================================================

# Unbiased r^2 (compute_ldscores_from_ld.compute_R2_unbiased).
.pf_r2_unbiased <- function(R, n) {
  R * R * ((n - 1) / (n - 2)) - 1 / (n - 2)
}

# chi2 filter of PolyFun.run_ldsc (MAX_CHI2 = 80).
.pf_chisq_keep <- function(chisq, n) {
  chisq < max(0.001 * n, 80)
}

# LDSC regression weights on the old_weights path: Hsq.weights at the
# aggregate h2 estimate with the null intercept, then sqrt, then scaled to
# sum to one (IRWLS._weight).
.pf_ldsc_weights <- function(ld_tot, w_ld, chisq, n, M_tot) {
  hsq <- M_tot * (mean(chisq) - 1) / mean(ld_tot * n)
  hsq <- min(max(hsq, 0), 1)
  ld  <- pmax(ld_tot, 1)
  wl  <- pmax(w_ld, 1)
  cc  <- hsq * n / M_tot
  w   <- 1 / (2 * (1 + cc * ld)^2) * (1 / wl)
  w   <- sqrt(w)
  if (any(w <= 0)) stop("Weights must be > 0", call. = FALSE)
  w / sum(w)
}

# sklearn.metrics.r2_score.
.pf_r2_score <- function(y, yhat) {
  1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
}

# Ridge solution with an unpenalised intercept in the last column.
.pf_ridge_solve <- function(XtX, Xty, lambda) {
  pen <- diag(lambda, nrow(XtX))
  pen[nrow(XtX), nrow(XtX)] <- 0
  as.numeric(solve(XtX + pen, Xty))
}

# Leave-one-chromosome-out choice of lambda (Jackknife_Ridge._find_best_lambda).
.pf_best_lambda <- function(x, y, chr, lambdas) {
  chrs <- sort(unique(chr))
  if (length(chrs) < 2L) stop("need at least two chromosomes to choose lambda",
                              call. = FALSE)
  XtX <- crossprod(x); Xty <- crossprod(x, y)
  pred <- matrix(NA_real_, length(y), length(lambdas))
  for (cc in chrs) {
    idx <- which(chr == cc)
    xc <- x[idx, , drop = FALSE]
    XtX_l <- XtX - crossprod(xc)
    Xty_l <- Xty - crossprod(xc, y[idx])
    taus <- vapply(lambdas, function(l) .pf_ridge_solve(XtX_l, Xty_l, l),
                   numeric(ncol(x)))
    pred[idx, ] <- xc %*% taus
  }
  scores <- apply(pred, 2L, function(p) .pf_r2_score(y, p))
  best <- which.max(scores)
  list(lambda = lambdas[best], score = scores[best])
}

# Steps 1: L2-regularised S-LDSC fitted separately on the even and odd sets.
# Returns per-set coefficient vectors (annotation columns only, divided by
# the column norms and mean N) and the chosen lambdas.
.pf_ridge_step <- function(L_all, wld, chisq, chr, n, M, num_lambdas) {
  keep <- .pf_chisq_keep(chisq, n)
  L <- L_all[keep, , drop = FALSE]; y <- chisq[keep]; ck <- chr[keep]
  w <- .pf_ldsc_weights(rowSums(L), wld[keep], y, n, sum(M))
  x <- cbind(L, 1)
  xm <- colMeans(x); xm[length(xm)] <- 0
  xc <- sweep(x, 2L, xm) * w
  yw <- y * w
  x_l2 <- sqrt(colSums(xc^2))
  if (any(x_l2 <= 0)) {
    stop("an annotation has no variation across the regression SNPs",
         call. = FALSE)
  }
  xs <- sweep(xc, 2L, x_l2, "/")
  mean_diag <- mean(colSums(xs^2))
  lambdas <- 10^seq(log10(mean_diag * 1e-8), log10(mean_diag * 1e2),
                    length.out = num_lambdas)
  sets <- list(even = sort(unique(ck[ck %% 2L == 0L])),
               odd  = sort(unique(ck[ck %% 2L == 1L])))
  tau <- list(); lambda <- c()
  for (s in names(sets)) {
    in_set <- ck %in% sets[[s]]
    bl <- .pf_best_lambda(xs[in_set, , drop = FALSE], yw[in_set], ck[in_set],
                          lambdas)
    est <- .pf_ridge_solve(crossprod(xs[in_set, , drop = FALSE]),
                           crossprod(xs[in_set, , drop = FALSE], yw[in_set]),
                           bl$lambda)
    est <- est / x_l2
    tau[[s]] <- est[seq_len(ncol(L))] / n
    lambda[s] <- bl$lambda
  }
  list(tau = tau, lambda = lambda)
}

# Step 3: Ckmedian.1d.dp into num_bins bins, then PolyFun's rebalancing.
.pf_bins <- function(snpvar, num_bins, min_bin_size) {
  ord <- order(snpvar)
  seg <- Ckmeans.1d.dp::Ckmedian.1d.dp(snpvar[ord], k = num_bins)
  sizes <- .pf_rebalance_bins(as.integer(seg$size), min_bin_size)
  bin <- integer(length(snpvar))
  bin[ord] <- rep(seq_along(sizes), sizes)
  bin
}

# PolyFun.create_df_bins: move SNPs from the previous bin into any bin
# (other than the first) holding fewer than min_bin_size SNPs.
.pf_rebalance_bins <- function(sizes, min_bin_size) {
  i <- length(sizes)
  if (i < 2L) return(sizes)
  repeat {
    if (sizes[i] >= min_bin_size) {
      i <- i - 1L
      if (i == 1L) break
      next
    }
    move <- min(min_bin_size - sizes[i], sizes[i - 1L])
    sizes[i] <- sizes[i] + move
    sizes[i - 1L] <- sizes[i - 1L] - move
    if (sizes[i - 1L] == 0L) {
      sizes <- sizes[-(i - 1L)]
      i <- i - 1L
    }
    if (sizes[i] >= min_bin_size) i <- i - 1L
    if (i == 1L) break
  }
  sizes
}

# Step 4: non-negative S-LDSC on the bin LD scores; for each target
# chromosome, fit on the other chromosomes of the same parity.
.pf_bins_step <- function(Lb_all, B_all, wld, chisq, chr, n, M_bins) {
  keep <- .pf_chisq_keep(chisq, n)
  Lb <- Lb_all[keep, , drop = FALSE]; y <- chisq[keep]; ck <- chr[keep]
  w <- .pf_ldsc_weights(rowSums(Lb), wld[keep], y, n, sum(M_bins))
  x <- cbind(Lb, 1) * w
  yw <- y * w
  snpvar <- numeric(nrow(B_all))
  for (r in sort(unique(chr))) {
    loco <- (ck %% 2L == r %% 2L) & (ck != r)
    if (!any(loco)) {
      stop("chromosome ", r, " has no other chromosome of the same parity",
           call. = FALSE)
    }
    est <- .nnls_lawson_hanson(x[loco, , drop = FALSE], yw[loco])
    taus <- est[seq_len(ncol(Lb))] / n
    rows <- which(chr == r)
    snpvar[rows] <- as.numeric(B_all[rows, , drop = FALSE] %*% taus)
  }
  snpvar
}

# Non-negative least squares, Lawson & Hanson (1974), as in
# scipy.optimize.nnls.
.nnls_lawson_hanson <- function(A, b, max_iter = 3L * ncol(A)) {
  m <- ncol(A)
  x <- numeric(m)
  passive <- logical(m)
  tol <- 10 * .Machine$double.eps * norm(A, "1") * max(dim(A))
  AtA <- crossprod(A); Atb <- as.numeric(crossprod(A, b))
  w <- Atb - as.numeric(AtA %*% x)
  iter <- 0L
  while (any(!passive) && max(w[!passive]) > tol) {
    cand <- which(!passive)
    passive[cand[which.max(w[cand])]] <- TRUE
    repeat {
      iter <- iter + 1L
      if (iter > max_iter) {
        warning("NNLS reached the iteration limit", call. = FALSE)
        return(x)
      }
      P <- which(passive)
      s <- numeric(m)
      s[P] <- solve(AtA[P, P, drop = FALSE], Atb[P])
      if (all(s[P] > tol)) break
      neg <- P[s[P] <= tol]
      alpha <- min(x[neg] / (x[neg] - s[neg]))
      x <- x + alpha * (s - x)
      passive[passive & x <= tol] <- FALSE
      x[!passive] <- 0
    }
    x <- s
    w <- Atb - as.numeric(AtA %*% x)
  }
  x
}

# Modified HESS (finemapper.estimate_h2_hess_wrapper, no min_h2 bound).
.pf_hess_h2 <- function(z, LD, n, prop_keep = 0.005, R_cutoff = 0.99,
                        n_samples = 100L) {
  pvals <- 2 * stats::pnorm(-abs(z))
  cutoff <- stats::quantile(pvals, prop_keep, type = 7, names = FALSE)
  if (cutoff == 0) cutoff <- min(pvals[pvals > 0])
  pot <- which(pvals < cutoff)
  if (length(pot) == 0L) return(0)
  Rp <- LD[pot, pot, drop = FALSE]
  adj <- abs(Rp) > R_cutoff
  diag(adj) <- FALSE
  a_all <- z[pot] / sqrt(n)
  h2 <- vapply(seq_len(n_samples), function(s) {
    inds <- sort(.pf_random_mis(adj))
    a <- a_all[inds]
    sum(a * solve(Rp[inds, inds, drop = FALSE], a)) - length(inds) / n
  }, numeric(1))
  mean(h2)
}

# A random maximal independent set (networkx.maximal_independent_set): start
# from a random node, then repeatedly add a random node that is not adjacent
# to the set.
.pf_random_mis <- function(adj) {
  m <- nrow(adj)
  start <- sample.int(m, 1L)
  indep <- start
  available <- setdiff(seq_len(m), c(start, which(adj[start, ])))
  while (length(available) > 0L) {
    node <- available[sample.int(length(available), 1L)]
    indep <- c(indep, node)
    available <- setdiff(available, c(node, which(adj[node, ])))
  }
  indep
}
