# =============================================================================
# scripts/analysis/iter004_lib.R
#
# Shared machinery for the Iteration 004 variable-importance analysis.
# Implements variable-importance-analysis.md §2.3 (strata), §6 (prep and
# assertions), §3.4b (variance components) and §7.2-7.3 (the noise-corrected
# decomposition).
#
# The hard part is §7.2. Everything else is bookkeeping.
#
# STATUS: written, never executed. No part of this has been run against data.
# =============================================================================

suppressWarnings(suppressMessages({
  if (requireNamespace("fmbenchmark", quietly = TRUE)) library(fmbenchmark)
}))

`%||%` <- function(x, y) if (is.null(x)) y else x

# The six non-poolable strata (§2.3). p_causal and enrichment_fold are NESTED,
# not crossed, so a single ANOVA over all the data is not well defined - and the
# nesting happens to align exactly with the never-pool rule, so one split solves
# both problems. `model` and `annotation_type` can therefore NEVER appear in an
# importance ranking: you cannot measure the importance of a variable you
# stratify on. They are the frame, not a result.
STRATA <- list(
  list(key = "sparse_none",       model = "sparse",     annot = "none"),
  list(key = "sparse_binary",     model = "sparse",     annot = "binary"),
  list(key = "sparse_cont",       model = "sparse",     annot = "continuous"),
  list(key = "sparseinf_none",    model = "sparse_inf", annot = "none"),
  list(key = "sparseinf_binary",  model = "sparse_inf", annot = "binary"),
  list(key = "sparseinf_cont",    model = "sparse_inf", annot = "continuous")
)

# Factors that live at the JOB level, i.e. are constant across a job's 250
# scenarios and therefore carry region-draw noise (§3.4). region_size is not a
# job-level factor but occupies the same error stratum, because two levels of
# region_size are always two different region draws.
JOB_LEVEL_FACTORS <- c("pc", "enrich")
WHOLE_PLOT_FACTORS <- c(JOB_LEVEL_FACTORS, "region_size")

# Factors swept WITHIN a job, on the same regions. Any term containing one of
# these annihilates the region draw exactly - its contrasts sum to zero within
# each block - so it sheds only iteration noise.
SUBPLOT_FACTORS <- c("S", "phi")

# m = the number of (S, phi) cells sharing one region draw. THIS FACTOR IS THE
# WHOLE POINT of §7.2: ubar is constant across the (S, phi) sweep, so its
# variance enters a whole-plot term's mean square multiplied by m. Correcting a
# whole-plot term with the iteration variance alone understates its noise content
# by two to three orders of magnitude - the classic split-plot error, and it
# inflates exactly the factors carrying the annotation conclusions.
M_SHARED <- 25L      # 5 S x 5 phi

N_ITER    <- 10L     # replicates per cell
N_REGIONS <- 2L      # same-size regions per replicate
N_FITS    <- N_ITER * N_REGIONS   # = 20 fits per cell


# ---------------------------------------------------------------------------
# §6 - prepare the analysis table
# ---------------------------------------------------------------------------

#' Filter to in-sample LD, then recode NA-as-level. IN THAT ORDER.
#'
#' aggregate(), split(), table() and sort(unique(.)) all drop NA silently, so
#' grouping naively deletes whole arms without warning, and
#' `subset(d, n_ref == NA)` returns zero rows rather than erroring. The filter
#' must therefore act on the RAW NA before anything is recoded.
#'
#' For Iteration 004 the n_ref filter is a NO-OP - the grid is in-sample by
#' construction - but the assertion still validates the row count, which is the
#' point of keeping it.
#'
#' @param d Data frame from the collect step, one row per (cell, method) or
#'   per (cell, method, replicate).
#' @param expect_jobs Expected number of distinct job_dir values after filtering.
prepare_analysis_table <- function(d, expect_jobs = 45L) {
  # 1. FILTER on the raw NA.
  if ("n_ref" %in% names(d)) d <- d[is.na(d$n_ref), , drop = FALSE]
  n_jobs <- length(unique(d$job_dir))
  if (n_jobs != expect_jobs) {
    stop(sprintf(
      "expected %d job_dirs after the in-sample filter, found %d. Either the grid changed or the filter dropped rows it should not have.",
      expect_jobs, n_jobs), call. = FALSE)
  }

  # 2. RECODE - NA is a MEANINGFUL LEVEL for these two columns, not missingness.
  #    enrichment_fold = NA is the `none` arm; p_causal = NA is the sparse model.
  d$enrich <- factor(ifelse(is.na(d$enrichment_fold), "none",
                            as.character(d$enrichment_fold)),
                     levels = c("none", "2.7", "5.4", "8.1", "10.8"))
  d$pc     <- factor(ifelse(is.na(d$p_causal), "sparse",
                            as.character(d$p_causal)),
                     levels = c("sparse", "0.2", "0.4", "0.6", "0.8"))
  d$S           <- factor(d$S)
  d$phi         <- factor(d$phi)
  d$region_size <- factor(d$region_size)

  # 3. ASSERT - a silent NA drop is the failure mode this guards against.
  key <- c("enrich", "pc", "S", "phi", "region_size")
  if (anyNA(d[key])) {
    stop("NA present in a grouping factor after recoding: ",
         paste(key[vapply(d[key], anyNA, logical(1))], collapse = ", "),
         call. = FALSE)
  }

  # Ordered factors get polynomial contrasts (§7.4): partitioning 4 df into
  # linear/quadratic/cubic/quartic makes "AP falls linearly in log phi" a
  # 1-df claim rather than a 4-df "the levels differ".
  for (v in c("S", "phi", "region_size", "enrich")) {
    if (nlevels(d[[v]]) > 2L) contrasts(d[[v]]) <- contr.poly(nlevels(d[[v]]))
  }
  d
}


#' Subset to one stratum and drop the factors that are constant within it.
#'
#' Returns NULL if the stratum is absent (e.g. asking for sparse_inf when only
#' the sparse arm has been collected).
stratum_subset <- function(d, stratum) {
  s <- d[d$model == stratum$model & d$annotation_type == stratum$annot, ,
         drop = FALSE]
  if (!nrow(s)) return(NULL)
  free <- character(0)
  for (v in c("pc", "enrich", "S", "phi", "region_size")) {
    if (nlevels(droplevels(s[[v]])) > 1L) {
      s[[v]] <- droplevels(s[[v]])
      if (nlevels(s[[v]]) > 2L) contrasts(s[[v]]) <- contr.poly(nlevels(s[[v]]))
      free <- c(free, v)
    }
  }
  list(data = s, free = free, key = stratum$key)
}


# ---------------------------------------------------------------------------
# §3.4b - variance components
# ---------------------------------------------------------------------------

#' Estimate sigma^2_eps and sigma^2_u from the FIT-level table.
#'
#' The design is split-plot with three error strata. A "job" is a row of
#' params_grid.csv and does NOT include S or phi; those are swept inside a job on
#' THE SAME REGIONS. Genotypes and annotation matrices are simulated once per job
#' and shared by all 250 of its scenarios, so an iteration resamples the causal
#' draw, the effect sizes and the phenotype noise, and nothing else.
#'
#' sigma^2_u is identified only because n_ref is excluded: the sole between-job
#' nuisance is then which loci and annotation matrices were drawn, and that IS
#' replicated inside a job (2 regions per size class x 5 classes x 45 jobs = 225
#' within-job region pairs). With n_ref in the design this would not work, since a
#' reference-panel realisation is drawn once per job and never replicated within.
#'
#'   sigma^2_eps       = within-region variance across the 10 iterations
#'   sigma^2_u(g)      = (1/2) Var(Ybar_r1 - Ybar_r2) - sigma^2_eps/10
#'   sigma^2_job       ~ sigma^2_u / 2      (a cell mean averages the 2 regions)
#'
#' Estimated PER SIZE CLASS, not pooled over all 225 pairs: sigma^2_u for a
#' p=500 draw and a p=2000 draw are unlikely to be equal.
#'
#' TWO ASSUMPTIONS TO STATE WHEN REPORTING (spec A1, A2):
#'   A1 exchangeability - regions are independent draws within a job, so a job's
#'      panel carries no component beyond its constituent region draws.
#'   A2 no region x scenario interaction - u_r is taken constant across S and phi.
#'      If it is not, the pair-difference estimator absorbs Var(w) into
#'      sigma^2_u: conservative for whole-plot terms, but the subplot correction
#'      becomes an UNDER-correction. Diagnose by computing sigma^2_u separately
#'      per (S, phi) and looking for heterogeneity - see `diagnose_A2()`.
#'
#' @param fits Fit-level table with columns job_dir, S, phi, iter, region_size,
#'   region_idx and the response.
#' @param response Name of the response column.
#' @return list(sigma2_eps, sigma2_u (named by size class), sigma2_u_bar,
#'   sigma2_job, n_pairs)
variance_components <- function(fits, response) {
  y <- fits[[response]]
  if (is.null(y)) stop("response '", response, "' not in the fit table", call. = FALSE)
  ok <- is.finite(y)
  f  <- fits[ok, , drop = FALSE]
  f$.y <- y[ok]

  # sigma^2_eps: variance across the 10 iterations WITHIN one region, pooled.
  # This is pure replicate error - the region, the loci and the annotations are
  # all literally identical across those 10 values.
  reg_key <- interaction(f$job_dir, f$region_size, f$region_idx,
                         f$S, f$phi, drop = TRUE)
  within  <- tapply(f$.y, reg_key, function(v) if (length(v) > 1L) var(v) else NA_real_)
  sigma2_eps <- mean(within, na.rm = TRUE)

  # Region means over the 10 iterations, then the within-(job, size) pair contrast.
  rm_key <- interaction(f$job_dir, f$region_size, f$region_idx, drop = TRUE)
  reg_mean <- tapply(f$.y, rm_key, mean, na.rm = TRUE)
  meta <- do.call(rbind, strsplit(names(reg_mean), ".", fixed = TRUE))
  rmdf <- data.frame(job = meta[, 1], size = meta[, 2], idx = meta[, 3],
                     ybar = as.numeric(reg_mean), stringsAsFactors = FALSE)

  sizes <- sort(unique(rmdf$size))
  s2u <- setNames(rep(NA_real_, length(sizes)), sizes)
  npair <- setNames(rep(0L, length(sizes)), sizes)
  for (g in sizes) {
    sub <- rmdf[rmdf$size == g, , drop = FALSE]
    sp  <- split(sub$ybar, sub$job)
    dif <- vapply(sp, function(v) if (length(v) == 2L) v[1] - v[2] else NA_real_,
                  numeric(1))
    dif <- dif[is.finite(dif)]
    npair[g] <- length(dif)
    if (length(dif) > 1L) {
      # Var(r1 - r2) = 2 sigma^2_u + 2 sigma^2_eps/10  =>  halve, then subtract.
      s2u[g] <- max(0, 0.5 * var(dif) - sigma2_eps / N_ITER)
    }
  }

  list(sigma2_eps = sigma2_eps,
       sigma2_u   = s2u,
       # Use the MEAN over size classes in lambda_T. Exact for whole-plot terms
       # that average over region_size, first-order for terms containing it; the
       # per-class values show how much that approximation matters.
       sigma2_u_bar = mean(s2u, na.rm = TRUE),
       sigma2_job   = mean(s2u, na.rm = TRUE) / 2,
       n_pairs      = npair)
}


#' Diagnostic for assumption A2 (§3.4b).
#'
#' Strong heterogeneity of sigma^2_u across (S, phi) means u_r is not constant
#' across scenarios, so the subplot correction is an under-correction.
diagnose_A2 <- function(fits, response) {
  sp <- split(fits, interaction(fits$S, fits$phi, drop = TRUE))
  out <- lapply(names(sp), function(k) {
    vc <- tryCatch(variance_components(sp[[k]], response), error = function(e) NULL)
    if (is.null(vc)) return(NULL)
    data.frame(cell = k, sigma2_u = vc$sigma2_u_bar, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(out)) return(NULL)
  out$ratio_to_median <- out$sigma2_u / median(out$sigma2_u, na.rm = TRUE)
  out[order(-out$ratio_to_median), ]
}


# ---------------------------------------------------------------------------
# §7.2-7.3 - the decomposition
# ---------------------------------------------------------------------------

#' Is an ANOVA term on the whole-plot margin?
#'
#' A term is whole-plot iff it contains NO subplot factor. Any term containing S
#' or phi annihilates the shared region draw exactly, leaving only iteration
#' noise; every other term additionally carries m * sigma^2_u / 2 per df.
is_whole_plot <- function(term_label) {
  parts <- trimws(strsplit(term_label, ":", fixed = TRUE)[[1]])
  !any(parts %in% SUBPLOT_FACTORS)
}


#' Noise-corrected variance shares on the cell-mean table (§7.2, §7.3).
#'
#' Fits the full factorial on CELL MEANS - one row per cell, zero residual df, so
#' the orthogonal SS partition the total exactly. There is deliberately no
#' residual line: every df belongs to a term, and noise is handled by explicit
#' correction, never by a residual.
#'
#' Working on cell means (rather than one replicate-level ANOVA with a per-term
#' MS_E) keeps every quantity on ONE scale. At replicate level each noise source
#' enters a mean square multiplied by the number of rows sharing the draw
#' (hundreds), so subtracting a cell-mean-scale sigma^2_u/2 there under-corrects
#' by orders of magnitude while LOOKING compliant.
#'
#'   lambda_T   = sigma^2_eps/20 + (m sigma^2_u / 2) * 1[T on whole-plot margin]
#'   omega^2_T  = (SS_T - df_T lambda_T) / SS_tot
#'
#' @param cells Cell-mean table for ONE stratum and ONE method.
#' @param response Response column name.
#' @param free Character vector of free factor names in this stratum.
#' @param vc Output of variance_components().
#' @return list(table, noise_floor, ss_tot, ...)
decompose <- function(cells, response, free, vc) {
  stopifnot(length(free) >= 1L)
  fml <- as.formula(paste(response, "~", paste(free, collapse = " * ")))
  fit <- aov(fml, data = cells)
  a   <- summary(fit)[[1]]
  # aov() pads term labels; strip.
  terms_lab <- trimws(rownames(a))
  df_T <- a[["Df"]]
  ss_T <- a[["Sum Sq"]]

  # There should be no residual line on a saturated design. If one appears the
  # design is not balanced/complete and the closed forms below do not apply.
  if (any(terms_lab == "Residuals")) {
    keep <- terms_lab != "Residuals"
    warning("residual df present: the stratum is not a complete balanced ",
            "factorial (check n_fits / n_failed). The closed-form correction ",
            "assumes 20 fits in every cell; use the lmer route instead.",
            call. = FALSE)
    terms_lab <- terms_lab[keep]; df_T <- df_T[keep]; ss_T <- ss_T[keep]
  }

  wp <- vapply(terms_lab, is_whole_plot, logical(1))
  lambda <- vc$sigma2_eps / N_FITS +
            ifelse(wp, M_SHARED * vc$sigma2_u_bar / 2, 0)

  ss_tot <- sum(ss_T)
  noise  <- df_T * lambda
  omega2 <- (ss_T - noise) / ss_tot

  # Effect span (§7.3 ii): df-free, in metric units, immune to the level-count
  # artefact. Only defined for main effects.
  span <- vapply(terms_lab, function(t) {
    if (grepl(":", t, fixed = TRUE) || !(t %in% free)) return(NA_real_)
    mu <- tapply(cells[[response]], cells[[t]], mean, na.rm = TRUE)
    if (all(is.na(mu))) NA_real_ else max(mu, na.rm = TRUE) - min(mu, na.rm = TRUE)
  }, numeric(1))

  order_of <- vapply(terms_lab, function(t)
    length(strsplit(t, ":", fixed = TRUE)[[1]]), integer(1))

  tab <- data.frame(
    term        = terms_lab,
    df          = df_T,
    order       = order_of,
    whole_plot  = wp,
    ss          = ss_T,
    lambda      = lambda,
    noise_ss    = noise,
    omega2      = omega2,
    span        = span,
    stringsAsFactors = FALSE
  )
  tab <- tab[order(-tab$omega2), ]

  noise_floor <- sum(noise) / ss_tot
  list(
    table = tab,
    ss_tot = ss_tot,
    noise_floor = noise_floor,
    noise_floor_iter   = sum(df_T * (vc$sigma2_eps / N_FITS)) / ss_tot,
    noise_floor_region = sum(df_T[wp] * (M_SHARED * vc$sigma2_u_bar / 2)) / ss_tot,
    # Sobol normalisation, for comparing across strata or metrics whose noise
    # floors differ (§7.3 iii).
    pi_T = setNames(tab$omega2 / (1 - noise_floor), tab$term)
  )
}


#' First-order and total-order (Sobol) indices from the CORRECTED components.
#'
#' S_Ti - S_i is interaction involvement. ~0 means the factor acts additively and
#' its marginal plot tells the whole story; large means its effect is conditional
#' and the marginal is misleading. Computed from corrected components so a noisy
#' high-order term is not read as interaction structure.
sobol_indices <- function(dec, free) {
  V <- pmax(dec$table$omega2, 0)     # clamp only HERE, for the index ratio
  names(V) <- dec$table$term
  Vtot <- sum(V)
  out <- lapply(free, function(f) {
    first <- V[[f]] %||% 0
    tot <- sum(V[vapply(names(V), function(t)
      f %in% trimws(strsplit(t, ":", fixed = TRUE)[[1]]), logical(1))])
    data.frame(factor = f, S_i = first / Vtot, S_Ti = tot / Vtot,
               interaction_involvement = (tot - first) / Vtot,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}


#' The mandatory caveats (§7.5). Print these beside every importance table.
importance_caveats <- function() {
  c(
    "RANGE DEPENDENCE: omega^2 is defined relative to the input distribution, here uniform over the levels swept. Widening a factor's range mechanically inflates its importance. Say 'X explains 34% of variance OVER THE RANGE SWEPT', never 'X is the most important factor'.",
    "UNIFORM WEIGHTING: marginal means weight every combination of the other factors equally. If deployment concentrates in some corner, reweight - it can substantially reorder the ranking.",
    "METRIC SCALE: omega^2 is not invariant to monotone transforms. Bounded metrics saturate - a factor operating where AP ~ 0.98 looks unimportant because there is no room to move. For bounded rates decompose on the logit scale.",
    "IN-SAMPLE ONLY: Iteration 004 has no reference-panel LD. Iteration 003 showed in-sample and panel LD invert method rankings, so none of these conclusions transfer.",
    "NEGATIVE omega^2 IS NOT CLAMPED: null factors can go below zero. Many strongly negative shares mean lambda-hat is too large - a diagnostic, not an embarrassment."
  )
}
