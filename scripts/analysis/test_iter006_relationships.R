#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/test_iter006_relationships.R
#
# >>> ITERATION 006 ONLY - TEMPORARY. See docs/autoresearch/iteration-006-REVERT.md
#
# Checks the Iteration 006 annotation -> causality relationships in
# R/simulate_phenotypes.R against brute-force re-implementations written
# independently here, and checks that nothing Iteration 005 relies on moved.
#
#   Rscript scripts/analysis/test_iter006_relationships.R
# =============================================================================
`%||%` <- function(a, b) if (is.null(a)) b else a
env <- new.env(); sys.source("R/simulate_phenotypes.R", envir = env)
n_pass <- 0L; n_fail <- 0L
check <- function(name, ok) {
  if (isTRUE(ok)) { n_pass <<- n_pass + 1L; cat("  PASS ", name, "\n") }
  else            { n_fail <<- n_fail + 1L; cat("  FAIL ", name, "\n") }
}
near <- function(a, b, tol = 1e-12) isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE))

# Brute force: loop over positions, renormalise over the neighbours that exist.
bf_avg <- function(x, w) {
  h <- (length(w) - 1) / 2; p <- length(x)
  vapply(seq_len(p), function(j) { i <- j + (-h:h); ok <- i >= 1 & i <= p
    sum(w[ok] * x[i[ok]]) / sum(w[ok]) }, 0)
}
w2 <- c(1, 2, 3, 2, 1) / 3
wc <- sapply(-10:10, function(d) exp(-abs(d) / 5) * (1 + cos(2 * pi * d / 10)) / 2)

cat("weights\n")
check("wavg2 weights are 1/3, 2/3, 1, 2/3, 1/3", near(env$.ITER006_W_WAVG2, w2))
check("cosine10 weights: 1 at d=0, 0 at |d|=5, exp(-2) at |d|=10",
      near(env$.ITER006_W_COS10[c(11, 6, 16, 1, 21)], c(1, 0, 0, exp(-2), exp(-2)), 1e-12))

cat("neighbour average against brute force (interior and both ends)\n")
set.seed(1)
for (p in c(7, 25, 1000)) {
  x <- rnorm(p)
  check(sprintf("wavg2, p = %d", p),    near(env$.iter006_neighbour_avg(x, w2), bf_avg(x, w2)))
  check(sprintf("cosine10, p = %d", p), near(env$.iter006_neighbour_avg(x, wc), bf_avg(x, wc)))
}

cat("scores against brute force\n")
A <- matrix(rnorm(1000 * 10), 1000, 10); I <- A[, 1:5]
thr <- function(a) a * (a > 1)
BF <- list(
  wavg2           = rowSums(apply(I, 2, bf_avg, w = w2)),
  valthresh       = rowSums(ifelse(I > 1, I, 0)),
  valthresh_wavg2 = rowSums(apply(thr(I), 2, bf_avg, w = w2)),
  square          = rowSums(I * I),
  cubic           = rowSums(I * I * I),
  cosine10        = rowSums(apply(I, 2, bf_avg, w = wc)))
for (rel in names(BF))
  check(sprintf("%s = brute force, via .causal_log_weights (lambda = log 5.4)", rel),
        near(env$.causal_log_weights(A, rel, log(5.4), 5L), log(5.4) * BF[[rel]]))
check("only the first n_informative columns enter",
      near(env$.causal_log_weights(cbind(A[, 1:5], 100 * A[, 6:10]), "square", 1, 5L),
           env$.causal_log_weights(A, "square", 1, 5L)))
check("an unknown relationship still stops",
      inherits(try(env$.causal_log_weights(A, "banana", 1, 5L), silent = TRUE), "try-error"))

cat("Iteration 005 relationships unchanged\n")
check("additive = lambda * sum", near(env$.causal_log_weights(A, "additive", 2, 5L), 2 * rowSums(I)))
check("null = zeros",            near(env$.causal_log_weights(A, "null", 2, 5L), rep(0, 1000)))
check("nonmono formula",         near(env$.causal_log_weights(A, "nonmono", 2, 5L), 2 * exp(-((rowSums(I) - 2)^2) / 2)))

# Iteration 006 runs at fold 2.7 (top-decile target 0.349, the user's choice);
# fold 5.4 (0.605, Iteration 005's level) is checked too.
top10 <- function(pr) sum(sort(pr, decreasing = TRUE)[1:100])
for (fold in c(2.7, 5.4)) {
  target <- env$.concentration_target(fold)
  cat(sprintf("strength matching through select_causal_variants, fold %.1f -> top-decile target %.3f\n", fold, target))
  for (rel in c("additive", env$.ITER006_RELATIONSHIPS)) {
    r <- env$select_causal_variants(p = 1000, S = 3, annotation_matrix = A,
                                    enrichment = c(rep(fold, 5), rep(1, 5)), n_annotations = 10,
                                    annotation_type = "continuous", relationship = rel, n_informative = 5L)
    check(sprintf("%-16s top-decile share %.3f within 0.01 of %.3f; 3 causal drawn", rel, top10(r$causal_probs), target),
          abs(top10(r$causal_probs) - target) < 0.01 && length(r$causal_indices) == 3L)
  }
}
r <- env$select_causal_variants(p = 1000, S = 3, annotation_matrix = A, enrichment = rep(1, 10),
                                n_annotations = 10, annotation_type = "continuous",
                                relationship = "null", n_informative = 5L)
check("null arm is uniform", near(r$causal_probs, rep(1 / 1000, 1000)))

cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
quit(status = if (n_fail) 1L else 0L)
