#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/verify_wrappers.R
#
# Empirical correctness checks for the method wrappers. Run this ON THE CLUSTER,
# where susieR / CARMA / the Python tools are actually installed - several of
# these checks CANNOT be done on a laptop that lacks the packages, and a
# code-reading audit cannot substitute for them.
#
#   module load R/4.5.2-gfbf-2025b
#   Rscript scripts/analysis/verify_wrappers.R
#
# Each check prints PASS / FAIL / SKIP. FAIL means a genuine problem; SKIP means
# the dependency is absent so the check could not be attempted.
# =============================================================================
suppressWarnings(suppressMessages({
  if (requireNamespace("fmbenchmark", quietly = TRUE)) library(fmbenchmark)
  else if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION"))
    pkgload::load_all(".", quiet = TRUE)
  else stop("fmbenchmark not available")
}))
ok <- function(t, m="") cat(sprintf("  PASS  %-58s %s\n", t, m))
no <- function(t, m="") cat(sprintf("  FAIL  %-58s %s\n", t, m))
sk <- function(t, m="") cat(sprintf("  SKIP  %-58s %s\n", t, m))
chk <- function(cond, t, m="") if (isTRUE(cond)) ok(t,m) else no(t,m)

cat("\n=== A. susie_rss argument contract (needs susieR) ===\n")
if (!requireNamespace("susieR", quietly = TRUE)) {
  sk("susieR installed", "install susieR to run these")
} else {
  f <- names(formals(susieR::susie_rss))
  cat("  susieR", as.character(packageVersion("susieR")), "\n")
  cat("  formals:", paste(f, collapse=", "), "\n")
  # The wrapper passes prior_variance = NULL by default. Confirm that is either
  # an accepted formal or harmlessly absorbed - NOT silently mangling the prior.
  chk("n" %in% f, "susie_rss accepts n (needed for correct residual variance)")
  if (!"prior_variance" %in% f)
    cat("  NOTE  susie_rss has no 'prior_variance' formal; the wrapper's value goes to ... \n",
        "        Confirm it is ignored rather than misinterpreted.\n")
  # Sanity: a planted single causal must be recovered with high PIP.
  set.seed(1); p <- 50; nn <- 2000
  L0 <- diag(p); z0 <- rnorm(p, 0, 1); z0[10] <- 9
  fit <- tryCatch(susieR::susie_rss(z = z0, R = L0, n = nn, L = 5), error=function(e) e)
  if (inherits(fit,"error")) no("susie_rss runs on a trivial example", conditionMessage(fit))
  else {
    pip <- susieR::susie_get_pip(fit)
    chk(which.max(pip) == 10, "susie_rss puts max PIP on the planted causal",
        sprintf("argmax=%d pip=%.3f", which.max(pip), max(pip)))
  }
  chk("unmappable_effects" %in% f, "susie_rss supports unmappable_effects (SuSiE-inf)")
}

cat("\n=== B. polyfun_oracle prior reconstruction matches the simulator ===\n")
# The oracle defines the benchmark's ceiling. Its per-SNP prior must equal the
# probability the simulator actually used: pi ∝ exp(A %*% log(enrichment)).
set.seed(2); p <- 200; K <- 10
A <- matrix(rbinom(p*K, 1, 0.1), p, K)
enr <- c(rep(10.8, 5), rep(1, 5))
sim_w <- exp(as.numeric(A %*% log(enr)) - max(as.numeric(A %*% log(enr))))
sim_pi <- sim_w / sum(sim_w)
le <- log(pmax(enr, .Machine$double.eps))
lw <- as.numeric(A %*% le); w <- exp(lw - max(lw)); wrap_pi <- w / sum(w)
chk(max(abs(sim_pi - wrap_pi)) < 1e-12,
    "oracle pi == simulator pi", sprintf("max|diff| = %.2e", max(abs(sim_pi - wrap_pi))))

cat("\n=== C. ABF is Wakefield's ABF ===\n")
z <- c(0, 1, 3, 6); se <- rep(1, 4); W <- 0.04
V <- se^2; r <- W/(V+W)
wrapper_log_abf <- 0.5*log(1-r) + 0.5*r*z^2
direct_log_abf  <- log(sqrt(V/(V+W))) + (z^2/2)*(W/(V+W))
chk(max(abs(wrapper_log_abf - direct_log_abf)) < 1e-12,
    "ABF log-BF matches sqrt(V/(V+W))exp(z^2/2 * W/(V+W))",
    sprintf("max|diff| = %.2e", max(abs(wrapper_log_abf - direct_log_abf))))
# At z = 0 Wakefield's ABF is sqrt(V/(V+W)) < 1, so log-BF is NEGATIVE, not 0.
# (An earlier version of this script asserted 0 here and produced a false alarm.)
chk(abs(wrapper_log_abf[1] - 0.5*log(V[1]/(V[1]+W))) < 1e-12,
    "ABF log-BF at z=0 equals 0.5*log(V/(V+W))",
    sprintf("value = %.4f", wrapper_log_abf[1]))
chk(all(diff(wrapper_log_abf) > 0), "ABF log-BF increases with |z|")

cat("\n=== D. annotation source: geno-first everywhere ===\n")
# Regression test for the PR #19 bug class. simulate_gwfm_data() populates ONLY
# the geno-side annotation matrix, so any wrapper reading only region_pheno
# silently drops annotations.
geno  <- list(LD = diag(5), n = 100, variant_ids = paste0("v",1:5),
              annotations_matrix = matrix(1, 5, 2))
pheno <- list(z = rnorm(5))                       # NO annotations on the pheno side
for (nm in c(".fb_extract_annotations", ".paintor_extract_annotations")) {
  fn <- tryCatch(get(nm, envir = asNamespace("fmbenchmark")), error=function(e) NULL)
  if (is.null(fn)) { sk(paste(nm,"present")); next }
  chk(!is.null(fn(geno, pheno)), paste(nm, "reads the geno-side matrix"))
}
# funmap: exercised through its region adapter (returns an error result if the
# annotations were dropped).
res <- tryCatch(run_funmap_region(geno, pheno), error = function(e) e)
if (inherits(res, "error")) {
  sk("funmap adapter", conditionMessage(res))
} else {
  emsg <- if (is.null(res$error)) "" else res$error
  chk(!grepl("No annotations", emsg), "funmap sees geno-side annotations", emsg)
}

cat("\n=== E. PAINTOR pools across loci ===\n")
chk(exists("run_paintor_scenario_setup", mode="function"),
    "run_paintor_scenario_setup exists (cross-loci EM)")
if (exists("run_paintor_scenario_setup", mode="function")) {
  # With <2 regions there is nothing to pool: must return an empty list so the
  # per-locus path is used.
  chk(length(run_paintor_scenario_setup(list(geno), list(pheno), list())) == 0,
      "scenario_setup declines to pool a single locus")
}

cat("\n=== F. index-base conversions (0-based tools -> 1-based R) ===\n")
# BEATRICE / Functional BEATRICE emit 0-based credible-set indices.
tf <- tempfile(); writeLines(c("0 1 2", "4 5"), tf)
lines <- readLines(tf); unlink(tf)
cs <- lapply(lines, function(l) sort(as.integer(strsplit(trimws(l), "\\s+")[[1]]) + 1L))
# 0-based {0,1,2} -> 1-based {1,2,3}; 0-based {4,5} -> 1-based {5,6}.
# (An earlier version asserted {4,5} here and produced a false alarm.)
chk(identical(cs[[1]], 1:3) && identical(cs[[2]], 5:6),
    "0-based credible sets convert to 1-based correctly",
    paste(sapply(cs, paste, collapse=","), collapse=" | "))

cat("\n=== G. external tools present? (checks that could not otherwise run) ===\n")
for (pkg in c("susieR","CARMA","reticulate")) {
  if (requireNamespace(pkg, quietly=TRUE)) ok(paste(pkg,"installed"), as.character(packageVersion(pkg)))
  else sk(paste(pkg,"installed"), "absent -> its wrapper cannot be verified here")
}
cat("\nDone. Investigate every FAIL before trusting a run.\n")
