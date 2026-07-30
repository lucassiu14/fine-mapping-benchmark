#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/verify_new_methods.R
#
# Plant known causal variants in a synthetic locus and check that each method
# actually recovers them. Run ON THE CLUSTER after installing the dependencies
# (see docs/autoresearch/adding-methods.md).
#
#   module load R/4.5.2-gfbf-2025b GSL/2.8-GCC-14.3.0 Python/3.12.3-GCCcore-13.3.0
#   Rscript scripts/analysis/verify_new_methods.R
#
# WHY THIS EXISTS. CARMA was registered in Iteration 001 and never installed, so
# it returned NA for every fit across three iterations without anyone noticing -
# a silently-absent method is indistinguishable from a badly-performing one in
# the collected results. Every new method must clear this before it is given
# array time.
#
# Reports PASS / FAIL / SKIP. SKIP = dependency absent, so nothing was tested.
# =============================================================================
suppressWarnings(suppressMessages({
  if (requireNamespace("fmbenchmark", quietly = TRUE)) library(fmbenchmark)
  else if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION"))
    pkgload::load_all(".", quiet = TRUE)
  else stop("fmbenchmark not available")
}))
TOOLS <- normalizePath("~/tools", mustWork = FALSE)
PY    <- file.path(TOOLS, "py-venv-runner.sh")

# ---- synthetic locus with planted causals -----------------------------------
set.seed(101)
p <- 150; n <- 3000; causal <- c(37, 91)
k <- 30; Lf <- matrix(rnorm(p * k), p, k) / sqrt(k)
Sg <- Lf %*% t(Lf) + diag(p) * 0.6
d  <- sqrt(diag(Sg)); R <- Sg / outer(d, d); R <- 0.97 * R + 0.03 * diag(p)
bt <- numeric(p); bt[causal] <- c(0.13, -0.11)
z  <- as.numeric(sqrt(n) * (R %*% bt) + t(chol(R)) %*% rnorm(p))
se <- rep(1 / sqrt(n), p); bh <- z * se
geno <- list(LD = R, n = n, maf = runif(p, 0.05, 0.5),
             variant_ids = paste0("rs", seq_len(p)))
phen <- list(z = z, beta_hat = bh, se = se)

cat(sprintf("\nSynthetic locus: p=%d n=%d, causal at %s\n\n",
            p, n, paste(causal, collapse = ", ")))

report <- function(name, fit, dep_ok = TRUE, dep_msg = "") {
  if (!dep_ok) { cat(sprintf("  SKIP  %-14s %s\n", name, dep_msg)); return(invisible(NULL)) }
  if (!is.null(fit$error)) {
    # Show the TAIL of the message: for external tools the real cause (a Python
    # traceback, a linker error) is at the end, and an earlier version of this
    # truncated to the first 90 characters and hid it.
    msg <- gsub("[\r\n]+", " | ", fit$error)
    msg <- if (nchar(msg) > 300) paste0("...", substr(msg, nchar(msg) - 299, nchar(msg))) else msg
    cat(sprintf("  FAIL  %-14s %s\n", name, msg)); return(invisible(NULL))
  }
  pip <- fit$pip
  if (length(pip) != p)      { cat(sprintf("  FAIL  %-14s pip length %d, expected %d\n", name, length(pip), p)); return(invisible(NULL)) }
  if (all(is.na(pip)))       { cat(sprintf("  FAIL  %-14s all PIPs are NA (method absent or crashed silently)\n", name)); return(invisible(NULL)) }
  top <- order(pip, decreasing = TRUE)[1:2]
  hit <- sum(causal %in% top)
  verdict <- if (hit >= 1) "PASS" else "FAIL"
  cat(sprintf("  %s  %-14s sum(PIP)=%6.2f  top2=%s  PIP@causal=%s  %s\n",
              verdict, name, sum(pip, na.rm = TRUE), paste(top, collapse = ","),
              paste(sprintf("%.2f", pip[causal]), collapse = ","),
              if (hit == 2) "(both causals in top 2)" else if (hit == 1) "(1 of 2)" else "(MISSED)"))
}

cat("=== CARMA (LD-discrepancy outlier detection) ===\n")
if (!requireNamespace("CARMA", quietly = TRUE)) {
  report("carma", NULL, FALSE, "CARMA not installed - see adding-methods.md")
} else {
  report("carma", tryCatch(run_carma_region(geno, phen),
                           error = function(e) list(error = conditionMessage(e))))
}

# CAVIAR is intentionally NOT verified or run: it was dropped on cost 2026-07-29
# (enumeration is combinatorial in region size - too slow already at p=150, and
# the grid runs to p=1000). The wrapper is retained; see R/wrapper_caviar.R.
cat("\n=== CAVIAR ===\n  SKIP  caviar         dropped on cost; see R/wrapper_caviar.R\n")

cat("\n=== FiniMOM ===\n")
if (!requireNamespace("finimom", quietly = TRUE)) {
  report("finimom", NULL, FALSE, "finimom not installed - note the CXX_STD trap")
} else {
  report("finimom", tryCatch(run_finimom_region(geno, phen),
                             error = function(e) list(error = conditionMessage(e))))
}

cat("\n=== FINEMAP-inf ===\n")
fmi <- file.path(TOOLS, "fine-mapping-inf")
if (!file.exists(file.path(fmi, "run_fine_mapping.py"))) {
  report("finemap_inf", NULL, FALSE, "fine-mapping-inf not cloned")
} else {
  report("finemap_inf",
         tryCatch(run_finemap_inf_region(geno, phen, finemap_inf_dir = fmi,
                                         python = if (file.exists(PY)) PY else "python"),
                  error = function(e) list(error = conditionMessage(e))))
}

cat("\n=== incumbents (sanity: these should already work) ===\n")
for (nm in c("susie", "finemap", "abf")) {
  fn <- paste0("run_", nm, "_region")
  if (!exists(fn, mode = "function")) { report(nm, NULL, FALSE, "not found"); next }
  args <- list(geno, phen)
  if (nm == "finemap") args <- c(args, list(finemap_path = file.path(TOOLS,
      "finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64")))
  report(nm, tryCatch(do.call(match.fun(fn), args),
                      error = function(e) list(error = conditionMessage(e))))
}

cat("\nA method must PASS here before it is given array time.\n")
