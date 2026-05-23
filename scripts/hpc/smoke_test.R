# =============================================================================
# scripts/hpc/smoke_test.R
#
# Quick smoke test to verify that all 6 fine-mapping methods work correctly
# on your HPC environment before submitting the full benchmark.
#
# Runs one tiny simulation (1 region, p=100, n=200, 2 iterations) and applies
# every method. If a method fails, the error is printed and the test continues.
# A pass/fail summary is printed at the end.
#
# Usage (from project root):
#   Rscript scripts/hpc/smoke_test.R
#
# Edit the paths in the Configuration block below to match your HPC setup.
# =============================================================================

# =============================================================================
# Configuration — must match run_benchmark_job.R
# =============================================================================

PROJECT_ROOT  <- "."

PYTHON        <- "python"   # full path to conda-env python if needed

FB_DIR        <- file.path(PROJECT_ROOT, "BEATRICE_annot_sparse")

BEATRICE_DIR  <- file.path(PROJECT_ROOT, "alt_methods", "Beatrice-Finemapping")

PAINTOR_PATH  <- "PAINTOR"

VCF_DIR       <- file.path(PROJECT_ROOT, "data", "gwfm_vcf")

# =============================================================================
# Load project
# =============================================================================

r_files <- list.files(file.path(PROJECT_ROOT, "R"),
                       pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
invisible(lapply(r_files, source))

# =============================================================================
# Tiny simulation with annotations
# =============================================================================

cat("Running smoke-test simulation (p=100, n=200, 1 region, 2 iterations)...\n")

sim <- run_simulation(
  n_regions     = 1L,
  n             = 200L,
  p             = 100L,
  n_iter        = 2L,
  S             = c(1L, 2L),
  phi           = c(0.1, 0.3),
  model         = "sparse",
  annotations   = "binary",
  n_annotations = 3L,
  enrichment    = 3.0,
  vcf_dir       = if (dir.exists(VCF_DIR)) VCF_DIR else NULL,
  seed          = 42L,
  verbose       = FALSE
)
cat(sprintf("  Simulation OK: %d scenarios.\n\n", length(sim$scenarios)))

# =============================================================================
# Test each method individually
# =============================================================================

METHODS <- c("susie", "susie_inf", "beatrice", "funmap",
             "functional_beatrice", "paintor")

method_args <- list(
  susie               = list(L = 5L, coverage = 0.95),
  susie_inf           = list(L = 5L, coverage = 0.95),
  beatrice            = list(beatrice_dir = BEATRICE_DIR, python = PYTHON,
                             max_iter = 500L, n_caus = 5L, sparse_concrete = 20L),
  funmap              = list(python = PYTHON),
  functional_beatrice = list(beatrice_dir = FB_DIR, python = PYTHON,
                             max_iter = 500L, n_caus = 5L, sparse_concrete = 20L),
  paintor             = list(paintor_path = PAINTOR_PATH, max_causal = 3L)
)

results <- list()
status  <- character(length(METHODS))
names(status) <- METHODS
errors  <- character(length(METHODS))
names(errors) <- METHODS

for (m in METHODS) {
  cat(sprintf("Testing %-22s ... ", m))
  t0 <- proc.time()
  out <- tryCatch(
    run_methods(simulation  = sim,
                methods     = m,
                method_args = method_args[m],
                save        = FALSE,
                verbose     = FALSE),
    error = function(e) e
  )
  elapsed <- as.numeric((proc.time() - t0)["elapsed"])

  if (inherits(out, "error")) {
    status[m] <- "FAIL"
    errors[m] <- conditionMessage(out)
    cat(sprintf("FAIL  (%.1f s)\n", elapsed))
    cat(sprintf("  Error: %s\n", errors[m]))
  } else {
    # Check that at least one result has PIPs
    any_pip <- any(sapply(out[[m]]$results, function(r) {
      !is.null(r$pip) && !all(is.na(r$pip))
    }))
    if (any_pip) {
      status[m] <- "PASS"
      cat(sprintf("PASS  (%.1f s)\n", elapsed))
    } else {
      status[m] <- "WARN"
      errors[m] <- "All PIPs are NA"
      cat(sprintf("WARN  (%.1f s) — all PIPs NA\n", elapsed))
    }
  }
  results[[m]] <- out
}

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SMOKE TEST SUMMARY\n")
cat(strrep("=", 60), "\n")
for (m in METHODS) {
  cat(sprintf("  %-24s : %s\n", m, status[m]))
  if (status[m] %in% c("FAIL", "WARN") && nchar(errors[m]) > 0) {
    cat(sprintf("    → %s\n", errors[m]))
  }
}
cat(strrep("=", 60), "\n")

n_pass <- sum(status == "PASS")
n_fail <- sum(status %in% c("FAIL", "WARN"))
cat(sprintf("\n%d / %d methods passed.\n", n_pass, length(METHODS)))

if (n_fail > 0) {
  cat("\nFix the failing methods before submitting the full benchmark.\n")
  cat("Common issues:\n")
  cat("  beatrice / functional_beatrice: wrong --python path or conda env not activated\n")
  cat("  funmap  : install with pip install git+https://github.com/LeeHITsz/Funmap.git, check PYTHON path\n")
  cat("  paintor : binary not on PATH — set PAINTOR_PATH to full path\n")
  quit(save = "no", status = 1)
} else {
  cat("\nAll methods passed. Ready to submit the benchmark:\n")
  cat("  bash scripts/hpc/submit_benchmark.sh\n")
}
