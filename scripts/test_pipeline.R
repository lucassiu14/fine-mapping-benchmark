# =============================================================================
# scripts/test_pipeline.R
#
# End-to-end test of the full benchmark pipeline:
#   1. Simulate genotypes + phenotypes (multiple regions, binary annotations)
#   2. Apply all five fine-mapping methods and save results
#
# Usage (from project root):
#   Rscript scripts/test_pipeline.R
# =============================================================================

source("R/simulate_genotypes.R")
source("R/simulate_phenotypes.R")
source("R/run_simulation.R")
source("R/run_methods.R")
source("R/evaluate.R")
source("R/plot_results.R")
source("R/wrappers/susie.R")
source("R/wrappers/susie_inf.R")
source("R/wrappers/abf.R")
source("R/wrappers/finemap.R")
source("R/wrappers/funmap.R")
source("R/wrappers/paintor.R")

# Path to FINEMAP binary (downloaded by setup_finemap()).
# On Apple Silicon without x86 libzstd this will fail gracefully.
FINEMAP_BIN <- tryCatch(
  setup_finemap(download = FALSE),
  error = function(e) "finemap"   # fall back to PATH search; will fail gracefully
)

# Path to PAINTOR binary. Install via: conda install -c bioconda paintor
# Falls back to "PAINTOR" (PATH search); will fail gracefully if not installed.
PAINTOR_BIN <- tryCatch(
  setup_paintor(),
  error = function(e) "PAINTOR"
)

OUTPUT_DIR <- "results/test_run"

# =============================================================================
# 1. Simulate
# =============================================================================

message("\n", strrep("=", 70))
message("STEP 1: Simulating data")
message(strrep("=", 70))

sim <- run_simulation(
  n_regions   = 3,
  n           = 300,
  p           = 100,
  n_iter      = 2,
  S           = c(1, 2),
  phi         = c(0.2, 0.4),
  model       = "sparse",
  annotations = "binary",
  n_annotations = 3,
  seed        = 42,
  save        = TRUE,
  output_dir  = OUTPUT_DIR,
  verbose     = TRUE
)

message(sprintf(
  "\nSimulation complete: %d scenarios x %d regions",
  length(sim$scenarios), length(sim$genotypes)
))

# Quick sanity checks
stopifnot(
  "3 regions"            = length(sim$genotypes) == 3,
  "annotations present"  = !is.null(sim$genotypes[[1]]$annotations_matrix),
  "scenarios have truth" = !is.null(sim$scenarios[[1]]$regions[[1]]$truth)
)
message("Sanity checks passed.")

# =============================================================================
# 2. Fine-mapping — run all methods
# =============================================================================

message("\n", strrep("=", 70))
message("STEP 2: Fine-mapping")
message(strrep("=", 70))

results <- run_methods(
  simulation  = sim,
  methods     = c("susie", "susie_inf", "abf", "finemap", "funmap", "paintor"),
  method_args = list(
    susie      = list(L = 10, coverage = 0.95),
    susie_inf  = list(L = 10, coverage = 0.95),
    abf        = list(prior_variance = 0.04),
    finemap    = list(finemap_path = FINEMAP_BIN, n_causal = 3),
    funmap     = list(python = "/opt/anaconda3/bin/python3", L = 10),
    paintor    = list(paintor_path = PAINTOR_BIN, max_causal = 2)
  ),
  save        = TRUE,
  output_dir  = OUTPUT_DIR,
  verbose     = TRUE
)

# =============================================================================
# 3. Summary
# =============================================================================

message("\n", strrep("=", 70))
message("SUMMARY")
message(strrep("=", 70))

for (m in results$methods_run) {
  res    <- results[[m]]
  n_ok   <- res$n_total - res$n_failed
  errors <- Filter(Negate(is.null), lapply(res$results, `[[`, "error"))
  cat(sprintf(
    "  %-12s %d/%d fits OK   %.1fs total",
    paste0(m, ":"), n_ok, res$n_total, res$total_runtime_seconds
  ))
  if (length(errors) > 0) {
    cat(sprintf("   [error: %s]", substr(errors[[1]], 1, 55)))
  }
  cat("\n")
}

message("\nAll outputs saved to: ", OUTPUT_DIR)

# =============================================================================
# 3. Evaluate
# =============================================================================

message("\n", strrep("=", 70))
message("STEP 3: Evaluation")
message(strrep("=", 70))

eval_out <- evaluate_methods(sim, results,
                             save       = TRUE,
                             output_dir = OUTPUT_DIR,
                             verbose    = TRUE)

# Print a compact metric table — one row per method
cat(sprintf(
  "\n%-12s  %6s  %6s  %6s  %6s  %6s  %6s\n",
  "Method", "AUPRC", "CS_cov", "CS_pow", "CS_med", "RT_mean", "n_fail"
))
cat(strrep("-", 62), "\n")

for (m in eval_out$methods_evaluated) {
  g <- eval_out[[m]]$global
  cat(sprintf(
    "%-12s  %6.3f  %6.3f  %6.3f  %6.1f  %6.1f  %6d\n",
    m,
    if (is.na(g$auprc))          NA_real_ else g$auprc,
    if (is.na(g$cs_coverage))    NA_real_ else g$cs_coverage,
    if (is.na(g$cs_power))       NA_real_ else g$cs_power,
    if (is.na(g$cs_size_median)) NA_real_ else g$cs_size_median,
    if (is.na(g$runtime_mean))   NA_real_ else g$runtime_mean,
    g$n_failed
  ))
}

# Stratified summary: AUPRC by S value for each method
cat("\nAUPRC by number of causal variants (S):\n")
all_S <- sort(unique(sapply(sim$scenarios, `[[`, "S")))
header <- sprintf("%-12s", "Method")
for (s in all_S) header <- paste0(header, sprintf("  S=%-3d", s))
cat(header, "\n", strrep("-", nchar(header)), "\n", sep = "")

for (m in eval_out$methods_evaluated) {
  row <- sprintf("%-12s", m)
  for (s in all_S) {
    key <- as.character(s)
    val <- eval_out[[m]]$by_S[[key]]$auprc
    row <- paste0(row, sprintf("  %-5s", if (is.null(val) || is.na(val)) "  NA " else sprintf("%.3f", val)))
  }
  cat(row, "\n")
}

message("\nEvaluation complete.")

# =============================================================================
# 4. Plot results
# =============================================================================

message("\n", strrep("=", 70))
message("STEP 4: Generating plots")
message(strrep("=", 70))

pdf_path <- file.path(OUTPUT_DIR, "results.pdf")
plot_results(eval_out, output_file = pdf_path, verbose = TRUE)
message("\nPlots saved to: ", pdf_path)
