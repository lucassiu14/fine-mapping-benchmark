L1 <- readRDS("results/iter004/combined_fit_metrics.rds")
keep <- grepl("^job_00[1-9]_sparse_an", L1$job_dir)
sp <- L1[keep, , drop = FALSE]
# combined_fit_metrics.rds is the POST-merge L1: it already carries the grid
# columns and the pc/enrich recodes. iter004_report.R merges the grid itself, so
# feeding this back produces p_causal.x / p_causal.y and L1$p_causal goes NULL.
# Strip everything the report will re-derive, leaving raw collect output.
drop <- c("model", "p_causal", "annotation_type", "enrichment_fold",
          "n_ref", "n", "n_iter", "pc", "enrich")
sp <- sp[, setdiff(names(sp), drop), drop = FALSE]
cat("sparse-arm rows:", nrow(sp), " cols:", ncol(sp), "\n")
cat("grid cols removed:", paste(intersect(drop, names(L1)), collapse=", "), "\n")
stopifnot(!any(c("p_causal","model") %in% names(sp)))
dir.create("/tmp/l1_sparse", showWarnings = FALSE)
unlink("/tmp/l1_sparse/*.rds")
saveRDS(sp, "/tmp/l1_sparse/L1_sparsearm.rds")
