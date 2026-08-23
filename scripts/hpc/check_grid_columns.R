#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/check_grid_columns.R
#
# Assert that a params_grid.csv carries every column run_benchmark_job.R reads.
#
# WHY THIS EXISTS. Iteration 005 introduced an alternative grid generator that
# emitted `enrichment_fold` but not `enrichment_values`. The worker reads the
# latter at run_benchmark_job.R:312:
#
#     enrichment_vec <- as.numeric(strsplit(job$enrichment_values, "\\|")[[1]])
#
# `job$enrichment_values` was NULL, strsplit() raised "non-character argument",
# and all 240 array tasks died there - five lines BEFORE the first dir.create(),
# so they left no output directory and no partial results to diagnose from. The
# design banner, the preflight and the method-set selection had all passed, so
# every check in place reported healthy right up to the failure.
#
# The missing column was the ONLY difference from the standing grid. So the
# invariant is simply: an alternative generator must emit at least what the
# standing generator emits. Extra columns are fine (Iteration 005 adds
# `relationship` and `n_informative`, both read behind is.null() guards).
#
#   Rscript scripts/hpc/check_grid_columns.R <candidate.csv> <canonical_generator.R>
#
# Exits non-zero and names the missing columns on failure.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("usage: check_grid_columns.R <candidate.csv> <canonical_generator.R>",
       call. = FALSE)
candidate <- args[1]; canonical_gen <- args[2]

if (!file.exists(candidate))     stop("grid not found: ", candidate, call. = FALSE)
if (!file.exists(canonical_gen)) stop("generator not found: ", canonical_gen, call. = FALSE)

cand_cols <- names(read.csv(candidate, nrows = 1L, stringsAsFactors = FALSE))

# Render the canonical grid to a throwaway path and read its header. Generating
# rather than hardcoding a list keeps this from drifting the moment the standing
# design changes - the same reason py_deps_check.py walks the tool sources
# instead of naming packages by hand.
tmp <- tempfile(fileext = ".csv")
on.exit(unlink(tmp), add = TRUE)
rc <- system2(file.path(R.home("bin"), "Rscript"), c(canonical_gen, tmp),
              stdout = FALSE, stderr = FALSE)
if (rc != 0L || !file.exists(tmp))
  stop("canonical generator failed to produce a grid: ", canonical_gen, call. = FALSE)
canon_cols <- names(read.csv(tmp, nrows = 1L, stringsAsFactors = FALSE))

missing <- setdiff(canon_cols, cand_cols)
extra   <- setdiff(cand_cols, canon_cols)

if (length(missing)) {
  cat(sprintf("GRID COLUMN CHECK FAILED\n  %s\n", candidate))
  cat(sprintf("  missing %d column(s) the worker reads: %s\n",
              length(missing), paste(missing, collapse = ", ")))
  cat("  compared against: ", canonical_gen, "\n", sep = "")
  quit(status = 1L)
}
cat(sprintf("grid column check ok - %d columns%s\n", length(cand_cols),
            if (length(extra))
              sprintf(" (%d beyond the standing design: %s)",
                      length(extra), paste(extra, collapse = ", ")) else ""))
