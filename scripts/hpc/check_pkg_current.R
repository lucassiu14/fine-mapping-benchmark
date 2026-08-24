#!/usr/bin/env Rscript
# =============================================================================
# scripts/hpc/check_pkg_current.R
#
# Assert that the INSTALLED fmbenchmark is not older than the checked-out source.
#
# WHY THIS EXISTS. run_benchmark_job.R loads the package with
#
#     ok <- requireNamespace("fmbenchmark", quietly = TRUE)
#     if (!ok) pkgload::load_all() else library(fmbenchmark)
#
# so on any machine where the package IS installed - which is every compute node
# here - the INSTALLED build wins and the source tree is ignored entirely. A
# `git pull` that changes anything under R/ therefore has NO effect on the nodes
# until the package is reinstalled.
#
# Iteration 005 lost all 240 tasks of array 3895282 to exactly this. The source
# had gained `relationship` and `n_informative` arguments on run_simulation();
# the installed build had not. Every task passed preflight, decoded its grid,
# printed its design banner, spent ~20 minutes reaching the simulation call and
# then died on "unused argument". Testing against the source tree - where
# pkgload::load_all() reads the new code - passes and proves nothing about what
# the nodes run.
#
# The check is a timestamp comparison rather than a signature test, so it
# catches ANY source change, not just new formals.
#
#   Rscript scripts/hpc/check_pkg_current.R [project_root]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L && nzchar(args[1])) args[1] else "."

src_files <- list.files(file.path(root, "R"), pattern = "[.][Rr]$", full.names = TRUE)
if (!length(src_files)) stop("no R/ sources found under ", root, call. = FALSE)

inst <- tryCatch(find.package("fmbenchmark"), error = function(e) NA_character_)
if (is.na(inst)) {
  # Not installed: the worker falls back to pkgload::load_all(), which reads the
  # source directly. That is always current, so there is nothing to check.
  cat("fmbenchmark is not installed - the worker will load_all() from source. ok\n")
  quit(status = 0L)
}

# Meta/package.rds is written at install time and is the most reliable stamp of
# when the build was produced; DESCRIPTION is a fallback for unusual layouts.
stamp_file <- file.path(inst, "Meta", "package.rds")
if (!file.exists(stamp_file)) stamp_file <- file.path(inst, "DESCRIPTION")
if (!file.exists(stamp_file))
  stop("cannot determine install time for fmbenchmark at ", inst, call. = FALSE)

built  <- file.mtime(stamp_file)
newest <- max(file.mtime(src_files))
stale  <- src_files[file.mtime(src_files) > built]

cat(sprintf("installed : %s  (built %s)\n", inst, format(built, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("source    : %d files, newest %s\n", length(src_files),
            format(newest, "%Y-%m-%d %H:%M:%S")))

if (length(stale)) {
  cat("\nPACKAGE IS STALE - the compute nodes would run the OLD code.\n")
  cat(sprintf("  %d source file(s) newer than the installed build:\n", length(stale)))
  for (f in head(sort(basename(stale)), 12L)) cat("    ", f, "\n", sep = "")
  if (length(stale) > 12L) cat(sprintf("    ... and %d more\n", length(stale) - 12L))
  cat("\nReinstall before submitting:\n")
  cat("  R -e 'install.packages(\".\", repos = NULL, type = \"source\")'\n")
  quit(status = 1L)
}
cat("package check ok - installed build is at least as new as the source\n")
