# =============================================================================
# utils.R
#
# Shared utilities used across the benchmark.
# =============================================================================

# Null-coalescing operator (base R has no `??` equivalent). Sourced once here
# so wrappers can rely on it without redefining.
`%||%` <- function(x, y) if (!is.null(x)) x else y


# -----------------------------------------------------------------------------
# fmb_extdata(): locate a file shipped under inst/extdata/
#
# When the package is installed, system.file() returns the install-tree path.
# During development (pkgload::load_all() / devtools::test() / running tests
# from tests/testthat/), system.file() can return "" because the package is
# not formally installed. In that case we fall back to candidate paths under
# inst/extdata/ relative to a few plausible working directories: the project
# root, one level below, and the testthat directory.
#
# Internal helper - not exported.
# -----------------------------------------------------------------------------
fmb_extdata <- function(filename) {
  path <- system.file("extdata", filename, package = "fmbenchmark")
  if (nzchar(path) && file.exists(path)) {
    return(path)
  }
  candidates <- c(
    file.path("inst", "extdata", filename),
    file.path("..", "inst", "extdata", filename),
    file.path("..", "..", "inst", "extdata", filename)
  )
  for (cand in candidates) {
    if (file.exists(cand)) {
      return(normalizePath(cand, mustWork = TRUE))
    }
  }
  stop(
    "Could not locate bundled file 'inst/extdata/", filename, "'. ",
    "If you installed the package via devtools::install() or remotes::install_github(), ",
    "this is a bug - please report it. If you are running from a source checkout, ",
    "ensure the file exists under inst/extdata/.",
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# Resolve a Python executable path for system2().
#
# Expands "~", falls back to the PATH, and makes the result absolute - but
# deliberately does NOT call normalizePath(), which follows symlinks. A
# virtualenv's bin/python is a symlink to the base interpreter, and resolving it
# bypasses the environment entirely: the base interpreter then cannot import
# torch and every call fails with status 1. Conda envs behave the same way.
# Returns the path unchanged when it cannot be found, so the caller's own
# error message is what the user sees.
#
# Internal helper - not exported.
# -----------------------------------------------------------------------------
fmb_resolve_python <- function(python, use_path = FALSE) {
  p <- path.expand(python)
  if (!file.exists(p)) {
    # Not a path on disk. Only look it up on PATH when the caller asks
    # (setup_*() does; the run_*() wrappers historically pass the string
    # straight through and record it verbatim in $params).
    if (!use_path) return(python)
    found <- Sys.which(python)
    if (nchar(found) == 0L) return(python)
    p <- unname(found)
  }
  if (!startsWith(p, "/")) p <- file.path(getwd(), p)
  p
}
