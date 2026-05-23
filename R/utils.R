# =============================================================================
# utils.R
#
# Shared utilities used across the benchmark.
# =============================================================================

# Null-coalescing operator (base R has no `??` equivalent). Sourced once here
# so wrappers can rely on it without redefining.
`%||%` <- function(x, y) if (!is.null(x)) x else y
