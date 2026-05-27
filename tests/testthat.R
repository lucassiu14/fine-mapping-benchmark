# Standard testthat entry point.
#
# Invoked by `R CMD check` automatically. During development, prefer
# `devtools::test()` or `testthat::test_dir("tests/testthat")` so that
# the package is loaded via pkgload::load_all() and the latest R/ files
# are picked up without requiring an install.

library(testthat)
library(fmbenchmark)

test_check("fmbenchmark")
