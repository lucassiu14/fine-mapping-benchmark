# Resolve a writable genetic-map cache directory for the test suite.
#
# Tests must work from several working directories: the project root
# (devtools::test()), tests/testthat/ (testthat::test_file()), and a
# detached R CMD check tempdir (where the source-tree data/ cache is
# absent because data/ is .Rbuildignore'd).
#
# Prefer an existing project-level cache so local runs reuse the already
# downloaded HapMap maps and avoid re-downloading. Otherwise fall back to a
# session tempdir, into which sim1000G downloads the ~1 MB chromosome map on
# first use. tempdir() always exists and is writable, so readGeneticMap()
# always receives a valid path.
fmb_test_map_dir <- function() {
  for (d in c("../../data/genetic_maps", "data/genetic_maps")) {
    if (dir.exists(d)) return(normalizePath(d))
  }
  d <- file.path(tempdir(), "fmb_genetic_maps")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}
