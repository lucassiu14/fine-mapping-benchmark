# =============================================================================
# inst/scripts/prepare_gwfm_vcfs.R
#
# Download VCF files for the 128 genome-wide benchmark regions defined in the
# bundled inst/extdata/gwfm_regions.csv.
#
# These regions are used by simulate_gwfm_data() and are spread across all
# 22 autosomes to provide a representative genome-wide sample of LD blocks.
# Regions are spaced to be approximately LD-independent of each other.
#
# Uses tabix to stream specific genomic windows from the 1000 Genomes Phase 3
# remote VCF files — no whole-chromosome download required.
#
# Source data:
#   1000 Genomes Project Phase 3 (GRCh37 / hg19)
#   2,504 individuals, all populations
#   http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
#
# Requirements:
#   - tabix and bgzip (install via: brew install htslib  OR  conda install -c bioconda htslib)
#
# Usage:
#   # From a source checkout (project root):
#   Rscript inst/scripts/prepare_gwfm_vcfs.R
#   # From an installed package:
#   Rscript "$(Rscript -e 'cat(system.file("scripts/prepare_gwfm_vcfs.R", package = "fmbenchmark"))')"
#
# Optional arguments (edit below):
#   VCF_DIR   where to save VCF files   (default: "data/gwfm_vcf")
#   REGIONS   path to regions CSV file  (default: bundled inst/extdata/gwfm_regions.csv)
#   OVERWRITE if TRUE, re-download existing files (default: FALSE)
#
# Output:
#   data/gwfm_vcf/gw001.vcf.gz   + gw001.vcf.gz.tbi
#   data/gwfm_vcf/gw002.vcf.gz   + gw002.vcf.gz.tbi
#   ...
#
# Each file is ~300 kb of genomic sequence (~300-600 SNPs before MAF filter).
# Total download: ~400 MB across all 128 regions.
# =============================================================================

# Locate a file bundled under inst/extdata/. Mirrors the package's internal
# fmb_extdata() helper, but inlined so this script runs standalone (no
# library(fmbenchmark) required). Resolution order: installed package, then
# the inst/extdata/ sibling of this script, then inst/extdata/ under cwd.
find_extdata <- function(filename) {
  p <- system.file("extdata", filename, package = "fmbenchmark")
  if (nzchar(p) && file.exists(p)) return(p)
  this_file  <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else getwd()
  candidates <- c(
    file.path(script_dir, "..", "extdata", filename),  # inst/scripts -> inst/extdata
    file.path("inst", "extdata", filename)              # run from project root
  )
  for (cand in candidates) if (file.exists(cand)) return(normalizePath(cand))
  stop("Could not locate bundled file 'inst/extdata/", filename, "'.", call. = FALSE)
}

VCF_DIR   <- "data/gwfm_vcf"
REGIONS   <- find_extdata("gwfm_regions.csv")
OVERWRITE <- FALSE

# 1000 Genomes Phase 3 remote VCF URL pattern (GRCh37, all populations)
REMOTE_VCF_PATTERN <-
  "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

# =============================================================================
# Setup
# =============================================================================

if (!dir.exists(VCF_DIR)) dir.create(VCF_DIR, recursive = TRUE)

# Check tabix and bgzip
tabix_path <- Sys.which("tabix")
bgzip_path <- Sys.which("bgzip")

if (tabix_path == "") {
  stop(
    "tabix not found on PATH.\n",
    "Install htslib via:\n",
    "  brew install htslib        (macOS)\n",
    "  conda install -c bioconda htslib  (conda)\n",
    "  sudo apt-get install tabix  (Linux)",
    call. = FALSE
  )
}
if (bgzip_path == "") {
  stop("bgzip not found on PATH. Install htslib (see above).", call. = FALSE)
}

# REGIONS is resolved by find_extdata() above, which already errors if the
# bundled CSV cannot be located, so no further existence check is needed.

regions <- read.csv(REGIONS, stringsAsFactors = FALSE)
stopifnot(all(c("region_id", "chrom", "start", "end") %in% names(regions)))

cat(sprintf(
  "Preparing %d genome-wide VCF regions → %s\n",
  nrow(regions), normalizePath(VCF_DIR)
))
cat(sprintf(
  "Estimated download: ~%.0f MB total (~3 MB per region)\n\n",
  nrow(regions) * 3
))

# =============================================================================
# Download loop
# =============================================================================

n_ok     <- 0L
n_skip   <- 0L
n_failed <- 0L
failed   <- character(0)

for (i in seq_len(nrow(regions))) {
  reg <- regions[i, ]
  rid <- reg$region_id
  chr <- as.character(reg$chrom)
  pos <- sprintf("%s:%d-%d", chr, reg$start, reg$end)

  out_vcf <- file.path(VCF_DIR, paste0(rid, ".vcf.gz"))
  out_tbi <- paste0(out_vcf, ".tbi")

  if (file.exists(out_vcf) && file.exists(out_tbi) && !OVERWRITE) {
    cat(sprintf("  [%3d/%d] %s %-32s SKIP (exists)\n",
                i, nrow(regions), rid, pos))
    n_skip <- n_skip + 1L
    next
  }

  cat(sprintf("  [%3d/%d] %s %-32s downloading...",
              i, nrow(regions), rid, pos))
  flush.console()

  remote_url <- sprintf(REMOTE_VCF_PATTERN, chr)

  # Stream the region with tabix, compress with bgzip
  # tabix -h includes the VCF header lines
  cmd_stream <- sprintf(
    "%s -h '%s' %s | %s -c > '%s'",
    tabix_path, remote_url, pos, bgzip_path, out_vcf
  )

  ret_stream <- system(cmd_stream, intern = FALSE, ignore.stderr = TRUE)

  if (ret_stream != 0 || !file.exists(out_vcf) || file.size(out_vcf) < 100) {
    cat(" FAILED (tabix stream)\n")
    if (file.exists(out_vcf)) file.remove(out_vcf)
    n_failed <- n_failed + 1L
    failed   <- c(failed, rid)
    next
  }

  # Index with tabix
  cmd_index <- sprintf("%s -p vcf '%s'", tabix_path, out_vcf)
  ret_index <- system(cmd_index, intern = FALSE, ignore.stderr = TRUE)

  if (ret_index != 0) {
    cat(" FAILED (tabix index)\n")
    file.remove(out_vcf)
    n_failed <- n_failed + 1L
    failed   <- c(failed, rid)
    next
  }

  size_kb <- round(file.size(out_vcf) / 1024)
  cat(sprintf(" OK (%d KB)\n", size_kb))
  n_ok <- n_ok + 1L
}

# =============================================================================
# Summary
# =============================================================================

cat(sprintf(
  "\nDone: %d downloaded, %d skipped (already exist), %d failed.\n",
  n_ok, n_skip, n_failed
))

if (length(failed) > 0) {
  cat("Failed regions:", paste(failed, collapse = ", "), "\n")
  cat("Re-run with OVERWRITE <- TRUE to retry failed regions.\n")
}

vcf_files <- list.files(VCF_DIR, pattern = "\\.vcf\\.gz$", full.names = TRUE)
vcf_files <- vcf_files[!grepl("\\.tbi$", vcf_files)]
if (length(vcf_files) > 0) {
  total_mb <- round(sum(file.size(vcf_files)) / 1e6, 1)
  cat(sprintf("Total VCF size on disk: %.1f MB\n", total_mb))
}

cat(sprintf(
  "\nPass vcf_dir = \"%s\" to simulate_gwfm_data() to use these regions.\n",
  VCF_DIR
))
