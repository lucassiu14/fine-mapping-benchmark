# =============================================================================
# scripts/prepare_vcfs.R
#
# Download VCF files for the 50 benchmark regions defined in data/regions.csv.
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
# Usage (from project root):
#   Rscript scripts/prepare_vcfs.R
#
# Optional arguments (set before sourcing or edit below):
#   VCF_DIR   where to save VCF files   (default: "data/vcf")
#   REGIONS   path to regions CSV file  (default: "data/regions.csv")
#   OVERWRITE if TRUE, re-download existing files (default: FALSE)
#
# Output:
#   data/vcf/r001.vcf.gz   + r001.vcf.gz.tbi
#   data/vcf/r002.vcf.gz   + r002.vcf.gz.tbi
#   ...
#
# Each file is ~300 kb of genomic sequence (~300-600 SNPs before MAF filter).
# Total download: ~150 MB across all 50 regions.
# =============================================================================

VCF_DIR   <- "data/vcf"
REGIONS   <- "data/regions.csv"
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

regions <- read.csv(REGIONS, stringsAsFactors = FALSE)
stopifnot(all(c("region_id", "chrom", "start", "end") %in% names(regions)))

cat(sprintf(
  "Preparing %d VCF regions → %s\n\n",
  nrow(regions), normalizePath(VCF_DIR)
))

# =============================================================================
# Download loop
# =============================================================================

n_ok     <- 0L
n_skip   <- 0L
n_failed <- 0L
failed   <- character(0)

for (i in seq_len(nrow(regions))) {
  reg  <- regions[i, ]
  rid  <- reg$region_id
  chr  <- as.character(reg$chrom)
  pos  <- sprintf("%s:%d-%d", chr, reg$start, reg$end)

  out_vcf <- file.path(VCF_DIR, paste0(rid, ".vcf.gz"))
  out_tbi <- paste0(out_vcf, ".tbi")

  if (file.exists(out_vcf) && file.exists(out_tbi) && !OVERWRITE) {
    cat(sprintf("  [%2d/%d] %s %-30s SKIP (exists)\n",
                i, nrow(regions), rid, pos))
    n_skip <- n_skip + 1L
    next
  }

  cat(sprintf("  [%2d/%d] %s %-30s downloading...",
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

total_mb <- round(sum(file.size(list.files(VCF_DIR, "*.vcf.gz", full.names = TRUE))) / 1e6, 1)
cat(sprintf("Total VCF size on disk: %.1f MB\n", total_mb))
cat(sprintf(
  "\nPass vcf_dir = \"%s\" to run_simulation() to use these regions.\n",
  VCF_DIR
))
