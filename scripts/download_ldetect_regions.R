# =============================================================================
# scripts/download_ldetect_regions.R
#
# Download the LDetect LD block partition and (optionally) the corresponding
# 1000 Genomes VCF slices for use with simulate_gwfm_data().
#
# LDetect (Berisa & Pickrell, 2016, Bioinformatics) partitions the autosome
# into ~1,600-1,700 approximately independent LD blocks, estimated separately
# for three continental ancestry groups:
#
#   EUR  European       ~1,703 blocks
#   AFR  African        ~1,445 blocks
#   ASN  East Asian     ~1,647 blocks
#
# Source: https://bitbucket.org/nygcresearch/ldetect-data
# Paper:  https://doi.org/10.1093/bioinformatics/btv546
#
# ---------------------------------------------------------------------------
# STEP 1 (this script): Download LDetect blocks → data/gwfm_regions_ldetect_{POP}.csv
#   - Fast (~seconds), ~100 KB download per population
#   - No large files stored on disk
#
# STEP 2 (optional, separate): Download VCF slices for those blocks
#   - Run with DOWNLOAD_VCFS <- TRUE below, OR run manually afterwards
#   - ~3 MB per block → ~5 GB total for all EUR blocks
#   - Uses tabix to stream only the required genomic window (no full chr download)
#   - ONLY run this if you have sufficient disk space and a stable connection
#
# ---------------------------------------------------------------------------
# Usage (from project root):
#   Rscript scripts/download_ldetect_regions.R
#
# Configuration (edit below):
#   POPULATION    ancestry group: "EUR", "AFR", or "ASN"
#   OUT_CSV       where to save the converted region CSV
#   DOWNLOAD_VCFS if TRUE, also download VCF slices (requires tabix/bgzip)
#   VCF_DIR       where to save VCF files (only used if DOWNLOAD_VCFS = TRUE)
#   MIN_BLOCK_KB  minimum block size to retain (removes very small blocks)
#   OVERWRITE     re-download even if output files already exist
# =============================================================================

POPULATION    <- "EUR"
OUT_CSV       <- file.path("data", sprintf("gwfm_regions_ldetect_%s.csv", POPULATION))
DOWNLOAD_VCFS <- FALSE          # set TRUE to also download VCF slices
VCF_DIR       <- file.path("data", sprintf("gwfm_vcf_ldetect_%s", POPULATION))
MIN_BLOCK_KB  <- 100            # drop blocks smaller than 100 kb
OVERWRITE     <- FALSE

# =============================================================================
# LDetect data URLs
# =============================================================================

LDETECT_BASE <-
  "https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc"

LDETECT_URLS <- list(
  EUR = file.path(LDETECT_BASE, "EUR", "fourier_ls-all.bed"),
  AFR = file.path(LDETECT_BASE, "AFR", "fourier_ls-all.bed"),
  ASN = file.path(LDETECT_BASE, "ASN", "fourier_ls-all.bed")
)

# 1000 Genomes Phase 3 remote VCF pattern (GRCh37, all populations)
REMOTE_VCF_PATTERN <-
  "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

# =============================================================================
# Validate inputs
# =============================================================================

if (!POPULATION %in% names(LDETECT_URLS)) {
  stop("POPULATION must be one of: ", paste(names(LDETECT_URLS), collapse = ", "),
       call. = FALSE)
}

if (!dir.exists("data")) dir.create("data", recursive = TRUE)

# =============================================================================
# STEP 1: Download and convert LDetect blocks
# =============================================================================

cat(sprintf("Downloading LDetect %s blocks from Bitbucket...\n", POPULATION))

url <- LDETECT_URLS[[POPULATION]]

raw <- tryCatch(
  readLines(url(url)),
  error = function(e) {
    stop(
      "Failed to download LDetect data.\n",
      "URL: ", url, "\n",
      "Error: ", conditionMessage(e), "\n",
      "Check your internet connection or download manually from:\n",
      "  https://bitbucket.org/nygcresearch/ldetect-data\n",
      call. = FALSE
    )
  }
)

# LDetect BED format: tab-separated, columns chr / start / stop
# Header line starts with "chr"
header_row <- grep("^chr\\s+start", raw, ignore.case = TRUE)
if (length(header_row) > 0) {
  raw <- raw[-header_row]
}

parsed <- do.call(rbind, strsplit(trimws(raw[nchar(trimws(raw)) > 0]), "\\s+"))
ldetect_df <- data.frame(
  chr   = parsed[, 1],
  start = as.integer(parsed[, 2]),
  stop  = as.integer(parsed[, 3]),
  stringsAsFactors = FALSE
)

cat(sprintf("  Downloaded %d raw blocks.\n", nrow(ldetect_df)))

# Strip "chr" prefix from chromosome column to match our convention (e.g. "1" not "chr1")
ldetect_df$chr <- sub("^chr", "", ldetect_df$chr)

# Keep only autosomes (1-22)
autosomes <- as.character(1:22)
ldetect_df <- ldetect_df[ldetect_df$chr %in% autosomes, ]
cat(sprintf("  %d blocks on autosomes.\n", nrow(ldetect_df)))

# Filter by minimum block size
block_kb <- (ldetect_df$stop - ldetect_df$start) / 1000
ldetect_df <- ldetect_df[block_kb >= MIN_BLOCK_KB, ]
cat(sprintf("  %d blocks >= %d kb.\n", nrow(ldetect_df), MIN_BLOCK_KB))

# Convert to our region CSV format
# We use the centre 300 kb of each block as the simulation window,
# consistent with the bundled gwfm_regions.csv.
# For blocks < 300 kb (after filtering), use the full block.

WINDOW_BP <- 300000L

ldetect_df$mid      <- as.integer((ldetect_df$start + ldetect_df$stop) / 2)
ldetect_df$reg_start <- pmax(ldetect_df$start, ldetect_df$mid - WINDOW_BP %/% 2L)
ldetect_df$reg_end   <- ldetect_df$reg_start + WINDOW_BP

# Clamp to block boundaries
ldetect_df$reg_start <- pmax(ldetect_df$reg_start, ldetect_df$start)
ldetect_df$reg_end   <- pmin(ldetect_df$reg_end,   ldetect_df$stop)

# Build region IDs: ld{POP}_{chr}_{index_within_chr}
region_df <- data.frame(
  region_id = character(nrow(ldetect_df)),
  chrom     = as.integer(ldetect_df$chr),
  start     = ldetect_df$reg_start,
  end       = ldetect_df$reg_end,
  block_start = ldetect_df$start,
  block_end   = ldetect_df$stop,
  notes     = sprintf("LDetect %s block (full block: %s:%d-%d)",
                      POPULATION, ldetect_df$chr,
                      ldetect_df$start, ldetect_df$stop),
  stringsAsFactors = FALSE
)

# Sort by chromosome then position
region_df <- region_df[order(region_df$chrom, region_df$start), ]

# Assign sequential region IDs
chr_counters <- integer(22)
for (k in seq_len(nrow(region_df))) {
  chr_k <- region_df$chrom[k]
  chr_counters[chr_k] <- chr_counters[chr_k] + 1L
  region_df$region_id[k] <- sprintf("ld%s_%02d_%04d",
                                     POPULATION, chr_k, chr_counters[chr_k])
}

cat(sprintf("\nConverted to %d regions (300 kb windows centred in each block).\n",
            nrow(region_df)))
cat(sprintf("Chromosomes covered: %s\n",
            paste(sort(unique(region_df$chrom)), collapse = ", ")))

# Write CSV
write.csv(region_df, OUT_CSV, row.names = FALSE, quote = FALSE)
cat(sprintf("Region file saved to: %s\n", OUT_CSV))

# =============================================================================
# Coverage summary
# =============================================================================

cat("\n--- Coverage summary ---\n")
cat(sprintf("  Total regions:   %d\n", nrow(region_df)))
cat(sprintf("  Chromosomes:     %d\n", length(unique(region_df$chrom))))
cat(sprintf("  Bundled set:     128 regions (~7%% of LDetect blocks)\n"))
cat(sprintf("  This set:        %d regions (~%.0f%% of LDetect blocks)\n",
            nrow(region_df), 100 * nrow(region_df) / nrow(ldetect_df)))
cat(sprintf("  Est. SNPs @ p=200: ~%s\n",
            formatC(nrow(region_df) * 200, format = "d", big.mark = ",")))
cat(sprintf("  Est. VCF download: ~%.1f GB (if DOWNLOAD_VCFS = TRUE)\n",
            nrow(region_df) * 3 / 1000))

cat("\nTo use these regions with simulate_gwfm_data():\n")
cat(sprintf('  regions <- read.csv("%s")\n', OUT_CSV))
cat('  sim <- simulate_gwfm_data(n = 1000, regions = regions,\n')
cat(sprintf('                             vcf_dir = "%s", ...)\n', VCF_DIR))

# =============================================================================
# STEP 2 (optional): Download VCF slices
# =============================================================================

if (!DOWNLOAD_VCFS) {
  cat(sprintf(
    "\nVCF download skipped (DOWNLOAD_VCFS = FALSE).\n",
    "To download VCF files for these regions, set DOWNLOAD_VCFS <- TRUE and re-run,\n",
    "or run scripts/prepare_gwfm_vcfs.R after pointing REGIONS at:\n  %s\n",
    OUT_CSV
  ))
  cat(sprintf(
    "Estimated disk space required: ~%.1f GB for all %d regions.\n",
    nrow(region_df) * 3 / 1000, nrow(region_df)
  ))
  quit(save = "no")
}

# --- VCF download -----------------------------------------------------------

cat(sprintf(
  "\n--- Downloading VCF slices for %d regions ---\n",
  nrow(region_df)
))
cat(sprintf(
  "WARNING: This will download approximately %.1f GB of data.\n",
  nrow(region_df) * 3 / 1000
))
cat("Each region is streamed directly from 1000 Genomes FTP via tabix.\n")
cat("Press Ctrl+C to cancel.\n\n")

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

if (!dir.exists(VCF_DIR)) dir.create(VCF_DIR, recursive = TRUE)

n_ok     <- 0L
n_skip   <- 0L
n_failed <- 0L
failed   <- character(0)

for (i in seq_len(nrow(region_df))) {
  reg     <- region_df[i, ]
  rid     <- reg$region_id
  chr     <- as.character(reg$chrom)
  pos     <- sprintf("%s:%d-%d", chr, reg$start, reg$end)
  out_vcf <- file.path(VCF_DIR, paste0(rid, ".vcf.gz"))
  out_tbi <- paste0(out_vcf, ".tbi")

  if (file.exists(out_vcf) && file.exists(out_tbi) && !OVERWRITE) {
    cat(sprintf("  [%4d/%d] %s  SKIP\n", i, nrow(region_df), rid))
    n_skip <- n_skip + 1L
    next
  }

  cat(sprintf("  [%4d/%d] %s (%s)  ...", i, nrow(region_df), rid, pos))
  flush.console()

  remote_url <- sprintf(REMOTE_VCF_PATTERN, chr)
  cmd_stream <- sprintf(
    "%s -h '%s' %s | %s -c > '%s'",
    tabix_path, remote_url, pos, bgzip_path, out_vcf
  )
  ret_stream <- system(cmd_stream, intern = FALSE, ignore.stderr = TRUE)

  if (ret_stream != 0 || !file.exists(out_vcf) || file.size(out_vcf) < 100) {
    cat(" FAILED\n")
    if (file.exists(out_vcf)) file.remove(out_vcf)
    n_failed <- n_failed + 1L
    failed   <- c(failed, rid)
    next
  }

  cmd_index <- sprintf("%s -p vcf '%s'", tabix_path, out_vcf)
  ret_index <- system(cmd_index, intern = FALSE, ignore.stderr = TRUE)

  if (ret_index != 0) {
    cat(" FAILED (index)\n")
    file.remove(out_vcf)
    n_failed <- n_failed + 1L
    failed   <- c(failed, rid)
    next
  }

  size_kb <- round(file.size(out_vcf) / 1024)
  cat(sprintf(" OK (%d KB)\n", size_kb))
  n_ok <- n_ok + 1L
}

cat(sprintf(
  "\nDone: %d downloaded, %d skipped, %d failed.\n",
  n_ok, n_skip, n_failed
))
if (length(failed) > 0) {
  cat("Failed:", paste(failed, collapse = ", "), "\n")
  cat("Re-run with OVERWRITE <- TRUE to retry.\n")
}

vcf_files <- list.files(VCF_DIR, "\\.vcf\\.gz$", full.names = TRUE)
vcf_files <- vcf_files[!grepl("\\.tbi$", vcf_files)]
if (length(vcf_files) > 0) {
  cat(sprintf(
    "Total VCF size on disk: %.1f GB\n",
    sum(file.size(vcf_files)) / 1e9
  ))
}
