# =============================================================================
# download_1kg_vcf.R
#
# Helper functions for downloading and extracting genomic regions from
# 1000 Genomes Phase 3 VCF files for use with simulate_genotypes().
#
# Requires: tabix (from htslib) installed and on PATH.
# =============================================================================

# --- Constants ----------------------------------------------------------------

#' @keywords internal
KG_BASE_URL <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

#' @keywords internal
KG_VCF_PATTERN <- "ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

#' @keywords internal
KG_PED_URL <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped"

#' @keywords internal
KG_EUR_POPS <- c("CEU", "TSI", "GBR", "FIN", "IBS")

#' @keywords internal
KG_SUPERPOPS <- list(
  EUR = c("CEU", "TSI", "GBR", "FIN", "IBS"),
  EAS = c("CHB", "JPT", "CHS", "CDX", "KHV"),
  AFR = c("YRI", "LWK", "GWD", "MSL", "ESN", "ACB", "ASW"),
  AMR = c("MXL", "PUR", "CLM", "PEL"),
  SAS = c("GIH", "PJL", "BEB", "STU", "ITU")
)


# --- Pre-defined genomic regions ----------------------------------------------

#' A curated set of independent genomic regions for benchmarking
#'
#' These regions are ~1 Mb each, drawn from different chromosomes to ensure
#' independence. They are chosen to contain moderate variant density and
#' diverse LD patterns. Coordinates are GRCh37/hg19.
#'
#' @return A data.frame with columns: chr, start, end, label
#'
#' @examples
#' \dontrun{
#' regions <- default_benchmark_regions()
#' head(regions)
#' }
#'
#' @export
default_benchmark_regions <- function() {
  data.frame(
    chr   = c( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
              11, 12, 13, 14, 15, 16, 17, 18, 19, 20),
    start = c(10000000, 30000000, 50000000, 40000000, 60000000,
              25000000, 70000000, 20000000, 35000000, 45000000,
              55000000, 15000000, 40000000, 30000000, 25000000,
              10000000, 35000000, 20000000, 5000000,  30000000),
    end   = c(11000000, 31000000, 51000000, 41000000, 61000000,
              26000000, 71000000, 21000000, 36000000, 46000000,
              56000000, 16000000, 41000000, 31000000, 26000000,
              11000000, 36000000, 21000000, 6000000,  31000000),
    label = paste0("region_chr", c(1:20)),
    stringsAsFactors = FALSE
  )
}


# --- Check for tabix ----------------------------------------------------------

#' Check whether tabix (from htslib) is available on the system PATH
#'
#' @return Logical. TRUE if tabix is found, FALSE otherwise.
#' @keywords internal
check_tabix <- function() {
  result <- tryCatch(
    system2("tabix", "--version", stdout = TRUE, stderr = TRUE),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  !is.null(result)
}

#' Check whether bcftools is available on the system PATH
#'
#' @return Logical. TRUE if bcftools is found, FALSE otherwise.
#' @keywords internal
check_bcftools <- function() {
  result <- tryCatch(
    system2("bcftools", "--version", stdout = TRUE, stderr = TRUE),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  !is.null(result)
}


# --- Download 1000 Genomes PED file -------------------------------------------

#' Download the 1000 Genomes Phase 3 sample information PED file
#'
#' Used to identify sample IDs belonging to specific populations or
#' super-populations for subsetting VCF files.
#'
#' @param output_dir Directory to save the PED file. Default: tempdir().
#' @param force Logical. Re-download even if file already exists. Default: FALSE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return Path to the downloaded PED file.
#'
#' @export
download_1kg_ped <- function(output_dir = tempdir(),
                             force = FALSE,
                             verbose = TRUE) {

  ped_file <- file.path(output_dir, "20130606_g1k.ped")

  if (file.exists(ped_file) && !force) {
    if (verbose) message("PED file already exists: ", ped_file)
    return(ped_file)
  }

  if (verbose) message("Downloading 1000 Genomes PED file...")
  download.file(KG_PED_URL, ped_file, mode = "w", quiet = !verbose)

  if (!file.exists(ped_file)) {
    stop("Failed to download PED file.", call. = FALSE)
  }

  if (verbose) message("Saved to: ", ped_file)
  ped_file
}


#' Get sample IDs for a given population or super-population
#'
#' @param population Character. Either a super-population code (e.g. "EUR",
#'   "EAS", "AFR", "AMR", "SAS") or a specific population code (e.g. "CEU",
#'   "GBR", "YRI"). Multiple codes can be provided. Default: "EUR".
#' @param ped_file Path to the 1000 Genomes PED file. If NULL, downloads it
#'   automatically. Default: NULL.
#' @param verbose Logical. Default: TRUE.
#'
#' @return Character vector of sample IDs.
#'
#' @examples
#' \dontrun{
#' # European samples
#' eur_ids <- get_1kg_sample_ids("EUR")
#'
#' # Specific populations
#' ceu_gbr_ids <- get_1kg_sample_ids(c("CEU", "GBR"))
#' }
#'
#' @export
get_1kg_sample_ids <- function(population = "EUR",
                               ped_file = NULL,
                               verbose = TRUE) {

  # Download PED if not provided
  if (is.null(ped_file)) {
    ped_file <- download_1kg_ped(verbose = verbose)
  }

  ped <- read.table(ped_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, comment.char = "")

  # Expand super-populations
  pop_codes <- character(0)
  for (p in population) {
    if (toupper(p) %in% names(KG_SUPERPOPS)) {
      pop_codes <- c(pop_codes, KG_SUPERPOPS[[toupper(p)]])
    } else {
      pop_codes <- c(pop_codes, toupper(p))
    }
  }
  pop_codes <- unique(pop_codes)

  # Filter
  ids <- ped$Individual.ID[ped$Population %in% pop_codes]

  if (length(ids) == 0) {
    stop(
      "No samples found for population(s): ",
      paste(population, collapse = ", "),
      "\nAvailable populations: ",
      paste(unique(ped$Population), collapse = ", "),
      call. = FALSE
    )
  }

  if (verbose) {
    message(
      sprintf(
        "Found %d samples for population(s): %s",
        length(ids), paste(pop_codes, collapse = ", ")
      )
    )
  }

  ids
}


# --- Download VCF region ------------------------------------------------------

#' Download a genomic region from 1000 Genomes Phase 3 VCF files
#'
#' Uses tabix to remotely extract a specific genomic region from the
#' 1000 Genomes FTP server. This avoids downloading entire chromosome
#' VCF files (which are several GB each).
#'
#' @param chr Integer or character. Chromosome number (1-22).
#' @param start Integer. Start position (bp, GRCh37/hg19).
#' @param end Integer. End position (bp, GRCh37/hg19).
#' @param output_dir Directory to save the extracted VCF. Default: tempdir().
#' @param population Character. Population(s) to subset to (e.g. "EUR",
#'   "CEU"). If NULL, all 2504 samples are included. Default: "EUR".
#' @param output_filename Character or NULL. Custom filename for the output
#'   VCF. If NULL, a default name is generated. Default: NULL.
#' @param force Logical. Re-download even if file already exists. Default: FALSE.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return Path to the downloaded and (optionally) population-filtered VCF file.
#'
#' @details
#' This function requires \code{tabix} (from htslib) to be installed and
#' available on the system PATH. If you also want to filter by population,
#' \code{bcftools} is additionally required.
#'
#' Coordinates are GRCh37/hg19 (the reference build used by 1000 Genomes
#' Phase 3).
#'
#' The function downloads the region directly from the EBI HTTP mirror of the
#' 1000 Genomes FTP site. No full chromosome VCF download is needed.
#'
#' @examples
#' \dontrun{
#' # Download a 1 Mb region on chr2 for European samples
#' vcf_path <- download_1kg_region(
#'   chr = 2, start = 30000000, end = 31000000,
#'   population = "EUR"
#' )
#'
#' # Use it with simulate_genotypes()
#' geno <- simulate_genotypes(
#'   n_regions = 1, n = 500, p = 300,
#'   vcf_files = vcf_path
#' )
#' }
#'
#' @export
download_1kg_region <- function(chr,
                                start,
                                end,
                                output_dir = tempdir(),
                                population = "EUR",
                                output_filename = NULL,
                                force = FALSE,
                                verbose = TRUE) {

  # --- Validate inputs --------------------------------------------------------

  chr <- as.character(chr)
  if (!chr %in% as.character(1:22)) {
    stop("chr must be an integer between 1 and 22.", call. = FALSE)
  }

  stopifnot(
    "start must be a positive integer" =
      is.numeric(start) && length(start) == 1 && start >= 1,
    "end must be a positive integer greater than start" =
      is.numeric(end) && length(end) == 1 && end > start
  )

  start <- as.integer(start)
  end   <- as.integer(end)

  # --- Check tools ------------------------------------------------------------

  has_tabix <- check_tabix()
  has_bcftools <- check_bcftools()

  if (!has_tabix && !has_bcftools) {
    stop(
      "Neither tabix nor bcftools found on PATH.\n",
      "Install htslib (provides tabix) and/or bcftools:\n",
      "  - conda install -c bioconda htslib bcftools\n",
      "  - brew install htslib bcftools  (macOS)\n",
      "  - apt-get install tabix bcftools  (Ubuntu/Debian)",
      call. = FALSE
    )
  }

  # --- Construct paths --------------------------------------------------------

  vcf_url <- sprintf("%s/%s", KG_BASE_URL, sprintf(KG_VCF_PATTERN, chr))
  region_str <- sprintf("%s:%d-%d", chr, start, end)

  if (is.null(output_filename)) {
    output_filename <- sprintf("1kg_chr%s_%d_%d.vcf.gz", chr, start, end)
  }
  output_path <- file.path(output_dir, output_filename)

  if (file.exists(output_path) && !force) {
    if (verbose) message("VCF already exists: ", output_path)
    return(output_path)
  }

  # --- Download region --------------------------------------------------------

  if (verbose) {
    message(sprintf("Downloading region %s from 1000 Genomes Phase 3...", region_str))
    message(sprintf("Source: %s", vcf_url))
  }

  # Strategy: use bcftools if we need population filtering, tabix otherwise
  if (!is.null(population) && has_bcftools) {

    # Get sample IDs for the population
    sample_ids <- get_1kg_sample_ids(population, verbose = verbose)

    # Write sample IDs to a temp file
    samples_file <- tempfile(fileext = ".txt")
    writeLines(sample_ids, samples_file)
    on.exit(unlink(samples_file), add = TRUE)

    # bcftools view with region and sample filtering, output bgzipped
    cmd <- sprintf(
      'bcftools view -r %s -S %s --force-samples -Oz -o %s %s',
      region_str,
      samples_file,
      shQuote(output_path),
      vcf_url
    )

    if (verbose) message("Running: ", cmd)
    exit_code <- system(cmd, ignore.stdout = !verbose, ignore.stderr = !verbose)

    if (exit_code != 0) {
      stop(
        "bcftools command failed (exit code ", exit_code, ").\n",
        "Check that the URL is accessible and bcftools is working.",
        call. = FALSE
      )
    }

  } else if (has_tabix) {

    # tabix can slice a region from a remote file
    # Output is uncompressed VCF; we'll need to bgzip it
    tmp_vcf <- tempfile(fileext = ".vcf")
    on.exit(unlink(tmp_vcf), add = TRUE)

    # First get the header
    header_cmd <- sprintf(
      'tabix -H %s > %s',
      vcf_url,
      shQuote(tmp_vcf)
    )

    if (verbose) message("Fetching VCF header...")
    system(header_cmd, ignore.stdout = !verbose, ignore.stderr = !verbose)

    # Then get the region and append
    region_cmd <- sprintf(
      'tabix %s %s >> %s',
      vcf_url,
      region_str,
      shQuote(tmp_vcf)
    )

    if (verbose) message("Fetching region ", region_str, "...")
    exit_code <- system(region_cmd, ignore.stdout = !verbose, ignore.stderr = !verbose)

    if (exit_code != 0) {
      stop(
        "tabix command failed (exit code ", exit_code, ").\n",
        "Check that the URL is accessible and tabix is working.",
        call. = FALSE
      )
    }

    # If we need population filtering but only have tabix (no bcftools),
    # warn the user
    if (!is.null(population) && !has_bcftools) {
      warning(
        "bcftools not found: cannot filter by population. ",
        "The VCF will contain all 2504 samples. ",
        "Install bcftools to enable population filtering.",
        call. = FALSE
      )
    }

    # Compress with bgzip if available, otherwise gzip
    if (has_bcftools) {
      # bcftools can convert
      system(sprintf("bcftools view -Oz -o %s %s",
                     shQuote(output_path), shQuote(tmp_vcf)),
             ignore.stdout = TRUE, ignore.stderr = TRUE)
    } else {
      # Try bgzip, fall back to gzip
      bgzip_available <- !is.null(tryCatch(
        system2("bgzip", "--version", stdout = TRUE, stderr = TRUE),
        error = function(e) NULL
      ))

      if (bgzip_available) {
        file.copy(tmp_vcf, sub("\\.gz$", "", output_path))
        system(sprintf("bgzip %s", shQuote(sub("\\.gz$", "", output_path))),
               ignore.stdout = TRUE, ignore.stderr = TRUE)
      } else {
        # Use R's gzip — note this won't be tabix-indexable, but sim1000G
        # doesn't need the index
        con_in  <- file(tmp_vcf, "rb")
        con_out <- gzfile(output_path, "wb")
        while (length(chunk <- readBin(con_in, "raw", n = 1e6)) > 0) {
          writeBin(chunk, con_out)
        }
        close(con_in)
        close(con_out)
      }
    }
  }

  # --- Verify -----------------------------------------------------------------

  if (!file.exists(output_path) || file.size(output_path) < 100) {
    stop(
      "Output VCF file is missing or empty. The download may have failed.\n",
      "Ensure you have internet access and the 1000 Genomes FTP is reachable.",
      call. = FALSE
    )
  }

  if (verbose) {
    size_mb <- round(file.size(output_path) / 1e6, 2)
    message(sprintf("Saved: %s (%.2f MB)", output_path, size_mb))
  }

  output_path
}


# --- Batch download multiple regions ------------------------------------------

#' Download multiple genomic regions from 1000 Genomes Phase 3
#'
#' Convenience wrapper around \code{\link{download_1kg_region}} that
#' downloads several regions in one call. Returns a vector of VCF paths
#' ready to pass to \code{\link{simulate_genotypes}}.
#'
#' @param regions A data.frame with columns \code{chr}, \code{start},
#'   \code{end}. Optionally a \code{label} column for naming. Use
#'   \code{\link{default_benchmark_regions}()} to get a curated set.
#' @param n_regions Integer. How many regions to download from the provided
#'   table. If NULL, download all rows. Default: NULL.
#' @param output_dir Directory to save VCF files. Default: "data/vcf" in the
#'   current working directory.
#' @param population Character. Population to subset to. Default: "EUR".
#' @param force Logical. Re-download even if files exist. Default: FALSE.
#' @param verbose Logical. Default: TRUE.
#'
#' @return Character vector of paths to downloaded VCF files.
#'
#' @examples
#' \dontrun{
#' # Download 5 default regions for European samples
#' vcf_paths <- download_1kg_regions(
#'   regions = default_benchmark_regions(),
#'   n_regions = 5,
#'   output_dir = "data/vcf"
#' )
#'
#' # Use them with simulate_genotypes()
#' geno <- simulate_genotypes(
#'   n_regions = 5, n = 1000, p = 500,
#'   vcf_files = vcf_paths
#' )
#' }
#'
#' @export
download_1kg_regions <- function(regions = default_benchmark_regions(),
                                 n_regions = NULL,
                                 output_dir = file.path("data", "vcf"),
                                 population = "EUR",
                                 force = FALSE,
                                 verbose = TRUE) {

  # Validate regions data.frame
  required_cols <- c("chr", "start", "end")
  missing_cols <- setdiff(required_cols, names(regions))
  if (length(missing_cols) > 0) {
    stop(
      "regions data.frame must have columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Subset if n_regions is specified
  if (!is.null(n_regions)) {
    if (n_regions > nrow(regions)) {
      warning(
        sprintf(
          "Requested %d regions but only %d available. Using all %d.",
          n_regions, nrow(regions), nrow(regions)
        ),
        call. = FALSE
      )
      n_regions <- nrow(regions)
    }
    regions <- regions[seq_len(n_regions), , drop = FALSE]
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) message("Created directory: ", output_dir)
  }

  # Download each region
  vcf_paths <- character(nrow(regions))
  for (i in seq_len(nrow(regions))) {

    label <- if ("label" %in% names(regions)) regions$label[i] else NULL
    filename <- if (!is.null(label)) {
      sprintf("%s.vcf.gz", label)
    } else {
      NULL
    }

    if (verbose) {
      message(sprintf(
        "\n--- Region %d/%d: chr%s:%d-%d ---",
        i, nrow(regions), regions$chr[i], regions$start[i], regions$end[i]
      ))
    }

    vcf_paths[i] <- download_1kg_region(
      chr = regions$chr[i],
      start = regions$start[i],
      end = regions$end[i],
      output_dir = output_dir,
      population = population,
      output_filename = filename,
      force = force,
      verbose = verbose
    )
  }

  vcf_paths
}
