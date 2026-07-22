# =============================================================================
# wrapper_fb_joint.R
#
# CROSS-REGION JOINT (Functional) BEATRICE  -  Iteration 003 model changes.
#
# Two methods that train ONE shared annotation prior head JOINTLY across all
# regions of a scenario (not the sequential warm-start chain, not per-locus):
#
#   fb_pooled   -> prior_head = 'linear'    (user idea #1, flagship):
#                  a shared logistic annotation->prior map, funmap/PolyFun-style.
#   fb_xregion  -> prior_head = 'lassonet'  (user idea #2):
#                  the LassoNet prior, shared and jointly trained across regions.
#
# Objective  L(phi,{psi_r}) = sum_r ELBO_r( psi_r ; p0_r = f_phi(v_r) ), with the
# shared head phi updated ONCE per step from the equal-weighted mean of the
# per-region gradients (EQUAL WEIGHT PER REGION; see docs/autoresearch/iteration-003.md).
#
# Mechanism. run_methods() calls a method's `run_<method>_scenario_setup(genotypes,
# regions, user_args)` ONCE per scenario (all regions visible) and merges the
# returned list into every region's args. We use that hook to run the joint Python
# trainer (BEATRICE_annot_sparse/beatrice_joint.py) over ALL regions in one process,
# parse each region's pip.csv/credible_set.txt, and return a cache keyed by a
# z-fingerprint. The thin per-region function looks up its own region's result.
#
# Falls back to per-region (single-locus) FB when a scenario has no annotations
# (the `none` arm -> plain BEATRICE) or when the joint run fails for a region;
# such results are tagged additional$joint_fallback = TRUE so they are traceable.
# =============================================================================


# --- deterministic per-region key from the z vector (always available) -------
# Computed identically in the scenario hook (regions[[i]]$z) and the per-region
# wrapper (region_pheno$z); those are the same in-memory vector, so the key is
# exact. Does not rely on region_id being set by the simulator.
.fb_fingerprint <- function(z) {
  z <- as.numeric(z)
  paste(length(z),
        format(sum(z),      digits = 15),
        format(sum(z * z),  digits = 15),
        format(z[1],        digits = 15),
        format(z[length(z)], digits = 15),
        sep = "|")
}


# --- write one region's .z/.ld/.annot into work_dir (mirrors FB wrapper I/O) --
.fb_joint_write_region <- function(work_dir, idx, z, LD, annotations, variant_ids) {
  p <- length(z)
  if (is.null(variant_ids)) variant_ids <- paste0("rs", seq(0, p - 1))
  variant_ids <- gsub("\\s+", "_", variant_ids)   # collapse embedded spaces (see FB wrapper)

  z_path     <- file.path(work_dir, sprintf("region%d.z", idx))
  ld_path    <- file.path(work_dir, sprintf("region%d.ld", idx))
  annot_path <- file.path(work_dir, sprintf("region%d.annot", idx))
  target_dir <- file.path(work_dir, sprintf("out%d", idx))
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  writeLines(paste(variant_ids, z, sep = " "), z_path)
  utils::write.table(round(LD, 8), ld_path,
                     quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
  annot_df <- cbind(variant_ids, as.data.frame(annotations))
  utils::write.table(annot_df, annot_path,
                     quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

  list(z = z_path, LD = ld_path, annot = annot_path, target = target_dir,
       variant_ids = variant_ids, p = p)
}


# --- parse one region's joint output (same contract as the FB wrapper) --------
.fb_joint_parse_output <- function(out_dir, variant_ids, p) {
  pip_path  <- file.path(out_dir, "pip.csv")
  cred_path <- file.path(out_dir, "credible_set.txt")
  cond_path <- file.path(out_dir, "conditional_credible_variants_probability.txt")
  if (!file.exists(pip_path)) return(NULL)

  pip_df <- tryCatch(utils::read.csv(pip_path, stringsAsFactors = FALSE),
                     error = function(e) NULL)
  if (is.null(pip_df) || !"pip" %in% names(pip_df)) return(NULL)

  if ("variant_names" %in% names(pip_df)) {
    ord <- match(variant_ids, pip_df$variant_names)
    if (anyNA(ord)) ord <- seq_len(p)
  } else {
    ord <- seq_len(p)
  }
  pip <- pmax(0, pmin(1, as.numeric(pip_df$pip[ord])))

  credible_sets <- list()
  if (file.exists(cred_path)) {
    cred_lines <- readLines(cred_path, warn = FALSE)
    cred_lines <- cred_lines[nchar(trimws(cred_lines)) > 0]
    if (length(cred_lines) > 0) {
      credible_sets <- lapply(cred_lines, function(line) {
        idx <- suppressWarnings(as.integer(strsplit(trimws(line), "\\s+")[[1]]))
        sort(idx[!is.na(idx)] + 1L)                # 0-based -> 1-based
      })
      credible_sets <- credible_sets[lengths(credible_sets) > 0]
    }
  }

  cs_pip <- NULL
  if (file.exists(cond_path)) {
    cond_lines <- readLines(cond_path, warn = FALSE)
    cond_lines <- cond_lines[nchar(trimws(cond_lines)) > 0]
    cs_pip <- lapply(cond_lines, function(line)
      suppressWarnings(as.numeric(strsplit(trimws(line), "\\s+")[[1]])))
  }

  list(pip = pip, credible_sets = credible_sets, cs_pip = cs_pip)
}


# --- the shared scenario-setup core ------------------------------------------
# Returns list(.fb_joint_cache = <fingerprint -> parsed result>,
#              .fb_joint_annotated = TRUE/FALSE, .fb_joint_head = head).
# Empty annotated-arm cache means the joint run failed -> per-region fallback.
.fb_joint_scenario_setup <- function(genotypes, regions, user_args, prior_head) {
  n_regions <- length(regions)

  # collect regions that have an annotation matrix (joint prior needs annotations)
  specs <- vector("list", n_regions)
  keys  <- character(n_regions)
  m_ann <- NULL
  for (i in seq_len(n_regions)) {
    A <- genotypes[[i]]$annotations_matrix
    if (is.null(A)) A <- regions[[i]]$annotations_matrix
    z <- regions[[i]]$z
    if (is.null(A) || !is.matrix(A) || is.null(z) || nrow(A) != length(z)) {
      return(list(.fb_joint_annotated = FALSE, .fb_joint_head = prior_head))
    }
    if (is.null(m_ann)) m_ann <- ncol(A) else if (ncol(A) != m_ann) {
      return(list(.fb_joint_annotated = FALSE, .fb_joint_head = prior_head))
    }
    specs[[i]] <- list(z = z, LD = genotypes[[i]]$LD, A = A,
                       variant_ids = genotypes[[i]]$variant_ids, n = genotypes[[i]]$n)
    keys[i] <- .fb_fingerprint(z)
  }

  beatrice_dir <- user_args$beatrice_dir
  python       <- user_args$python
  script       <- file.path(beatrice_dir, "beatrice_joint.py")
  if (is.null(beatrice_dir) || is.null(python) || !file.exists(script)) {
    return(list(.fb_joint_annotated = TRUE, .fb_joint_head = prior_head,
                .fb_joint_cache = list()))
  }
  py <- if (file.exists(path.expand(python))) normalizePath(path.expand(python)) else python

  work_dir <- tempfile(pattern = "fbjoint_")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  # write regions + manifest
  written <- vector("list", n_regions)
  man <- c("z\tLD\tannot\ttarget\tN")
  for (i in seq_len(n_regions)) {
    w <- .fb_joint_write_region(work_dir, i, specs[[i]]$z, specs[[i]]$LD,
                                specs[[i]]$A, specs[[i]]$variant_ids)
    written[[i]] <- w
    man <- c(man, sprintf("%s\t%s\t%s\t%s\t%d",
                          w$z, w$LD, w$annot, w$target, as.integer(specs[[i]]$n)))
  }
  manifest_path <- file.path(work_dir, "manifest.tsv")
  writeLines(man, manifest_path)

  gv <- function(k, d) { v <- user_args[[k]]; if (is.null(v)) d else v }
  args <- c(
    script,
    "--manifest",             manifest_path,
    "--prior_head",           prior_head,
    "--max_iter",             as.character(as.integer(gv("max_iter", 1500))),
    "--n_caus",               as.character(as.integer(gv("n_caus", 5))),
    "--sigma_sq",             as.character(gv("sigma_sq", 0.05)),
    "--gamma_coverage",       as.character(gv("gamma_coverage", 0.95)),
    "--sparse_concrete",      as.character(as.integer(gv("sparse_concrete", 50))),
    "--prior_regularisation", as.character(gv("prior_regularisation", 1.0)),
    "--lambda_l1",            as.character(gv("lambda_l1", 0.01)),
    "--hierarchy_M",          as.character(gv("hierarchy_M", 10.0))
  )

  run_output <- tryCatch(
    system2(py, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) structure(conditionMessage(e), class = "fb_joint_error"))

  if (inherits(run_output, "fb_joint_error")) {
    message(sprintf("    WARNING: fb_joint (%s) system2 failed: %s",
                    prior_head, as.character(run_output)))
    return(list(.fb_joint_annotated = TRUE, .fb_joint_head = prior_head,
                .fb_joint_cache = list()))
  }

  # parse per-region outputs into the fingerprint-keyed cache
  cache <- list()
  for (i in seq_len(n_regions)) {
    parsed <- .fb_joint_parse_output(written[[i]]$target,
                                     written[[i]]$variant_ids, written[[i]]$p)
    if (!is.null(parsed)) cache[[keys[i]]] <- parsed
  }
  if (length(cache) == 0L) {
    message(sprintf("    WARNING: fb_joint (%s) produced no parseable region output.",
                    prior_head))
  }
  list(.fb_joint_cache = cache, .fb_joint_annotated = TRUE, .fb_joint_head = prior_head)
}


# --- public scenario-setup hooks (looked up by name in run_methods) ----------
#' @export
run_fb_pooled_scenario_setup <- function(genotypes, regions, user_args) {
  .fb_joint_scenario_setup(genotypes, regions, user_args, prior_head = "linear")
}
#' @export
run_fb_xregion_scenario_setup <- function(genotypes, regions, user_args) {
  .fb_joint_scenario_setup(genotypes, regions, user_args, prior_head = "lassonet")
}


# --- the shared per-region core ----------------------------------------------
.fb_joint_region <- function(region_geno, region_pheno, method_name,
                             .fb_joint_cache = NULL, .fb_joint_annotated = TRUE,
                             .fb_joint_head = NULL, ...) {
  p <- length(region_pheno$z)
  key <- .fb_fingerprint(region_pheno$z)
  entry <- if (!is.null(.fb_joint_cache)) .fb_joint_cache[[key]] else NULL

  if (!is.null(entry)) {
    return(list(
      pip             = entry$pip,
      credible_sets   = entry$credible_sets,
      method          = method_name,
      input_type      = "summary",
      params          = list(prior_head = .fb_joint_head, joint = TRUE),
      runtime_seconds = NA_real_,
      additional      = list(cs_pip = entry$cs_pip, joint_fallback = FALSE)
    ))
  }

  # Cache miss: no annotations (none arm) or joint failed -> per-region FB.
  # On the none arm this is plain BEATRICE; on an annotated arm it is single-locus
  # FB. Tagged joint_fallback = TRUE so it is distinguishable in analysis.
  fb <- tryCatch(
    run_functional_beatrice_region(region_geno = region_geno,
                                   region_pheno = region_pheno, ...),
    error = function(e) list(pip = rep(NA_real_, p), credible_sets = list(),
                             method = method_name, input_type = "summary",
                             params = list(), runtime_seconds = NA_real_,
                             additional = list(), error = conditionMessage(e)))
  fb$method <- method_name
  if (is.null(fb$additional)) fb$additional <- list()
  fb$additional$joint_fallback <- TRUE
  fb
}

#' @export
run_fb_pooled_region <- function(region_geno, region_pheno, ...) {
  .fb_joint_region(region_geno, region_pheno, method_name = "fb_pooled", ...)
}
#' @export
run_fb_xregion_region <- function(region_geno, region_pheno, ...) {
  .fb_joint_region(region_geno, region_pheno, method_name = "fb_xregion", ...)
}
