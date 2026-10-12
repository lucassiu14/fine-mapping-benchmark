#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/build_main_simulation_deposit.R
#
# Stage the main simulation's results for deposit, at ONE grain: a fit, meaning
# one method applied to one region of one iteration. Per-cell means, false
# discovery rates and calibration curves are all aggregates of these files,
# computed by rules the data dictionary states, and are deliberately not
# deposited: holding several grains in one folder is what makes a results
# directory hard to read.
#
# Identified by simulation FACTORS, not by how the study was run. job_dir and
# scenario_id were a job-scheduler directory name and an index into one row's
# 250 scenarios; neither is a property of the simulation. They are dropped, and
# a fit is keyed by
#     model, p_causal, annotation_type, enrichment_fold   (the parameter row)
#     S, phi, iter, region_id, method                     (within it)
# which is unique over all 2,031,230 fits. The four row-level factors name the
# file rather than repeating down every row of it.
#
# Also dropped: `pc` and `enrich`, character duplicates of p_causal and
# enrichment_fold carrying "sparse"/"none" where those hold NA; `region_idx`,
# derivable from region_id; `n_ref`, empty because this study used in-sample LD
# throughout; and `n` and `n_iter`, constants recorded in simulation-grid.csv.
#
#   Rscript scripts/analysis/build_main_simulation_deposit.R [SRC] [DEST]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
SRC  <- if (length(args) >= 1) args[1] else "results/iter004"
DEST <- if (length(args) >= 2) args[2] else "zenodo/main-simulation"
stopifnot(dir.exists(SRC))

unlink(DEST, recursive = TRUE)
for (d in c("", "posterior-inclusion-probabilities", "credible-sets"))
  dir.create(file.path(DEST, d), recursive = TRUE, showWarnings = FALSE)

fld <- function(recs, f, empty) vapply(recs, function(r) {
  v <- r[[f]]; if (is.null(v) || !length(v)) empty else v[1] }, empty)

# The parameter row's filename, built from its factors.
stem <- function(model, p_causal, annotation_type, enrichment_fold) {
  m <- if (identical(model, "sparse_inf")) sprintf("sparse-inf_rho%s", p_causal) else "sparse"
  a <- if (identical(annotation_type, "none")) "none"
       else sprintf("%s_fold%s", annotation_type, enrichment_fold)
  paste(m, a, sep = "_")
}

ROW_KEY <- c("model", "p_causal", "annotation_type", "enrichment_fold")
FIT_KEY <- c("S", "phi", "iter", "region_id", "method")

# --- the map from how it ran to what it is ----------------------------------
L1 <- readRDS(file.path(SRC, "combined_fit_metrics.rds"))
rows <- unique(L1[, c("job_dir", ROW_KEY)])
rows$file <- vapply(seq_len(nrow(rows)), function(i)
  stem(rows$model[i], rows$p_causal[i], rows$annotation_type[i], rows$enrichment_fold[i]), "")
stopifnot(nrow(rows) == 45L, !anyDuplicated(rows$file))
# scenario_id is an index into a row's 250 scenarios; (S, phi, iter) is what it means.
scen <- unique(L1[, c("job_dir", "scenario_id", "S", "phi", "iter")])
cat("45 parameter rows, named by factor; scenario_id resolves to (S, phi, iter)\n\n")

per_fit <- vector("list", 45L); scen_rt <- vector("list", 45L); fit_rt <- vector("list", 45L)

for (i in seq_len(45L)) {
  jd <- rows$job_dir[i]; out <- rows$file[i]
  sc <- scen[scen$job_dir == jd, ]
  resolve <- function(sid) sc[match(sid, sc$scenario_id), c("S", "phi", "iter")]

  # --- posterior inclusion probabilities -> one row per retained variant -----
  p  <- readRDS(file.path(SRC, "piptail", paste0("piptail_", jd, ".rds")))
  tl <- p$tail; n <- vapply(tl, function(r) length(r$pip), 0L)
  sid <- rep(fld(tl, "scenario_id", 0L), n)
  pip <- cbind(resolve(sid), data.frame(
    region_id = rep(fld(tl, "region_id", 0L), n),
    method    = rep(fld(tl, "method", ""),    n),
    variant   = unlist(lapply(tl, `[[`, "idx"), use.names = FALSE),
    pip       = unlist(lapply(tl, `[[`, "pip"), use.names = FALSE),
    is_causal = as.integer(unlist(lapply(tl, `[[`, "is_causal"), use.names = FALSE)),
    stringsAsFactors = FALSE))
  rownames(pip) <- NULL
  stopifnot(nrow(pip) == sum(n), !anyNA(pip$S))
  saveRDS(pip, file.path(DEST, "posterior-inclusion-probabilities", paste0(out, ".rds")))

  pf <- cbind(resolve(fld(tl, "scenario_id", 0L)), data.frame(
    region_id = fld(tl, "region_id", 0L), method = fld(tl, "method", ""),
    pip_n_retained       = n,
    pip_n_below_floor    = fld(tl, "n_below", 0L),
    pip_mass_below_floor = fld(tl, "sum_pip_below", 0), stringsAsFactors = FALSE))
  rm(p, tl, pip)

  # --- credible sets -> one row per variant per reported set -----------------
  cs <- readRDS(file.path(SRC, "credsets", paste0("cs_", jd, ".rds")))
  rw <- cs$rows; nset <- fld(rw, "n_sets", 0L); keep <- which(nset > 0L)
  # Vectorised: building one data frame per fit and rbinding 42,000 of them is
  # quadratic and takes hours over the 45 rows.
  kr <- rw[keep]
  sz <- lapply(kr, function(r) vapply(r$sets, length, 0L))
  nv <- vapply(sz, sum, 0L)
  sets <- if (length(kr)) cbind(
    resolve(rep(fld(kr, "scenario_id", 0L), nv)),
    data.frame(
      region_id = rep(fld(kr, "region_id", 0L), nv),
      method    = rep(fld(kr, "method", ""),    nv),
      set = unlist(lapply(kr, function(r) {
        s <- r$sets
        l <- if (is.null(names(s))) as.character(seq_along(s)) else names(s)
        rep(l, vapply(s, length, 0L))
      }), use.names = FALSE),
      variant = unlist(lapply(kr, function(r) unlist(r$sets, use.names = FALSE)),
                       use.names = FALSE),
      stringsAsFactors = FALSE)) else NULL
  if (is.null(sets))
    sets <- data.frame(S = integer(0), phi = numeric(0), iter = integer(0),
                       region_id = integer(0), method = character(0),
                       set = character(0), variant = integer(0))
  rownames(sets) <- NULL
  saveRDS(sets, file.path(DEST, "credible-sets", paste0(out, ".rds")))

  csf <- cbind(resolve(fld(rw, "scenario_id", 0L)), data.frame(
    region_id = fld(rw, "region_id", 0L), method = fld(rw, "method", ""),
    native_sets_reported = fld(rw, "reported", NA),
    native_n_sets = nset, stringsAsFactors = FALSE))
  rm(cs, rw, sets)

  # --- runtimes --------------------------------------------------------------
  a <- readRDS(file.path(SRC, "aux", paste0("aux_", jd, ".rds")))
  rt <- a$runtime; scope <- fld(rt, "scope", "")
  fr <- rt[scope == "fit"]
  fit_rt[[i]] <- cbind(resolve(fld(fr, "scenario_id", 0L)), data.frame(
    region_id = fld(fr, "region_id", 0L), method = fld(fr, "method", ""),
    runtime_seconds = fld(fr, "runtime_seconds", NA_real_), stringsAsFactors = FALSE))
  st <- rt[scope == "scenario_total"]
  scen_rt[[i]] <- cbind(rows[i, ROW_KEY], resolve(fld(st, "scenario_id", 0L)),
    data.frame(method = fld(st, "method", ""),
               runtime_seconds = fld(st, "runtime_seconds", NA_real_),
               n_total = fld(st, "n_total", NA_integer_),
               n_failed = fld(st, "n_failed", NA_integer_), stringsAsFactors = FALSE))
  rm(a, rt, fr, st)

  per_fit[[i]] <- cbind(job_dir = jd, merge(pf, csf, by = FIT_KEY, all = TRUE))
  rm(pf, csf)
  cat(sprintf("  [%2d/45] %-42s -> %s\n", i, jd, out))
}

# --- fit-metrics -------------------------------------------------------------
DROP <- c("job_dir", "scenario_id", "region_idx", "pc", "enrich", "n_ref", "n", "n_iter")
add <- do.call(rbind, per_fit)
rt  <- do.call(rbind, Map(function(d, j) cbind(job_dir = j, d), fit_rt, rows$job_dir))
rm(per_fit, fit_rt); invisible(gc())

before <- nrow(L1)
K <- c("job_dir", FIT_KEY)
L1 <- merge(L1, rt,  by = K, all.x = TRUE)
L1 <- merge(L1, add, by = K, all.x = TRUE)
if (nrow(L1) != before) stop("merge changed the row count: a key is not unique", call. = FALSE)

L1 <- L1[, setdiff(names(L1), DROP)]
front <- c(ROW_KEY, FIT_KEY, "region_size")
L1 <- L1[, c(front, setdiff(names(L1), front))]
L1 <- L1[do.call(order, unname(as.list(L1[c(ROW_KEY, FIT_KEY)]))), ]
rownames(L1) <- NULL
key <- do.call(paste, c(L1[c(ROW_KEY, FIT_KEY)], sep = "\r"))
stopifnot(!anyDuplicated(key))

cat(sprintf("\nfit-metrics: %s rows, %d columns (runtime on %.1f%%, PIP counts on %.1f%%)\n",
            format(nrow(L1), big.mark = ","), ncol(L1),
            100 * mean(!is.na(L1$runtime_seconds)), 100 * mean(!is.na(L1$pip_n_retained))))
saveRDS(L1, file.path(DEST, "fit-metrics.rds"))
gz <- gzfile(file.path(DEST, "fit-metrics.csv.gz"), "w")
write.csv(L1, gz, row.names = FALSE); close(gz)

sr <- do.call(rbind, scen_rt); rownames(sr) <- NULL
saveRDS(sr, file.path(DEST, "scenario-runtimes.rds"))
cat(sprintf("scenario-runtimes: %s rows\n", format(nrow(sr), big.mark = ",")))

# --- the design, with the filename each row maps to -------------------------
g <- read.csv("scripts/hpc/params_grid.csv", stringsAsFactors = FALSE)
g$job_dir <- sprintf("job_%03d_%s", g$job_id, g$label)
g <- merge(g, rows[, c("job_dir", "file")], by = "job_dir")
g <- g[, c("file", ROW_KEY, "annotation_correlation", "n_annotations", "n_regions",
           "n", "n_iter", "S_values", "phi_values", "p_values", "enrichment_values")]
g <- g[order(g$file), ]
write.csv(g, file.path(DEST, "simulation-grid.csv"), row.names = FALSE)
cat("simulation-grid.csv:", nrow(g), "rows\n")

for (d in c("README.md", "DATA-DICTIONARY.md")) {
  s <- file.path("docs/deposit", d)
  if (!file.exists(s)) stop("missing ", s, call. = FALSE)
  file.copy(s, file.path(DEST, d), overwrite = TRUE)
}
cat("\nstaged to ", DEST, "\n", sep = "")

# --- package exactly what gets uploaded -------------------------------------
# Zenodo stores a flat list of files, not folders, so the structure survives only
# inside the archive. The four small files sit beside it so the record page is
# legible without a 1.4 GB download; each is a copy of a file already inside.
UP <- file.path(dirname(DEST), "upload")
unlink(UP, recursive = TRUE); dir.create(UP, recursive = TRUE)
tar <- file.path(UP, paste0(basename(DEST), ".tar.gz"))
system2("tar", c("-czf", shQuote(tar), "-C", shQuote(dirname(DEST)), shQuote(basename(DEST))))
for (f in c("README.md", "DATA-DICTIONARY.md", "simulation-grid.csv", "fit-metrics.csv.gz"))
  file.copy(file.path(DEST, f), file.path(UP, f), overwrite = TRUE)

cat("\nupload set (", UP, "):\n", sep = "")
inf <- file.info(list.files(UP, full.names = TRUE))
for (i in order(-inf$size))
  cat(sprintf("  %-28s %8.1f MB\n", basename(rownames(inf))[i], inf$size[i] / 1024^2))
cat("\nsha256 of the archive:\n  ",
    system2("shasum", c("-a", "256", shQuote(tar)), stdout = TRUE), "\n")
