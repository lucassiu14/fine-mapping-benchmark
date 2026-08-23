# =============================================================================
# >>> ITERATION 005 - EXPLORATORY FORK. DELETE WITH THE ITERATION. <<<
#
# An UNMODIFIED copy of iter004_collect.R. Collect needed no changes for
# Iteration 005 - it stores the job label verbatim and parses nothing - but it
# is forked anyway so that every script Iteration 005 runs is named iter005_*.
# The point is to remove any reason to edit an iter004_* file: the standard's
# scripts must stay byte-identical to what produced the standing results.
# =============================================================================
#!/usr/bin/env Rscript
# =============================================================================
# scripts/analysis/iter004_collect.R
#
# STAGE A of the Iteration 004 analysis pipeline: re-collect from the stored PIPs
# into the L1 fit-level table specified in variable-importance-analysis.md §5.5.
#
# WHY RE-COLLECT AT ALL. The benchmark already wrote evaluation.rds per scenario,
# but those numbers cannot support the analysis:
#   - AP was read off the 0.005 FDR grid, which imports a phi-dependent bias into
#     the estimator of a metric whose phi dependence we are measuring (§5.1);
#   - AP was computed on the POOLED ranking of a scenario's two same-size regions,
#     letting one locus's false positives outrank the other's true positive (§5.1b);
#   - the replicate dimension was summed away, so there is no paired SE for the
#     Sense D gate and no within-cell noise floor for the §7 correction;
#   - there is no region identifier, so sigma^2_u is unidentifiable and p_causal
#     and enrichment_fold have NO VALID ERROR TERM AT ALL (§3.4b).
#
# This reads results.rds (raw PIPs) and sim.rds (truth) and recomputes everything
# from scratch. NO REFITTING - the expensive part is already done.
#
# ARRAYED OVER ROWS. One task per params_grid row (45), because sim.rds is cached
# per row and holds the truth for all 250 of that row's scenarios.
#
# Usage:
#   Rscript scripts/analysis/iter004_collect.R <array_idx> <bench_root> <out_dir>
#
# STATUS: written, never executed.
# =============================================================================

suppressWarnings(suppressMessages({
  if (requireNamespace("fmbenchmark", quietly = TRUE)) library(fmbenchmark)
}))
source(file.path("R", "evaluate_extras.R"))

args       <- commandArgs(trailingOnly = TRUE)
array_idx  <- as.integer(args[1])
bench_root <- args[2]
out_dir    <- args[3]
stopifnot(!is.na(array_idx), dir.exists(bench_root))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Decision thresholds (§4.6 P1) and calibration bands (§4.3 Q2). Counts are stored,
# never per-fit RATES: a rate on 0-2 selections is unestimable, and rates must be
# formed from counts POOLED to at least region level (§3.2).
THRESH <- c(50, 80, 90, 95, 99) / 100
BANDS  <- list(lo = c(0, 0.1), mid = c(0.1, 0.5), hi = c(0.5, 0.9),
               top = c(0.9, 1 + 1e-9))

job_dirs <- sort(list.dirs(bench_root, recursive = FALSE))
job_dirs <- job_dirs[grepl("/job_[0-9]+", job_dirs)]
if (array_idx > length(job_dirs)) {
  message("array index ", array_idx, " > ", length(job_dirs), " rows; nothing to do")
  quit(save = "no", status = 0)
}
job_dir <- job_dirs[array_idx]
label   <- basename(job_dir)
message("row ", array_idx, ": ", label)

sim_file <- file.path(job_dir, "sim.rds")
if (!file.exists(sim_file)) {
  stop("no sim.rds in ", job_dir, " - the truth vectors are unavailable and AP ",
       "cannot be recomputed from PIPs.", call. = FALSE)
}
sim <- readRDS(sim_file)
p_nominal <- sim$params$p            # nominal region size, indexed by region_id

#' region_idx: which of the 2 same-size regions this is (§5.5, mandatory).
#'
#' P_VECTOR is 2 regions per size class in order, so region_id 1,2 are the first
#' size class, 3,4 the second, and so on. Without this the two regions in a class
#' cannot be told apart and sigma^2_u is unidentifiable.
region_idx_of <- function(region_id) ((region_id - 1L) %% 2L) + 1L

#' Everything computed from one (pip, truth) pair. This is the ONLY place raw
#' PIPs are touched; every downstream statistic is a function of these columns.
fit_row <- function(pip, causal_idx, n_var) {
  y <- logical(n_var); y[causal_idx] <- TRUE
  S <- sum(y)

  rs <- .rank_set_metrics(pip, y, alpha = c(0.5, 0.95), S = S)

  # Counts at the decision thresholds. Cumulative in the usual sense: nsel(t) is
  # the number of variants with pip >= t.
  cnt <- lapply(THRESH, function(t) {
    sel  <- pip >= t
    nsel <- sum(sel)
    c(tp = sum(y & sel), nsel = nsel, sum_pip_sel = sum(pip[sel]))
  })
  names(cnt) <- sprintf("%02d", THRESH * 100)

  # Calibration band counts. Bands, not the full 26-bin table: the bin table is
  # rebuilt at L3 where there are enough variants for it to mean anything.
  bnd <- lapply(BANDS, function(b) {
    sel <- pip >= b[1] & pip < b[2]
    c(n = sum(sel), c = sum(y & sel), sum_pip = sum(pip[sel]))
  })

  out <- list(
    ap = .compute_ap_exact(pip, y),
    n_variants = n_var, n_causal = S, sum_pip_total = sum(pip)
  )
  out <- c(out, rs)
  for (nm in names(cnt)) {
    out[[paste0("tp_at_",          nm)]] <- unname(cnt[[nm]]["tp"])
    out[[paste0("nsel_at_",        nm)]] <- unname(cnt[[nm]]["nsel"])
    out[[paste0("sum_pip_sel_at_", nm)]] <- unname(cnt[[nm]]["sum_pip_sel"])
  }
  for (nm in names(bnd)) {
    out[[paste0("n_band_",       nm)]] <- unname(bnd[[nm]]["n"])
    out[[paste0("c_band_",       nm)]] <- unname(bnd[[nm]]["c"])
    out[[paste0("sum_pip_band_", nm)]] <- unname(bnd[[nm]]["sum_pip"])
  }
  out
}

scen_dirs <- sort(list.dirs(job_dir, recursive = FALSE))
scen_dirs <- scen_dirs[grepl("/scenario_[0-9]+", scen_dirs)]
message("  scenarios on disk: ", length(scen_dirs))

rows <- vector("list", length(scen_dirs) * 200L)
k <- 0L
n_skipped <- 0L
n_badpip  <- 0L   # fits carrying non-finite PIPs (see the guard below)
n_errored <- 0L   # fits whose metric computation threw despite the guard

for (sd in scen_dirs) {
  # THE SCENARIO INDEX COMES FROM THE DIRECTORY NAME, NOT FROM THE FIT.
  # run_benchmark_job.R subsets sim to a single scenario and resets
  # `mini$scenarios[[1]]$scenario_id <- 1L`, because evaluate_methods() uses
  # scenario_id as a LIST INDEX into the subset. So every saved fit carries
  # scenario_id = 1 and the true index survives only in scenario_<sc>/.
  # Indexing sim$scenarios by fit$scenario_id would silently give every fit
  # scenario 1's truth - wrong S, wrong phi, wrong causal indices.
  sc_idx <- suppressWarnings(as.integer(sub("^scenario_", "", basename(sd))))
  if (is.na(sc_idx) || sc_idx < 1L || sc_idx > length(sim$scenarios)) {
    n_skipped <- n_skipped + 1L; next
  }
  sc <- sim$scenarios[[sc_idx]]

  rf <- file.path(sd, "results.rds")
  if (!file.exists(rf)) { n_skipped <- n_skipped + 1L; next }
  res <- tryCatch(readRDS(rf), error = function(e) NULL)
  if (is.null(res)) { n_skipped <- n_skipped + 1L; next }

  for (m in names(res)) {
    v <- res[[m]]
    # results.rds carries non-method top-level entries (methods_run,
    # run_timestamp, simulation_params); skip anything without $results.
    if (!is.list(v) || is.null(v$results)) next

    for (fit in v$results) {
      if (is.null(fit$pip) || is.null(fit$scenario_id) || is.null(fit$region_id)) next
      failed <- !is.null(fit$error)

      tr <- sc$regions[[fit$region_id]]$truth
      if (is.null(tr)) next

      n_var <- length(fit$pip)
      # is.finite() rather than is.na(): exp() overflow yields Inf, and Inf/sum
      # yields NaN, so both must be caught. A zero-length PIP vector is equally
      # unusable and is folded into the same flag.
      # Scoped to fits that did NOT report an error, so bad_pip means exactly
      # "the method claimed success and returned non-finite output" - the
      # alarming case. A fit that already declared failure is accounted for by
      # `failed`; letting it also set bad_pip would bury the overflow signal
      # under structural absences like funmap on the unannotated arm.
      bad_pip <- !failed && (n_var == 0L || !all(is.finite(fit$pip)))
      if (bad_pip) n_badpip <- n_badpip + 1L
      rsz <- if (!is.null(p_nominal) && length(p_nominal) >= fit$region_id)
               as.integer(p_nominal[fit$region_id]) else n_var

      base <- list(
        job_dir     = label,
        scenario_id = sc_idx,
        S           = tr$S,
        phi         = tr$phi,
        iter        = sc$iter %||% NA_integer_,
        region_id   = fit$region_id,
        region_idx  = region_idx_of(fit$region_id),
        region_size = rsz,
        method      = m,
        failed      = failed,
        bad_pip     = bad_pip
      )

      # A failed fit contributes an all-NA metric row, NOT a dropped row. Its
      # absence must stay visible downstream - funmap on the `none` arm is
      # structural, not missing at random, and an NA-dropping aggregation would
      # quietly change the method set per cell.
      #
      # NON-FINITE PIPs COUNT AS FAILURE (Iteration 004). The phi=0.6 arm added
      # this iteration produced fits carrying NaN/Inf PIPs - numerical overflow
      # at high heritability - which is what killed evaluate_methods() on 389
      # scenarios with "'vec' must be sorted non-decreasingly and not contain
      # NAs". The PARTIALLY non-finite vector is the dangerous case: the old
      # all(is.na()) test returns FALSE for it, .compute_ap_exact then returns
      # NA via a length-recycled comparison (with only a warning), and
      # .rank_set_metrics throws outright.
      #
      # Such a posterior has undefined total mass and an undefined ranking for
      # the affected variants, so NO metric derived from it is trustworthy: the
      # fit is marked failed, not repaired. Imputing the non-finite entries to 0
      # would invent posterior mass and would flatter precisely the method that
      # overflowed. bad_pip is carried as its own column so this loss is
      # countable per method x phi rather than being conflated with `failed`.
      met <- if (failed || bad_pip) {
        proto <- fit_row(rep(0, max(n_var, 1L)), integer(0), max(n_var, 1L))
        lapply(proto, function(z) NA_real_)
      } else {
        # Belt and braces: one pathological fit must not kill a 250-scenario task.
        tryCatch(fit_row(fit$pip, tr$causal_indices, n_var), error = function(e) {
          n_errored <<- n_errored + 1L
          proto <- fit_row(rep(0, max(n_var, 1L)), integer(0), max(n_var, 1L))
          lapply(proto, function(z) NA_real_)
        })
      }

      k <- k + 1L
      rows[[k]] <- c(base, met)
    }
  }
}

if (k == 0L) {
  message("  no fits collected for this row")
  quit(save = "no", status = 0)
}
rows <- rows[seq_len(k)]
L1 <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))

out_file <- file.path(out_dir, sprintf("L1_%s.rds", label))
saveRDS(L1, out_file)
message(sprintf("  wrote %s  (%d fit rows, %d scenarios skipped)",
                basename(out_file), nrow(L1), n_skipped))
# Surfaced per row so a systematic overflow shows up in the Stage A logs rather
# than only as thinned cells much later in the decomposition.
if (n_badpip > 0L || n_errored > 0L) {
  message(sprintf("  NON-FINITE PIPs: %d fits marked failed; %d metric errors caught",
                  n_badpip, n_errored))
  bp <- L1[L1$bad_pip, , drop = FALSE]
  if (nrow(bp)) {
    tb <- table(bp$method, bp$phi)
    message("  bad_pip fits by method x phi:")
    for (i in seq_len(nrow(tb))) {
      message(sprintf("    %-22s %s", rownames(tb)[i],
                      paste(sprintf("phi=%s:%d", colnames(tb), tb[i, ]), collapse = "  ")))
    }
  }
}
