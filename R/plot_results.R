# =============================================================================
# plot_results.R
#
# Visualisation module for the fine-mapping benchmark evaluation output.
#
# Usage:
#   plot_results(eval_out, output_file = "results/my_run/results.pdf")
#
# Output PDF sections:
#   1. Global  — PR curve (all methods), PIP calibration (faceted by method),
#                summary table with ± SE
#   2. By S    — PR grid (method × S), calibration grid (method × S),
#                metric line plots vs S with error bars
#   3. By phi  — same structure as By S
#   4. By p_causal — same structure (sparse_inf model only)
# =============================================================================


# =============================================================================
# Color palette and theme
# =============================================================================

# Colorblind-friendly palette (Wong 2011, extended)
.FM_COLORS <- c(
  susie     = "#0072B2",
  susie_inf = "#56B4E9",
  abf       = "#009E73",
  carma     = "#D55E00",
  funmap    = "#E69F00",
  paintor   = "#CC79A7",
  beatrice              = "#F0E442",
  functional_beatrice   = "#E76BF3",
  finemap               = "#999999",
  # marginal_z is the baseline (neutral grey-blue);
  # polyfun_* methods get a related-but-distinct hue family to signal
  # they're an annotation-aware family.
  marginal_z            = "#7F7F7F",   # mid grey — baseline / floor
  polyfun_oracle        = "#1A535C",   # deep teal — the ceiling
  polyfun_est           = "#4ECDC4",   # lighter teal — the (naive) realistic version
  polyfun_ldsc          = "#2E86AB",   # deep blue — the LD-score-corrected version
  sbayesrc              = "#8B5A2B",   # bronze — mixture-of-normals Bayesian
  # Deep wine — visually distinct from the BEATRICE
  # magenta and PAINTOR pink so the variational methods don't blur together.
  sparsepro             = "#6B2D5C"
)

.method_color_scale <- function(methods) {
  known   <- intersect(methods, names(.FM_COLORS))
  unknown <- setdiff(methods, names(.FM_COLORS))
  cols    <- .FM_COLORS[known]
  if (length(unknown) > 0L) {
    extras        <- grDevices::hcl.colors(length(unknown), palette = "Dark 3")
    extras        <- setNames(extras, unknown)
    cols          <- c(cols, extras)
  }
  cols[methods]
}

.fm_theme <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = base_size + 1),
      legend.position = "bottom",
      strip.text      = ggplot2::element_text(size = base_size - 1, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

.METRIC_LABELS <- c(
  auprc          = "AUPRC",
  cs_coverage    = "CS Coverage",
  cs_power       = "CS Power",
  cs_size_median = "Median CS Size",
  runtime_mean   = "Runtime (s)"
)


# =============================================================================
# Data-extraction helpers
# =============================================================================

# Long data frame of PR curve rows; stratum_val = NA for global
.extract_pr <- function(eval_out, methods, stratum = "global") {
  dfs <- lapply(methods, function(m) {
    if (stratum == "global") {
      curve <- eval_out[[m]]$global$fdr_power_curve
      if (is.null(curve)) return(NULL)
      curve$method      <- m
      curve$stratum_val <- NA_real_
      curve
    } else {
      strat <- eval_out[[m]][[stratum]]
      if (is.null(strat)) return(NULL)
      sub <- lapply(names(strat), function(v) {
        curve <- strat[[v]]$fdr_power_curve
        if (is.null(curve)) return(NULL)
        curve$method      <- m
        curve$stratum_val <- as.numeric(v)
        curve
      })
      do.call(rbind, Filter(Negate(is.null), sub))
    }
  })
  do.call(rbind, Filter(Negate(is.null), dfs))
}

# Long data frame of calibration rows
.extract_cal <- function(eval_out, methods, stratum = "global") {
  dfs <- lapply(methods, function(m) {
    if (stratum == "global") {
      cal <- eval_out[[m]]$global$pip_calibration
      if (is.null(cal)) return(NULL)
      cal$method      <- m
      cal$stratum_val <- NA_real_
      cal
    } else {
      strat <- eval_out[[m]][[stratum]]
      if (is.null(strat)) return(NULL)
      sub <- lapply(names(strat), function(v) {
        cal <- strat[[v]]$pip_calibration
        if (is.null(cal)) return(NULL)
        cal$method      <- m
        cal$stratum_val <- as.numeric(v)
        cal
      })
      do.call(rbind, Filter(Negate(is.null), sub))
    }
  })
  do.call(rbind, Filter(Negate(is.null), dfs))
}

# Long data frame of scalar metrics (value + se)
.extract_scalars <- function(eval_out, methods, stratum,
                              metrics = names(.METRIC_LABELS)) {
  rows <- list()
  for (m in methods) {
    if (stratum == "global") {
      source_list <- list("global" = eval_out[[m]]$global)
    } else {
      source_list <- eval_out[[m]][[stratum]]
    }
    if (is.null(source_list)) next
    for (key in names(source_list)) {
      g <- source_list[[key]]
      for (metric in metrics) {
        val <- g[[metric]]
        se  <- g[[paste0(metric, "_se")]]
        rows[[length(rows) + 1L]] <- data.frame(
          method      = m,
          stratum_val = if (stratum == "global") NA_real_ else as.numeric(key),
          metric      = metric,
          value       = if (is.null(val) || length(val) == 0L) NA_real_ else as.numeric(val[[1L]]),
          se          = if (is.null(se)  || length(se)  == 0L) NA_real_ else as.numeric(se[[1L]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}


# =============================================================================
# Individual plot functions
# =============================================================================

# -- Precision-recall curve (global: all methods on one panel) ----------------

.plot_pr_global <- function(df, title = "Global Precision-Recall Curves") {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  df <- df[!is.na(df$recall) & !is.na(df$precision), ]
  if (nrow(df) == 0L) return(NULL)

  cols    <- .method_color_scale(unique(df$method))
  has_se  <- "precision_se" %in% names(df) && any(!is.na(df$precision_se))

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = recall, y = precision, color = method, fill = method, group = method
  ))

  p +
    ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = cols, name = "Method") +
    ggplot2::scale_fill_manual(values = cols,  name = "Method") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(title = title, x = "Recall (Power)", y = "Precision (1 - FDR)") +
    .fm_theme()
}


# -- PR curves faceted by stratum (rows = methods, cols = stratum values) -----

.plot_pr_facet <- function(df, stratum_label, title) {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  df <- df[!is.na(df$recall) & !is.na(df$precision), ]
  if (nrow(df) == 0L) return(NULL)

  cols   <- .method_color_scale(unique(df$method))

  sv_sorted <- sort(unique(df$stratum_val))
  df$col_label <- factor(
    paste0(stratum_label, " = ", df$stratum_val),
    levels = paste0(stratum_label, " = ", sv_sorted)
  )

  ggplot2::ggplot(df, ggplot2::aes(
    x = recall, y = precision, color = method, group = method
  )) +
    ggplot2::geom_line(linewidth = 0.6, na.rm = TRUE) +
    ggplot2::facet_grid(method ~ col_label) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(title = title, x = "Recall", y = "Precision") +
    .fm_theme(base_size = 9) +
    ggplot2::theme(legend.position = "none")
}


# -- PIP calibration (global: faceted by method) ------------------------------

.plot_cal_global <- function(df, title = "Global PIP Calibration") {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  df <- df[!is.na(df$mean_pip), ]
  if (nrow(df) == 0L) return(NULL)

  cols    <- .method_color_scale(unique(df$method))
  has_se  <- "frac_causal_se" %in% names(df) && any(!is.na(df$frac_causal_se))
  n_meth  <- length(unique(df$method))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_pip, y = frac_causal, color = method))

  if (has_se) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(0, frac_causal - frac_causal_se),
        ymax = pmin(1, frac_causal + frac_causal_se)
      ),
      width = 0.025, linewidth = 0.5, alpha = 0.7, na.rm = TRUE
    )
  }

  p +
    ggplot2::geom_point(size = 2.5, na.rm = TRUE) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey50", linewidth = 0.5) +
    ggplot2::facet_wrap(~ method, ncol = min(3L, n_meth)) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(title = title,
                  x = "Mean PIP (expected)", y = "Fraction causal (observed)") +
    .fm_theme() +
    ggplot2::theme(legend.position = "none")
}


# -- PIP calibration faceted by stratum (rows = methods, cols = values) -------

.plot_cal_facet <- function(df, stratum_label, title) {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  df <- df[!is.na(df$mean_pip), ]
  if (nrow(df) == 0L) return(NULL)

  cols   <- .method_color_scale(unique(df$method))
  has_se <- "frac_causal_se" %in% names(df) && any(!is.na(df$frac_causal_se))

  sv_sorted <- sort(unique(df$stratum_val))
  df$col_label <- factor(
    paste0(stratum_label, " = ", df$stratum_val),
    levels = paste0(stratum_label, " = ", sv_sorted)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_pip, y = frac_causal, color = method))

  if (has_se) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(0, frac_causal - frac_causal_se),
        ymax = pmin(1, frac_causal + frac_causal_se)
      ),
      width = 0.025, linewidth = 0.4, alpha = 0.7, na.rm = TRUE
    )
  }

  p +
    ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey50", linewidth = 0.4) +
    ggplot2::facet_grid(method ~ col_label) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    ggplot2::labs(title = title,
                  x = "Mean PIP (expected)", y = "Fraction causal (observed)") +
    .fm_theme(base_size = 9) +
    ggplot2::theme(legend.position = "none")
}


# -- Metrics line plots (metric vs stratum variable, one line per method) -----

.plot_metrics_line <- function(df, x_label, title,
                                metrics = names(.METRIC_LABELS)) {
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  df <- df[!is.na(df$stratum_val) & df$metric %in% metrics, ]
  if (nrow(df) == 0L) return(NULL)

  cols   <- .method_color_scale(unique(df$method))
  has_se <- any(!is.na(df$se))

  # Nice facet labels
  df$metric_label <- factor(
    ifelse(df$metric %in% names(.METRIC_LABELS),
           .METRIC_LABELS[df$metric], df$metric),
    levels = .METRIC_LABELS[metrics[metrics %in% df$metric]]
  )

  # Error bar width: 3% of the x range
  x_range <- diff(range(df$stratum_val, na.rm = TRUE))
  eb_width <- max(x_range * 0.03, 0.02)

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = stratum_val, y = value, color = method, group = method
  ))

  if (has_se) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = value - se, ymax = value + se),
      width = eb_width, linewidth = 0.5, alpha = 0.7, na.rm = TRUE
    )
  }

  p +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.5, na.rm = TRUE) +
    ggplot2::facet_wrap(~ metric_label, scales = "free_y", ncol = 3L) +
    ggplot2::scale_color_manual(values = cols, name = "Method") +
    ggplot2::labs(title = title, x = x_label, y = NULL) +
    .fm_theme()
}


# -- Summary table (global metrics, all methods, formatted with ± SE) ---------

.plot_summary_table <- function(eval_out, methods) {
  fmt <- function(val, se, digits = 3L) {
    if (is.null(val) || length(val) == 0L || is.na(val)) return("NA")
    s <- formatC(as.numeric(val), digits = digits, format = "f")
    if (!is.null(se) && length(se) > 0L && !is.na(se))
      s <- paste0(s, " \u00b1 ", formatC(as.numeric(se), digits = digits, format = "f"))
    s
  }

  rows <- lapply(methods, function(m) {
    g <- eval_out[[m]]$global
    data.frame(
      Method          = m,
      AUPRC           = fmt(g$auprc,          g$auprc_se),
      `CS Coverage`   = fmt(g$cs_coverage,    g$cs_coverage_se),
      `CS Power`      = fmt(g$cs_power,       g$cs_power_se),
      `Med CS Size`   = fmt(g$cs_size_median, g$cs_size_median_se, digits = 1L),
      `Runtime (s)`   = fmt(g$runtime_mean,   g$runtime_mean_se,   digits = 2L),
      `n_fits`        = as.character(g$n_fits),
      `n_failed`      = as.character(g$n_failed),
      check.names     = FALSE,
      stringsAsFactors = FALSE
    )
  })

  tbl <- do.call(rbind, rows)

  tt <- gridExtra::tableGrob(
    tbl, rows = NULL,
    theme = gridExtra::ttheme_default(
      core    = list(fg_params = list(fontsize = 8.5, hjust = 0.5),
                     bg_params = list(fill = c("white", "#F5F5F5"))),
      colhead = list(fg_params = list(fontsize = 9.5, fontface = "bold", hjust = 0.5),
                     bg_params = list(fill = "#DDEEFF"))
    )
  )

  title_g <- grid::textGrob(
    "Global Evaluation Summary (mean \u00b1 SE across replicates)",
    gp = grid::gpar(fontsize = 13, fontface = "bold")
  )

  gridExtra::arrangeGrob(title_g, tt, ncol = 1L,
                          heights = grid::unit(c(0.06, 0.94), "npc"))
}


# =============================================================================
# Section helpers (print one PDF section)
# =============================================================================

.print_section <- function(eval_out, valid_methods, all_methods,
                            stratum, stratum_label, x_label, section_title,
                            metrics = names(.METRIC_LABELS)) {

  # Skip if no data exists for this stratum
  has_data <- any(sapply(valid_methods, function(m) {
    !is.null(eval_out[[m]][[stratum]])
  }))
  if (!has_data) return(invisible(NULL))

  # PR curve grid (methods × stratum values)
  pr_df <- .extract_pr(eval_out, valid_methods, stratum)
  p_pr  <- .plot_pr_facet(pr_df, stratum_label,
                           paste0(section_title, ": Precision-Recall Curves"))
  if (!is.null(p_pr)) print(p_pr)

  # Calibration grid
  cal_df <- .extract_cal(eval_out, valid_methods, stratum)
  p_cal  <- .plot_cal_facet(cal_df, stratum_label,
                             paste0(section_title, ": PIP Calibration"))
  if (!is.null(p_cal)) print(p_cal)

  # Metrics vs stratum variable (include all methods for completeness)
  sc_df <- .extract_scalars(eval_out, all_methods, stratum, metrics)
  p_sc  <- .plot_metrics_line(sc_df, x_label,
                               paste0(section_title, ": Metrics"))
  if (!is.null(p_sc)) print(p_sc)

  invisible(NULL)
}


# =============================================================================
# Main entry point
# =============================================================================

#' Plot benchmark evaluation results to a PDF
#'
#' Generates a multi-page PDF covering global and stratified (by S, phi, and
#' optionally p_causal) precision-recall curves, PIP calibration plots, and
#' metric summary panels with ±1 SE error bands/bars.
#'
#' @param eval_out  Output of \code{\link{evaluate_methods}}.
#' @param output_file Character or NULL. Full path of the PDF to write.
#'   Takes precedence over \code{output_dir} when specified. If NULL
#'   (default), the file is written as \code{evaluation.pdf} inside
#'   \code{output_dir}.
#' @param output_dir Character. Directory in which to save the PDF when
#'   \code{output_file} is NULL and \code{save = TRUE}. Created
#'   automatically if it does not exist. Default: \code{"results"}.
#' @param save Logical. If FALSE, the PDF is not written to disk (useful
#'   for checking that the function runs without producing a file).
#'   Default: TRUE.
#' @param methods Character vector or NULL. Methods to include. NULL = all
#'   evaluated methods.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return Invisibly returns the path of the PDF that was (or would have
#'   been) written.
#' @export
plot_results <- function(eval_out,
                          output_file = NULL,
                          output_dir  = "results",
                          save        = TRUE,
                          methods     = NULL,
                          verbose     = TRUE) {

  for (pkg in c("ggplot2", "gridExtra", "grid")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required. Install with install.packages('%s').",
                   pkg, pkg), call. = FALSE)
    }
  }

  # --- Resolve output path ----------------------------------------------------

  stopifnot(
    "output_dir must be a single character string" =
      is.character(output_dir) && length(output_dir) == 1L
  )

  if (is.null(output_file)) {
    output_file <- file.path(output_dir, "evaluation.pdf")
  }

  # --- Validate methods -------------------------------------------------------

  all_methods <- eval_out$methods_evaluated
  if (is.null(methods)) methods <- all_methods
  methods <- intersect(methods, all_methods)
  if (length(methods) == 0L) stop("No valid methods specified.", call. = FALSE)

  # Methods with at least some successful fits (non-NULL fdr_power_curve)
  valid_methods <- Filter(function(m) {
    !is.null(eval_out[[m]]$global$fdr_power_curve)
  }, methods)

  has_p_causal <- any(sapply(valid_methods, function(m) {
    !is.null(eval_out[[m]]$by_p_causal)
  }))

  if (!save) {
    if (verbose) message("  save = FALSE: skipping PDF output.")
    return(invisible(output_file))
  }

  # Create parent dir if needed
  out_dir <- dirname(output_file)
  if (nchar(out_dir) > 0L && out_dir != ".") {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  grDevices::pdf(output_file, width = 13, height = 9)
  on.exit(grDevices::dev.off(), add = TRUE)

  metrics <- names(.METRIC_LABELS)

  # ============================================================
  # Section 1 — Global
  # ============================================================
  if (verbose) message("  Plotting: Global section")

  # Page 1: PR curves
  pr_global <- .extract_pr(eval_out, valid_methods, "global")
  p1 <- .plot_pr_global(pr_global, "Global Precision-Recall Curves")
  if (!is.null(p1)) print(p1)

  # Page 2: PIP calibration
  cal_global <- .extract_cal(eval_out, valid_methods, "global")
  p2 <- .plot_cal_global(cal_global, "Global PIP Calibration")
  if (!is.null(p2)) print(p2)

  # Page 3: Summary table
  tbl <- .plot_summary_table(eval_out, methods)
  if (!is.null(tbl)) {
    grid::grid.newpage()
    grid::grid.draw(tbl)
  }

  # ============================================================
  # Section 2 — By S
  # ============================================================
  if (verbose) message("  Plotting: By S section")
  .print_section(eval_out, valid_methods, methods,
                 stratum        = "by_S",
                 stratum_label  = "S",
                 x_label        = "Number of Causal Variants (S)",
                 section_title  = "By S",
                 metrics        = metrics)

  # ============================================================
  # Section 3 — By phi
  # ============================================================
  if (verbose) message("  Plotting: By phi section")
  .print_section(eval_out, valid_methods, methods,
                 stratum        = "by_phi",
                 stratum_label  = "phi",
                 x_label        = "Proportion of Variance Explained (phi)",
                 section_title  = "By phi",
                 metrics        = metrics)

  # ============================================================
  # Section 4 — By p_causal (sparse_inf only)
  # ============================================================
  if (has_p_causal) {
    if (verbose) message("  Plotting: By p_causal section")
    .print_section(eval_out, valid_methods, methods,
                   stratum        = "by_p_causal",
                   stratum_label  = "p_causal",
                   x_label        = "Sparse Proportion (p_causal)",
                   section_title  = "By p_causal",
                   metrics        = metrics)
  }

  if (verbose) message(sprintf("  PDF saved: %s", output_file))
  invisible(output_file)
}
