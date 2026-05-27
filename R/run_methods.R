# =============================================================================
# run_methods.R
#
# Unified entry point for applying fine-mapping methods to simulation output.
#
# Usage:
#   results <- run_methods(
#     simulation  = sim,
#     methods     = c("susie"),
#     method_args = list(susie = list(L = 5, coverage = 0.9)),
#     save        = TRUE,
#     output_dir  = "results"
#   )
#
# Adding a new method
# -------------------
# 1. Create R/wrappers/<method>.R implementing:
#      run_<method>_region(region_geno, region_pheno, ...)
#    The function must return a list with at minimum:
#      pip, credible_sets, method, input_type, params, runtime_seconds,
#      additional (list of method-specific outputs)
# 2. Add one entry to .FM_REGISTRY below.
# =============================================================================


# =============================================================================
# Method registry
# =============================================================================

# Maps a method name (lowercase) to the name of its region-level adapter
# function. The adapter must have the signature:
#   run_<method>_region(region_geno, region_pheno, ...)
#
# Methods may optionally define `run_<method>_scenario_setup(genotypes,
# regions, user_args)` to compute scenario-wide state (e.g. pooled
# per-annotation coefficients). When present, run_methods() calls it once
# per scenario and merges the returned list into method_args for that
# scenario's per-region calls.
.FM_REGISTRY <- list(
  susie     = "run_susie_region",
  susie_inf = "run_susie_inf_region",
  finemap   = "run_finemap_region",
  abf       = "run_abf_region",
  funmap    = "run_funmap_region",
  paintor   = "run_paintor_region",
  beatrice              = "run_beatrice_region",
  functional_beatrice   = "run_functional_beatrice_region",
  carma                 = "run_carma_region",
  # Baselines and annotation-aware comparators added in Phase 1.
  marginal_z            = "run_marginal_z_region",
  polyfun_oracle        = "run_polyfun_oracle_region",
  polyfun_est           = "run_polyfun_est_region",
  # Phase 4: modern variational comparator (without annotations).
  sparsepro             = "run_sparsepro_region"
)


# =============================================================================
# run_methods
# =============================================================================

#' Apply fine-mapping methods to simulation output
#'
#' Takes the output of \code{\link{run_simulation}} and applies one or more
#' fine-mapping methods to every scenario x region combination. Each method
#' is run with the user-supplied arguments (or defaults where not specified).
#'
#' @param simulation List. Output of \code{\link{run_simulation}}, containing
#'   \code{genotypes}, \code{scenarios}, and \code{params}.
#' @param methods Character vector. One or more method names. Currently
#'   supported: \code{"susie"}. Case-insensitive.
#' @param method_args Named list. Per-method argument overrides. Each element
#'   is itself a named list of arguments passed to the method's region runner.
#'   Arguments not listed here use the method's own defaults. Example:
#'   \code{list(susie = list(L = 5, coverage = 0.9))}.
#' @param save Logical. If TRUE, save results to \code{output_dir}. Each
#'   method's results are saved as \code{<method>.rds} inside a subdirectory
#'   named after the simulation parameters. A \code{run_metadata.rds} file
#'   is also written. Default: FALSE.
#' @param output_dir Character. Root directory for saved results when
#'   \code{save = TRUE}. Created automatically if it does not exist.
#'   Default: \code{"results"}.
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A list with one element per method (named by method), plus
#'   top-level metadata:
#' \describe{
#'   \item{<method>}{A list containing:
#'     \describe{
#'       \item{results}{Flat list of per-fit results (one per scenario x
#'         region). Each element is the standardised output of the method's
#'         region runner, plus metadata fields \code{scenario_id},
#'         \code{region_id}, \code{S}, \code{phi}, \code{p_causal},
#'         \code{iter}.}
#'       \item{method_args}{The argument list actually used (user overrides
#'         only; method defaults are embedded in each fit's \code{params}).}
#'       \item{n_total}{Integer. Total number of fits attempted.}
#'       \item{n_failed}{Integer. Number of fits that errored.}
#'       \item{total_runtime_seconds}{Numeric.}
#'     }
#'   }
#'   \item{methods_run}{Character vector of methods that were run.}
#'   \item{simulation_params}{The \code{params} list from the simulation.}
#'   \item{run_timestamp}{POSIXct. When the run started.}
#' }
#'
#' @examples
#' \dontrun{
#' sim <- run_simulation(n_regions = 2, n = 200, p = 100, seed = 1)
#'
#' # Run SuSiE with defaults
#' out <- run_methods(sim, methods = "susie")
#'
#' # Override SuSiE arguments
#' out <- run_methods(
#'   sim,
#'   methods     = "susie",
#'   method_args = list(susie = list(L = 5, coverage = 0.9))
#' )
#'
#' # Access results
#' out$susie$results[[1]]$pip
#' out$susie$results[[1]]$credible_sets
#' out$susie$results[[1]]$additional$alpha
#'
#' # Save to disk
#' out <- run_methods(sim, methods = "susie", save = TRUE, output_dir = "results")
#' }
#'
#' @export
run_methods <- function(simulation,
                        methods     = "susie",
                        method_args = list(),
                        save        = FALSE,
                        output_dir  = "results",
                        verbose     = TRUE) {

  # --- Validate inputs --------------------------------------------------------

  methods <- tolower(trimws(methods))

  unknown <- setdiff(methods, names(.FM_REGISTRY))
  if (length(unknown) > 0) {
    stop(sprintf(
      "Unknown method(s): %s.\nAvailable methods: %s",
      paste(unknown, collapse = ", "),
      paste(sort(names(.FM_REGISTRY)), collapse = ", ")
    ), call. = FALSE)
  }

  if (!is.list(method_args)) {
    stop("method_args must be a named list (e.g. list(susie = list(L = 5))).",
         call. = FALSE)
  }

  unrecognised_arg_keys <- setdiff(names(method_args), methods)
  if (length(unrecognised_arg_keys) > 0) {
    warning(sprintf(
      "method_args contains entries for method(s) not being run: %s",
      paste(unrecognised_arg_keys, collapse = ", ")
    ))
  }

  if (is.null(simulation$scenarios) || is.null(simulation$genotypes)) {
    stop(
      "simulation must be the output of run_simulation(), containing ",
      "'genotypes' and 'scenarios'.",
      call. = FALSE
    )
  }

  n_scenarios <- length(simulation$scenarios)
  n_regions   <- length(simulation$genotypes)
  n_total     <- n_scenarios * n_regions

  run_timestamp <- Sys.time()

  # --- Run each method --------------------------------------------------------

  output <- vector("list", length(methods))
  names(output) <- methods

  for (method in methods) {

    run_fn    <- match.fun(.FM_REGISTRY[[method]])
    user_args <- if (!is.null(method_args[[method]])) method_args[[method]] else list()

    if (verbose) {
      message(sprintf(
        "\n=== %s | %d scenario(s) x %d region(s) = %d fit(s) ===",
        toupper(method), n_scenarios, n_regions, n_total
      ))
      if (length(user_args) > 0) {
        message(sprintf(
          "    Args: %s",
          paste(names(user_args), user_args, sep = " = ", collapse = ", ")
        ))
      }
    }

    # Optional scenario-level setup hook.
    # If a method defines `run_<method>_scenario_setup(genotypes, regions,
    # user_args)`, it is called once per scenario before that scenario's
    # per-region calls. The hook returns a named list whose entries are
    # merged into `user_args` for every region of that scenario. This lets
    # methods compute shared state (e.g. pooled per-annotation coefficients
    # in polyfun_est) without needing to see all regions inside the
    # per-region wrapper.
    #
    # Methods that do not define this hook behave exactly as before.
    scenario_setup_name <- paste0("run_", method, "_scenario_setup")
    has_scenario_setup  <- exists(scenario_setup_name, mode = "function")
    if (has_scenario_setup) {
      scenario_setup_fn <- match.fun(scenario_setup_name)
    }

    batch_start    <- proc.time()
    method_results <- vector("list", n_total)
    n_failed       <- 0L
    idx            <- 0L

    for (sc in seq_len(n_scenarios)) {
      scenario <- simulation$scenarios[[sc]]

      # Run the scenario-level setup hook (if any) to obtain extra args.
      effective_args <- user_args
      if (has_scenario_setup) {
        extra_args <- tryCatch(
          scenario_setup_fn(
            genotypes = simulation$genotypes,
            regions   = scenario$regions,
            user_args = user_args
          ),
          error = function(e) {
            if (verbose) {
              message(sprintf(
                "    WARNING: scenario_setup for %s failed (scenario %d): %s",
                method, scenario$scenario_id, conditionMessage(e)
              ))
            }
            list()
          }
        )
        if (is.list(extra_args) && length(extra_args) > 0L) {
          effective_args <- utils::modifyList(user_args, extra_args)
        }
      }

      for (rg in seq_len(n_regions)) {
        idx <- idx + 1L

        if (verbose) {
          message(sprintf(
            "  [%d/%d] Scenario %d (S=%d, phi=%.2f, iter=%d), Region %d",
            idx, n_total,
            scenario$scenario_id, scenario$S, scenario$phi, scenario$iter, rg
          ))
        }

        region_geno  <- simulation$genotypes[[rg]]
        region_pheno <- scenario$regions[[rg]]

        fit <- tryCatch(
          do.call(run_fn, c(
            list(region_geno = region_geno, region_pheno = region_pheno),
            effective_args
          )),
          error = function(e) list(
            pip             = rep(NA_real_, nrow(region_geno$LD)),
            credible_sets   = list(),
            method          = method,
            input_type      = NA_character_,
            params          = effective_args,
            runtime_seconds = NA_real_,
            additional      = list(),
            error           = conditionMessage(e)
          )
        )

        if (!is.null(fit$error)) {
          n_failed <- n_failed + 1L
          if (verbose) {
            message(sprintf("    WARNING: %s failed — %s", method, fit$error))
          }
        }

        # Attach scenario / region metadata
        fit$scenario_id <- scenario$scenario_id
        fit$region_id   <- rg
        fit$S           <- scenario$S
        fit$phi         <- scenario$phi
        fit$p_causal    <- scenario$p_causal
        fit$iter        <- scenario$iter

        method_results[[idx]] <- fit
      }
    }

    batch_elapsed <- as.numeric((proc.time() - batch_start)["elapsed"])

    if (verbose) {
      message(sprintf(
        "    %d/%d fits successful (%.1f s total).",
        n_total - n_failed, n_total, batch_elapsed
      ))
    }

    output[[method]] <- list(
      results               = method_results,
      method_args           = user_args,
      n_total               = n_total,
      n_failed              = n_failed,
      total_runtime_seconds = batch_elapsed
    )
  }

  # --- Assemble top-level result ----------------------------------------------

  result <- c(
    output,
    list(
      methods_run       = methods,
      simulation_params = simulation$params,
      run_timestamp     = run_timestamp
    )
  )

  # --- Save to disk -----------------------------------------------------------

  if (save) {
    stopifnot(
      "output_dir must be a single character string" =
        is.character(output_dir) && length(output_dir) == 1
    )

    sim_label <- .make_sim_label(simulation$params)
    run_dir   <- file.path(output_dir, sim_label)

    if (!dir.exists(run_dir)) {
      dir.create(run_dir, recursive = TRUE)
    }

    for (method in methods) {
      fpath <- file.path(run_dir, paste0(method, ".rds"))
      saveRDS(result[[method]], file = fpath)
      if (verbose) message(sprintf("    Saved %s results: %s", method, fpath))
    }

    meta <- list(
      methods_run       = methods,
      method_args       = method_args,
      simulation_params = simulation$params,
      run_timestamp     = run_timestamp
    )
    meta_path <- file.path(run_dir, "run_metadata.rds")
    saveRDS(meta, file = meta_path)
    if (verbose) message(sprintf("    Saved run metadata: %s", meta_path))
  }

  result
}


# =============================================================================
# Internal helpers
# =============================================================================

# Build a human-readable label from simulation params for use as a directory
# name. Mirrors the filename convention used in run_simulation().
.make_sim_label <- function(params) {
  seed_tag <- if (!is.null(params$seed)) paste0("seed", params$seed) else "noseed"
  p_vals   <- params$p %||% "p?"
  p_tag    <- if (length(unique(p_vals)) == 1) as.character(p_vals[1]) else
                paste(p_vals, collapse = "-")
  S_tag    <- paste(params$S_values, collapse = "-")

  sprintf(
    "%s_n%d_p%s_S%s_iter%d_%s",
    params$model %||% "model",
    params$n     %||% 0L,
    p_tag,
    S_tag,
    params$n_iter %||% 0L,
    seed_tag
  )
}

# %||% is defined in R/utils.R.
