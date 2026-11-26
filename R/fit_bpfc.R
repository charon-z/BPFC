#' Fit paired trajectories with a choice of cluster-count selection
#'
#' `selection = "automatic"` runs a sparse overfitted mixture once at `J_max`
#' and retains components by their cumulative posterior mean weight.
#' `selection = "bic"` fits each value in `J_grid` and selects the smallest
#' posterior-mean plug-in BIC. The same supplied observations are used throughout.
#'
#' @param x Data accepted by [run_mcmc_binary()].
#' @param selection Either `"automatic"` (default) or `"bic"`.
#' @param J_max Upper bound for automatic selection. Required in automatic mode.
#' @param J_grid Candidate positive integer cluster counts for BIC selection.
#' @param auto_control Named list of settings for automatic selection:
#'   `e0`, `mass_threshold`, `occupied_min_count`, and `credible_prob`.
#'   e0 and mass_threshold are required; the remaining settings have the
#'   defaults documented in [run_overfitted_binary()].
#' @param ... Arguments to [run_mcmc_binary()], including `times`, `niter`,
#'   `seed`, `layout`, and `priors`. Do not supply `J`.
#' @return A `bpfc_selection_result` with `selection`, `selected_K`,
#'   `clustering`, `cluster_prob`, `cluster_uncertainty`, `fit`, and
#'   `selection_details`. In automatic mode, `fit` is the full `J_max` fit;
#'   the top-level allocations and probabilities retain all J_max components;
#'   selecting dominant mass does not reassign tail features.
#'   In BIC mode, `fit` is the winning fixed-J fit; `selection_details`
#'   retains the BIC table and all candidate fits. Exact BIC ties favor smaller J.
#' @export
fit_bpfc <- function(x, selection = c("automatic", "bic"), J_max = NULL,
                     J_grid = NULL, auto_control = list(), ...) {
  selection <- match.arg(selection)
  dots <- list(...)
  if ("J" %in% names(dots)) {
    stop("Use J_max or J_grid for selection; use run_mcmc_binary() for a fixed J.",
         call. = FALSE)
  }
  if (!is.list(auto_control) || (length(auto_control) &&
      (is.null(names(auto_control)) || any(!nzchar(names(auto_control))) ||
       anyDuplicated(names(auto_control))))) {
    stop("auto_control must be a named list with unique names.", call. = FALSE)
  }
  controls <- c("e0", "mass_threshold", "occupied_min_count", "credible_prob")
  if (any(!names(auto_control) %in% controls)) {
    stop("Unknown auto_control option.", call. = FALSE)
  }
  if (selection == "automatic") {
    if (is.null(J_max)) stop("Set J_max for automatic selection.", call. = FALSE)
    if (!is.null(J_grid)) stop("J_grid is used only with selection = 'bic'.", call. = FALSE)
    details <- do.call(run_overfitted_binary, c(list(x = x, J_max = J_max), auto_control, dots))
    fit <- details$fit
    selected_K <- details$selected_K
    allocations <- details
  } else {
    if (!is.null(J_max) || length(auto_control)) {
      stop("J_max and auto_control are used only with selection = 'automatic'.", call. = FALSE)
    }
    if (is.null(J_grid)) stop("Set J_grid for BIC selection.", call. = FALSE)
    details <- do.call(fit_many_J, c(list(x = x, J_grid = J_grid), dots))
    if (any(!is.finite(details$eval$BIC))) {
      stop("A candidate BIC is not finite; inspect the candidate fits.", call. = FALSE)
    }
    winner <- order(details$eval$BIC, details$eval$J)[1L]
    fit <- details$fits[[winner]]
    selected_K <- as.integer(details$eval$J[winner])
    details$selected_row <- winner
    allocations <- fit
  }
  result <- list(selection = selection, selected_K = selected_K,
                 clustering = allocations$clustering,
                 cluster_prob = allocations$cluster_prob,
                 cluster_confidence = 1-allocations$cluster_uncertainty,
                 cluster_uncertainty = allocations$cluster_uncertainty,
                 fit = fit, selection_details = details)
  class(result) <- "bpfc_selection_result"
  result
}

#' Print a BPFC cluster-count selection result
#' @param x A `bpfc_selection_result`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.bpfc_selection_result <- function(x, ...) {
  cat("BPFC cluster-count selection\n")
  cat("  Method:", x$selection, "\n")
  cat("  Selected K:", x$selected_K, "\n")
  if (x$selection == "bic") print(x$selection_details$eval, row.names = FALSE)
  invisible(x)
}
