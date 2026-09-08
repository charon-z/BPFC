#' Select a functional-module resolution in one overfitted-mixture run
#'
#' Fits one finite mixture at an upper bound `J_max` using a sparse symmetric
#' Dirichlet prior. The post-burn-in draws are summarized by the posterior
#' distribution of the occupied-component count and by the posterior-mean
#' mixing proportions, matching Equations (20)--(23) of the accompanying
#' manuscript. No grid of fixed-J models is fitted by this function.
#'
#' The returned `selected_K` is the smallest number of components whose sorted
#' posterior-mean weights reach `mass_threshold`. `K_occ_mode` is reported
#' separately because a large dataset can assign a few observations to a
#' negligible-weight component even under a sparse prior. The full `J_max` fit
#' and both diagnostics are retained, so the selection is auditable rather
#' than hidden behind a single number.
#'
#' @param x Data accepted by [run_mcmc_binary()].
#' @param J_max Integer upper bound on the number of mixture components.
#' @param e0 Positive symmetric Dirichlet concentration. Small values encourage
#'   redundant components to empty; the default reproduces the sparse setting
#'   used by the manuscript's overfitted-mixture analyses.
#' @param mass_threshold Number in `(0, 1]`. The working resolution is the
#'   smallest number of sorted posterior-mean component weights whose
#'   cumulative mass reaches this threshold.
#' @param occupied_min_count Minimum number of allocated features required for
#'   a component to count as occupied. The manuscript definition uses `1`.
#' @param credible_prob Probability used for the equal-tail credible interval
#'   of the occupied-component count.
#' @param ... Arguments passed to [run_mcmc_binary()], such as `times`, `niter`,
#'   `thin`, `burnin_frac`, `seed`, `two_phi`, and other prior settings through
#'   `priors`. Supplying `J` is an error; use `J_max`.
#'
#' @return An object of class `bpfc_overfit_result`. Important fields are
#'   `selected_K`, `K_occ_mode`, `K_occ_interval`, `occupied_posterior`,
#'   `component_spectrum`, `clustering`, and `fit`. `fit` is the complete
#'   `J_max` [run_mcmc_binary()] result. `clustering` reassigns features to the
#'   retained dominant components using their conditional posterior
#'   probabilities and numbers them by decreasing posterior mean weight.
#' @export
#' @examples
#' \donttest{
#' data(example_binary)
#' overfit <- run_overfitted_binary(
#'   example_binary$y, J_max = 6, times = example_binary$times,
#'   niter = 300, seed = 1
#' )
#' overfit$selected_K
#' overfit$occupied_posterior
#' plot_overfitted_selection(overfit)
#' }
run_overfitted_binary <- function(
    x,
    J_max,
    e0 = 0.02,
    mass_threshold = 0.90,
    occupied_min_count = 1L,
    credible_prob = 0.95,
    ...
) {
  if (!is.numeric(J_max) || length(J_max) != 1L || is.na(J_max) ||
      J_max != as.integer(J_max) || J_max < 2L) {
    stop("J_max must be one integer greater than or equal to 2.", call. = FALSE)
  }
  J_max <- as.integer(J_max)
  if (!is.numeric(e0) || length(e0) != 1L || !is.finite(e0) || e0 <= 0) {
    stop("e0 must be one positive finite number.", call. = FALSE)
  }
  .validate_overfit_summary_args(mass_threshold, occupied_min_count, credible_prob)

  dots <- list(...)
  if ("J" %in% names(dots)) {
    stop("Do not supply J; run_overfitted_binary() uses J_max.", call. = FALSE)
  }
  priors <- dots$priors
  if (is.null(priors)) priors <- list()
  if (!is.list(priors)) stop("priors must be a list.", call. = FALSE)
  priors$e0 <- e0
  dots$priors <- priors

  fit <- do.call(
    run_mcmc_binary,
    c(list(x = x, J = J_max), dots)
  )
  summary_fields <- .summarize_overfitted_fit(
    fit,
    mass_threshold = mass_threshold,
    occupied_min_count = occupied_min_count,
    credible_prob = credible_prob
  )
  result <- c(
    list(
      fit = fit,
      J_max = J_max,
      e0 = e0,
      mass_threshold = mass_threshold,
      occupied_min_count = as.integer(occupied_min_count),
      credible_prob = credible_prob,
      selection_rule = paste0(
        "smallest K with cumulative sorted posterior-mean mixing mass >= ",
        format(mass_threshold, trim = TRUE)
      )
    ),
    summary_fields
  )
  class(result) <- "bpfc_overfit_result"
  result
}

#' @keywords internal
#' @noRd
.validate_overfit_summary_args <- function(mass_threshold, occupied_min_count,
                                           credible_prob) {
  if (!is.numeric(mass_threshold) || length(mass_threshold) != 1L ||
      !is.finite(mass_threshold) || mass_threshold <= 0 || mass_threshold > 1) {
    stop("mass_threshold must be in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(occupied_min_count) || length(occupied_min_count) != 1L ||
      is.na(occupied_min_count) || occupied_min_count != as.integer(occupied_min_count) ||
      occupied_min_count < 1L) {
    stop("occupied_min_count must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(credible_prob) || length(credible_prob) != 1L ||
      !is.finite(credible_prob) || credible_prob <= 0 || credible_prob >= 1) {
    stop("credible_prob must be in (0, 1).", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.summarize_overfitted_fit <- function(fit, mass_threshold = 0.90,
                                      occupied_min_count = 1L,
                                      credible_prob = 0.95) {
  .validate_overfit_summary_args(mass_threshold, occupied_min_count, credible_prob)
  if (!inherits(fit, "mcmcgraph_result")) {
    stop("fit must be a mcmcgraph_result.", call. = FALSE)
  }
  z <- as.matrix(fit$posterior_samples$z)
  p <- as.matrix(fit$posterior_samples$p)
  J_max <- as.integer(fit$model_info$J)
  if (!nrow(z) || !ncol(z) || nrow(p) != nrow(z) || ncol(p) != J_max) {
    stop("fit does not contain compatible post-burn-in z and p draws.", call. = FALSE)
  }

  component_counts <- vapply(
    seq_len(J_max),
    function(j) rowSums(z == j),
    numeric(nrow(z))
  )
  if (J_max == 1L) component_counts <- matrix(component_counts, ncol = 1L)
  occupied_count <- rowSums(component_counts >= as.integer(occupied_min_count))
  occupied_tab <- tabulate(occupied_count, nbins = J_max)
  occupied_posterior <- data.frame(
    K = seq_len(J_max),
    probability = occupied_tab / length(occupied_count)
  )
  K_occ_mode <- which.max(occupied_tab)
  tail_prob <- (1 - credible_prob) / 2
  K_occ_interval <- as.integer(stats::quantile(
    occupied_count, probs = c(tail_prob, 1 - tail_prob),
    names = FALSE, type = 1
  ))
  names(K_occ_interval) <- c("lower", "upper")

  posterior_mean_weight <- unname(colMeans(p))
  component_order <- order(posterior_mean_weight, decreasing = TRUE)
  sorted_weight <- posterior_mean_weight[component_order]
  cumulative_mass <- cumsum(sorted_weight)
  selected_K <- which(cumulative_mass + sqrt(.Machine$double.eps) >= mass_threshold)[1L]
  dominant_components <- component_order[seq_len(selected_K)]

  map_counts <- tabulate(fit$clustering, nbins = J_max)
  component_spectrum <- data.frame(
    rank = seq_len(J_max),
    component = component_order,
    posterior_mean_weight = sorted_weight,
    cumulative_mass = cumulative_mass,
    map_count = map_counts[component_order],
    map_mass = map_counts[component_order] / sum(map_counts),
    dominant = seq_len(J_max) <= selected_K
  )

  raw_dominant_probability <- fit$cluster_prob[, dominant_components, drop = FALSE]
  retained_probability <- rowSums(raw_dominant_probability)
  conditional_probability <- raw_dominant_probability /
    pmax(retained_probability, .Machine$double.xmin)
  colnames(conditional_probability) <- paste0("cluster", seq_len(selected_K))
  clustering <- max.col(conditional_probability, ties.method = "first")
  cluster_uncertainty <- 1 - apply(conditional_probability, 1L, max)

  list(
    selected_K = as.integer(selected_K),
    K_occ_mode = as.integer(K_occ_mode),
    K_occ_interval = K_occ_interval,
    occupied_count = as.integer(occupied_count),
    occupied_posterior = occupied_posterior,
    component_spectrum = component_spectrum,
    dominant_components = as.integer(dominant_components),
    clustering = as.integer(clustering),
    cluster_prob = conditional_probability,
    cluster_uncertainty = cluster_uncertainty,
    retained_probability = retained_probability
  )
}

#' Print an overfitted BPFC selection result
#'
#' @param x A `bpfc_overfit_result`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.bpfc_overfit_result <- function(x, ...) {
  cat("BPFC single-run overfitted finite mixture\n")
  cat("  J_max:", x$J_max, "  e0:", format(x$e0), "\n")
  cat("  selected_K (mass rule):", x$selected_K,
      sprintf("[threshold %.1f%%]\n", 100 * x$mass_threshold))
  cat("  posterior mode K_occ:", x$K_occ_mode,
      sprintf("[%g%% interval %d--%d]\n", 100 * x$credible_prob,
              x$K_occ_interval[1L], x$K_occ_interval[2L]))
  invisible(x)
}

#' Plot diagnostics from one overfitted-mixture run
#'
#' Produces three base-R panels: the posterior distribution of the occupied
#' component count, sorted posterior-mean mixing weights, and their cumulative
#' mass. The working resolution selected by `mass_threshold` is highlighted.
#'
#' @param result A `bpfc_overfit_result` from [run_overfitted_binary()].
#' @param dominant_col,tail_col Colours for retained and tail components.
#' @return Invisibly, `result$component_spectrum`.
#' @export
plot_overfitted_selection <- function(result, dominant_col = "#D55E00",
                                      tail_col = "grey75") {
  if (!inherits(result, "bpfc_overfit_result")) {
    stop("result must be from run_overfitted_binary().", call. = FALSE)
  }
  spectrum <- result$component_spectrum
  posterior <- result$occupied_posterior
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.4, 0.8))

  graphics::barplot(
    posterior$probability, names.arg = posterior$K,
    col = ifelse(posterior$K == result$K_occ_mode, dominant_col, tail_col),
    border = NA, xlab = expression(K[occ]), ylab = "Posterior probability",
    main = "Occupied components"
  )

  component_col <- ifelse(spectrum$dominant, dominant_col, tail_col)
  graphics::barplot(
    spectrum$posterior_mean_weight, names.arg = spectrum$rank,
    col = component_col, border = NA,
    xlab = "Component rank", ylab = expression(bar(pi)[j]),
    main = "Sorted component mass"
  )

  graphics::plot(
    spectrum$rank, spectrum$cumulative_mass,
    type = "b", pch = 21, bg = component_col, col = "grey30",
    xlab = "Number of components retained", ylab = "Cumulative mass",
    ylim = c(0, 1), main = "Working resolution"
  )
  graphics::abline(h = result$mass_threshold, lty = 2, col = "grey45")
  graphics::abline(v = result$selected_K, lty = 2, col = dominant_col)
  graphics::points(result$selected_K,
                   spectrum$cumulative_mass[result$selected_K],
                   pch = 21, bg = dominant_col, cex = 1.4)
  invisible(spectrum)
}
