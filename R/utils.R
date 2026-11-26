#' Prepare complete paired observations in Method2 time order
#' @param x Numeric wide matrix, or long data frame.
#' @param format Input format, detected from long column names when auto.
#' @param times Original nonnegative strictly increasing observation times.
#'   Required explicitly for wide input; never replaced by visit indices.
#' @param layout Wide column order: interleaved (Eq. 1), or concat for old
#'   BPFC matrices with state1(all times), state2(all times).
#' @param id_col,ind_col,time_col,value_col Long input column names. State
#'   order follows factor levels, or sorted character state names.
#' @return List y (interleaved), d_single, times_single, state_levels,
#'   input_layout, format.
#' @export
as_binary_data <- function(x, format = c("auto", "wide", "long"), times = NULL,
                           layout = c("interleaved", "concat"), id_col = "snp",
                           ind_col = "ind", time_col = "time", value_col = "y") {
  format <- match.arg(format); layout <- match.arg(layout)
  if (format == "auto") format <- if (is.data.frame(x) &&
    all(c(ind_col, time_col, value_col) %in% names(x))) "long" else "wide"
  state_levels <- c("state1", "state2")
  if (format == "long") {
    need <- c(id_col, ind_col, time_col, value_col)
    if (!all(need %in% names(x)) || !nrow(x)) stop("Invalid long input columns.")
    if (anyNA(x[, need])) stop("Long input contains NA.")
    if (!is.numeric(x[[time_col]]) || !is.numeric(x[[value_col]]))
      stop("Long times and observations must be numeric.")
    if (anyDuplicated(x[, c(id_col, ind_col, time_col)]))
      stop("Duplicate feature/state/time observations are not supported.")
    states <- x[[ind_col]]
    state_levels <- if (is.factor(states)) levels(droplevels(states)) else
      sort(unique(as.character(states)))
    if (length(state_levels) != 2L) stop("Exactly two states are required.")
    observed_times <- sort(unique(x[[time_col]]))
    if (is.null(times)) times <- observed_times
    if (!identical(as.numeric(times), as.numeric(observed_times)))
      stop("times must match all long-format observation times.")
    ids <- unique(x[[id_col]])
    y <- matrix(NA_real_, length(ids), 2L*length(times),
                dimnames = list(as.character(ids), NULL))
    rows <- match(x[[id_col]], ids)
    cols <- 2L*(match(x[[time_col]], times)-1L) + match(as.character(states), state_levels)
    y[cbind(rows, cols)] <- x[[value_col]]
  } else {
    y <- as.matrix(x)
    if (!is.numeric(y)) stop("Wide observations must be numeric.")
    if (ncol(y) < 4L || ncol(y) %% 2L) stop("Wide input requires 2*T columns, T >= 2.")
    if (is.null(times)) stop("Supply original observation times explicitly.")
    if (length(times) != ncol(y)/2L) stop("times length must equal T.")
    if (layout == "concat") y <- y[, as.vector(rbind(seq_along(times),
      length(times)+seq_along(times))), drop = FALSE]
  }
  if (!is.numeric(times) || length(times) < 2L || any(!is.finite(times)) ||
      any(times < 0) || any(diff(times) <= 0))
    stop("times must be finite, nonnegative and strictly increasing, with T >= 2.")
  if (!nrow(y) || any(!is.finite(y))) stop("Complete finite observations are required (no NA/Inf).")
  storage.mode(y) <- "double"
  list(y = y, d_single = as.integer(length(times)), times_single = as.numeric(times),
       state_levels = state_levels, input_layout = if (format == "long") "long" else layout,
       format = format)
}

calc_Z_legendre <- function(times, order = 4L) {
  if (order != 4L) stop("Method2 fixes polynomial order at four.")
  x <- 2*(times-min(times))/(max(times)-min(times))-1
  cbind(1, x, (3*x^2-1)/2, (5*x^3-3*x)/2, (35*x^4-30*x^2+3)/8)
}
make_Z0_binary <- function(times_single) {
  Z <- calc_Z_legendre(times_single)
  B <- matrix(0, 2L*nrow(Z), 10L)
  B[seq.int(1L, nrow(B), 2L), 1:5] <- Z
  B[seq.int(2L, nrow(B), 2L), 6:10] <- Z
  B
}
log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x-m)))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
