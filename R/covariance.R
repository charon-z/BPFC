#' Construct the bivariate SAD(1) covariance of Method2
#'
#' No stationarity or spectral-radius constraint is imposed. Method2 specifies
#' a finite process starting at r1 = eta1. Visit gaps do not change A;
#' original times enter rho(t) only.
#' @param times Original nonnegative increasing times.
#' @param phi Two within-state lag parameters in (-1,1).
#' @param psi Two cross-state lag parameters in (-1,1), ordered A12, A21.
#' @param lambda Two positive parameters, without an ordering constraint.
#' @param v_sq,u_sq Positive innovation variances, one per observation time.
#' @return List A, G, Omega, R, rho in time-interleaved order (Eqs. 6-8).
#' @export
sad_covariance <- function(times, phi, psi, lambda, v_sq, u_sq) {
  T <- length(times)
  if (!is.numeric(times) || T < 2L || any(!is.finite(times)) ||
      any(times < 0) || any(diff(times) <= 0)) stop("Invalid original times.")
  for (a in list(phi, psi)) if (length(a) != 2L || any(!is.finite(a)) ||
    any(abs(a) >= 1)) stop("phi and psi require two values strictly inside (-1,1).")
  if (length(lambda) != 2L || any(!is.finite(lambda)) || any(lambda <= 0))
    stop("lambda requires two positive values.")
  for (v in list(v_sq, u_sq)) if (length(v) != T || any(!is.finite(v)) ||
    any(v <= 0)) stop("Innovation variances require T positive values.")
  A <- rbind(c(phi[1], psi[1]), c(psi[2], phi[2]))
  rho <- exp(-lambda[1]*times)-exp(-lambda[2]*times)
  if (any(abs(rho) >= 1)) stop("Innovation correlation is numerically singular.")
  G <- Omega <- matrix(0, 2L*T, 2L*T)
  powers <- list(diag(2))
  for (k in seq_len(T-1L)) powers[[k+1L]] <- powers[[k]] %*% A
  for (t in seq_len(T)) {
    ii <- c(2L*t-1L, 2L*t)
    off <- rho[t]*sqrt(v_sq[t]*u_sq[t])
    Omega[ii, ii] <- rbind(c(v_sq[t], off), c(off, u_sq[t]))
    for (s in seq_len(t)) G[ii, c(2L*s-1L, 2L*s)] <- powers[[t-s+1L]]
  }
  list(A = A, G = G, Omega = Omega, R = G %*% Omega %*% t(G), rho = rho)
}
