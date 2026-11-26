#' Generate paired observations from the exact Method2 model
#'
#' This utility is for later simulation work; no scientific scenario is preset.
#' @param n Number of features.
#' @param times Original observation times.
#' @param beta J by 10 Legendre coefficients: five for each state.
#' @param p Positive mixture weights summing to one.
#' @param phi,psi,lambda,v_sq,u_sq Covariance parameters as in [sad_covariance()].
#' @param seed Random seed.
#' @return List y (time-interleaved), times, truth (labels and all generating
#'   parameters), layout. The first residual is the first innovation.
#' @export
simulate_bpfc <- function(n,times,beta,p,phi,psi,lambda,v_sq,u_sq,seed=123L) {
  if(length(n)!=1L || !is.finite(n) || n<1 || n!=floor(n)) stop("n must be a positive integer.")
  if(!is.matrix(beta) || ncol(beta)!=10L || !nrow(beta) || any(!is.finite(beta)))
    stop("beta requires a finite J x 10 matrix.")
  J <- nrow(beta)
  if(length(p)!=J || any(!is.finite(p)) || any(p<=0) || abs(sum(p)-1)>1e-8)
    stop("p requires J positive weights summing to one.")
  cv <- sad_covariance(times,phi,psi,lambda,v_sq,u_sq)
  mu <- beta %*% t(make_Z0_binary(times))
  set.seed(seed)
  z <- sample.int(J,n,replace=TRUE,prob=p)
  y <- t(vapply(z,function(j) rBivSAD(1L,mu[j,],phi,psi,lambda,v_sq,u_sq,times),
                numeric(2L*length(times))))
  list(y=y,times=times,layout="interleaved",
       truth=list(J=J,z=z,beta=beta,p=p,phi=phi,psi=psi,lambda=lambda,
                  v_sq=v_sq,u_sq=u_sq,R=cv$R),seed=seed)
}
