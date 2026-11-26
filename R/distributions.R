# Method2 Eqs. (5)-(7), (19). Time-interleaved order throughout.
# The innovation transform has unit determinant, so no extra Jacobian.
dBivSAD <- nimble::nimbleFunction(
  run = function(x = double(1), mean = double(1), phi = double(1),
                 psi = double(1), lambda = double(1), v_sq = double(1),
                 u_sq = double(1), times = double(1),
                 log = integer(0, default = 0L)) {
    returnType(double(0))
    q <- 0.0
    prev1 <- 0.0
    prev2 <- 0.0
    for (t in 1:length(times)) {
      if (v_sq[t] <= 0.0 | u_sq[t] <= 0.0) {
        if (log) return(-Inf) else return(0.0)
      }
      rho <- exp(-lambda[1]*times[t]) - exp(-lambda[2]*times[t])
      h <- 1.0 - rho*rho
      if (h <= 0.0) {
        if (log) return(-Inf) else return(0.0)
      }
      r1 <- x[2*t-1] - mean[2*t-1]
      r2 <- x[2*t] - mean[2*t]
      e1 <- r1 - phi[1]*prev1 - psi[1]*prev2
      e2 <- r2 - psi[2]*prev1 - phi[2]*prev2
      a <- e1/sqrt(v_sq[t])
      b <- e2/sqrt(u_sq[t])
      q <- q + 2.0*log(2.0*3.141592653589793) + log(v_sq[t]) +
        log(u_sq[t]) + log(h) + (a*a - 2.0*rho*a*b + b*b)/h
      prev1 <- r1
      prev2 <- r2
    }
    ans <- -0.5*q
    if (log) return(ans) else return(exp(ans))
  }
)

rBivSAD <- nimble::nimbleFunction(
  run = function(n = integer(0, default = 1L), mean = double(1),
                 phi = double(1), psi = double(1), lambda = double(1),
                 v_sq = double(1), u_sq = double(1), times = double(1)) {
    returnType(double(1))
    if (n != 1L) stop("rBivSAD supports n=1 only")
    out <- numeric(2*length(times), init = TRUE)
    prev1 <- 0.0
    prev2 <- 0.0
    for (t in 1:length(times)) {
      rho <- exp(-lambda[1]*times[t]) - exp(-lambda[2]*times[t])
      a <- rnorm(1, 0, 1)
      b <- rnorm(1, 0, 1)
      r1 <- phi[1]*prev1 + psi[1]*prev2 + sqrt(v_sq[t])*a
      r2 <- psi[2]*prev1 + phi[2]*prev2 +
        sqrt(u_sq[t])*(rho*a + sqrt(1.0-rho*rho)*b)
      out[2*t-1] <- mean[2*t-1] + r1
      out[2*t] <- mean[2*t] + r2
      prev1 <- r1
      prev2 <- r2
    }
    return(out)
  }
)
