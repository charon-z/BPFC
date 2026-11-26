#' Refit a grid of fixed-J models for supplementary BIC comparison
#' @param x Paired observations.
#' @param J_grid Positive integer component counts.
#' @param ... Arguments to [run_mcmc_binary()].
#' @return List of fits and a table of posterior point-estimate BIC values.
#' @export
fit_many_J <- function(x, J_grid=2:6, ...) {
  if (!is.numeric(J_grid) || !length(J_grid) || any(!is.finite(J_grid)) ||
      any(J_grid < 1 | J_grid != floor(J_grid))) stop("J_grid requires positive integers.")
  J_grid <- sort(unique(as.integer(J_grid)))
  fits <- lapply(J_grid,function(j) run_mcmc_binary(x,J=j,...))
  ev <- lapply(fits,function(f) eval_bic(f,x,format=f$data_info$format,
                                        layout=f$data_info$input_layout))
  names(fits) <- paste0("J",J_grid)
  list(fits=fits,eval=do.call(rbind,ev))
}

#' Evaluate Method2 BIC at ECR-aligned posterior means
#' @param fit A new Method2 fit. Old covariance fits must be refitted.
#' @param x The original observations, in the original layout.
#' @param format Input format; defaults to the fitted format.
#' @param layout Input layout; defaults to the fitted layout.
#' @return Data frame J, loglik, q, BIC. q = 10*J + (J-1) + 2*T + 6.
#'   The observed mixture likelihood is used, not the complete-data likelihood.
#' @export
eval_bic <- function(fit,x,format=NULL,layout=NULL) {
  if (!identical(fit$model_info$covariance,"BIVARIATE_SAD1_GOmegaGt"))
    stop("Refit with the Method2 bivariate covariance before evaluating BIC.")
  format <- format %||% fit$data_info$format
  layout <- layout %||% fit$data_info$input_layout
  dat <- as_binary_data(x,format=format,times=fit$data_info$times_single,
                        layout=if(layout=="long") "interleaved" else layout)
  y <- dat$y
  if (nrow(y)!=fit$model_info$n || dat$d_single!=fit$model_info$d_single)
    stop("BIC requires the original fitted data dimensions.")
  pm <- fit$posterior_mean; J <- fit$model_info$J
  mu <- pm$beta %*% t(make_Z0_binary(dat$times_single))
  loglik <- sum(vapply(seq_len(nrow(y)),function(i) {
    lc <- vapply(seq_len(J),function(j) log(pm$p[j])+dBivSAD(y[i,],mu[j,],
      pm$phi,pm$psi,pm$lambda,pm$v_sq,pm$u_sq,dat$times_single,log=1L),numeric(1))
    log_sum_exp(lc)
  },numeric(1)))
  q <- 10L*J + J-1L + 2L*dat$d_single + 6L
  data.frame(J=J,loglik=loglik,q=q,BIC=-2*loglik+q*log(nrow(y)))
}
