.bpfc_single_diagnostics <- function(samples) {
  ess <- if (requireNamespace("coda",quietly=TRUE)) vapply(seq_len(ncol(samples)),function(k) {
    x <- samples[,k]
    if (length(x)<10L || stats::var(x)==0) return(NA_real_)
    tryCatch(as.numeric(coda::effectiveSize(x)),error=function(e) NA_real_)
  },numeric(1)) else rep(NA_real_,ncol(samples))
  data.frame(parameter=colnames(samples),mean=colMeans(samples),
    sd=apply(samples,2L,stats::sd),ESS=ess,row.names=NULL)
}

#' Compare covariance convergence diagnostics across independent fits
#' @param fits List of at least two fits with identical model and time grid.
#' @return coda Gelman-Rubin PSRF and effective sample sizes for shared
#'   covariance parameters. Cluster means and weights require a separate
#'   across-chain ECR alignment and are deliberately excluded here.
#' @export
mcmc_diagnostics <- function(fits) {
  if (!requireNamespace("coda",quietly=TRUE)) stop("Install coda for chain diagnostics.")
  if (length(fits)<2L) stop("At least two independent fits are required for R-hat.")
  ref <- fits[[1]]
  for(f in fits) if (!identical(f$model_info$covariance,ref$model_info$covariance) ||
    !identical(f$model_info$J,ref$model_info$J) ||
    !identical(f$data_info$times_single,ref$data_info$times_single)) stop("Fits must share model and times.")
  cols <- grep("^(phi|psi|lambda|v_sq|u_sq)\\[",colnames(ref$posterior_samples$all_params),value=TRUE)
  chains <- lapply(fits,function(f) coda::mcmc(f$posterior_samples$all_params[,cols,drop=FALSE]))
  if (length(unique(vapply(chains,nrow,integer(1))))!=1L) stop("Chains require equal retained lengths.")
  ml <- coda::mcmc.list(chains)
  list(PSRF=coda::gelman.diag(ml,autoburnin=FALSE,multivariate=FALSE)$psrf,
       ESS=coda::effectiveSize(ml))
}
