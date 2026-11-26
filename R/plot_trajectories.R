#' Plot paired posterior mean functional-module trajectories
#' @param result A Method2 mcmcgraph_result.
#' @param components Components to plot; default all.
#' @param ncol Number of panel columns.
#' @param colours Two state colours.
#' @return Invisibly the paired posterior mean matrix, in Eq. 3 time order.
#' @export
plot_paired_trajectories <- function(result,components=NULL,ncol=3L,
                                     colours=c("#0072B2","#D55E00")) {
  if (!identical(result$model_info$covariance,"BIVARIATE_SAD1_GOmegaGt"))
    stop("A Method2 fit is required.")
  J <- result$model_info$J
  if(is.null(components)) components <- seq_len(J)
  if(!length(components) || any(!is.finite(components)) ||
     any(components!=floor(components)) || any(components<1|components>J))
    stop("Invalid components.")
  if(length(ncol)!=1L || !is.finite(ncol) || ncol<1 || ncol!=floor(ncol)) stop("Invalid ncol.")
  mu <- result$posterior_mean$mu[components,,drop=FALSE]
  times <- result$data_info$times_single; TT <- length(times)
  old <- graphics::par(no.readonly=TRUE);on.exit(graphics::par(old),add=TRUE)
  graphics::par(mfrow=c(ceiling(length(components)/ncol),ncol))
  for(k in seq_along(components)) {
    yy <- cbind(mu[k,seq.int(1,2*TT,2)],mu[k,seq.int(2,2*TT,2)])
    graphics::matplot(times,yy,type="l",lty=1:2,col=colours,
                     xlab="Original time",ylab="Posterior mean",
                     main=paste("Component",components[k]))
    graphics::legend("topright",legend=result$data_info$state_levels,
                     col=colours,lty=1:2,bty="n",cex=.8)
  }
  invisible(mu)
}
