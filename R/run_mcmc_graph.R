#' Fit the bivariate SAD(1) paired finite mixture in Method2
#'
#' @param x Paired data accepted by [as_binary_data()].
#' @param J Fixed component count, between one and the number of features.
#' @param format Input format.
#' @param times Original observation times; required for wide data.
#' @param layout Wide input order; select concat explicitly for old BPFC data.
#' @param niter Total iterations including burn-in.
#' @param thin Posterior thinning interval, applied after burn-in.
#' @param burnin_frac Fraction of iterations discarded; default 25 percent.
#' @param priors Named list: alpha (scalar or J), beta_mean (length 10 or
#'   J by 10), beta_cov (10 by 10 or J by 10 by 10), a_v,b_v,a_u,b_u
#'   (scalar or T), mu_phi,var_phi,mu_psi,var_psi (scalars shared within
#'   each prior family), a_lambda,b_lambda (scalar or length 2). Gamma uses shape/rate; inverse-Gamma has density
#'   proportional to x^(-a-1)*exp(-b/x). ALL settings must be supplied
#'   explicitly: Method2 specifies prior families but no numerical values.
#' @param init Optional list; z alone is sufficient. Complete or partial
#'   overrides may include beta, p, v_sq, u_sq, phi, psi, lambda.
#' @param seed Random seed.
#' @return mcmcgraph_result with ECR-aligned samples, full posterior assignment
#'   probabilities, confidence, paired posterior mean trajectories, covariance,
#'   sampler configuration and numerical diagnostics. Short runs verify software
#'   operation and do not establish convergence.
#' @export
run_mcmc_binary <- function(x, J, format = c("auto", "wide", "long"),
                            times = NULL, layout = c("interleaved", "concat"),
                            niter = 3000L, thin = 1L, burnin_frac = 0.25,
                            priors, init = NULL, seed = 123L) {
  format <- match.arg(format); layout <- match.arg(layout)
  dat <- as_binary_data(x, format, times, layout)
  y <- dat$y; n <- nrow(y); TT <- dat$d_single; d <- ncol(y); P <- 10L
  check_int <- function(a, lo, nm) {
    if (length(a) != 1L || !is.finite(a) || a != floor(a) || a < lo)
      stop(nm, " must be an integer >= ", lo, call. = FALSE)
  }
  check_int(J, 1, "J"); check_int(niter, 2, "niter"); check_int(thin, 1, "thin")
  if (J > n) stop("J cannot exceed the feature count.")
  if (length(burnin_frac) != 1L || !is.finite(burnin_frac) ||
      burnin_frac < 0 || burnin_frac >= 1) stop("burnin_frac must lie in [0,1).")
  nburnin <- floor(niter*burnin_frac)
  if (floor((niter-nburnin)/thin) < 2L) stop("Retain at least two posterior samples.")
  J <- as.integer(J); niter <- as.integer(niter); thin <- as.integer(thin)
  pr <- .bpfc_priors(priors, J, TT)
  set.seed(seed)
  Z0 <- make_Z0_binary(dat$times_single)
  ini <- .bpfc_init(y, J, Z0, dat$times_single, init)
  .bpfc_register_bivariate()
  code <- nimble::nimbleCode({
    if (J > 1) p[1:J] ~ ddirch(alpha[1:J,1])
    for (j in 1:J) {
      beta[j, 1:P] ~ dmnorm(beta_mean[j, 1:P], cov = beta_cov[j, 1:P, 1:P])
      for (t in 1:d) mu[j,t] <- inprod(beta[j,1:P], Z0[t,1:P])
    }
    for (t in 1:TT) {
      v_sq[t] ~ dinvgamma(a_v[t], scale=b_v[t])
      u_sq[t] ~ dinvgamma(a_u[t], scale=b_u[t])
    }
    for (k in 1:2) {
      phi[k] ~ T(dnorm(mu_phi[k], var = var_phi[k]), -1, 1)
      psi[k] ~ T(dnorm(mu_psi[k], var = var_psi[k]), -1, 1)
      lambda[k] ~ dgamma(a_lambda[k], rate = b_lambda[k])
    }
    for (i in 1:n) {
      if (J > 1) {
        z[i] ~ dcat(p[1:J])
        y[i,1:d] ~ dBivSAD(mu[z[i],1:d], phi[1:2], psi[1:2], lambda[1:2],
                           v_sq[1:TT], u_sq[1:TT], original_times[1:TT])
      } else {
        # Avoid NIMBLE's unsupported dynamic indexing into a singleton axis.
        y[i,1:d] ~ dBivSAD(mu[1,1:d], phi[1:2], psi[1:2], lambda[1:2],
                           v_sq[1:TT], u_sq[1:TT], original_times[1:TT])
      }
    }
  })
  constants <- c(list(J=J, P=P, n=n, TT=TT, d=d, Z0=Z0,
                     original_times=dat$times_single), pr)
  constants$alpha <- matrix(pr$alpha, nrow=J, ncol=1L)
  dims <- list(alpha=c(J,1L),beta_mean=c(J,P),beta_cov=c(J,P,P),
               beta=c(J,P),mu=c(J,d),Z0=c(d,P))
  if (J > 1L) { dims$p <- J; dims$z <- n } else {
    # The one-component Dirichlet/categorical distributions are point masses.
    constants$alpha <- NULL; dims$alpha <- NULL; ini$p <- NULL; ini$z <- NULL
  }
  model <- nimble::nimbleModel(code, constants=constants, data=list(y=y),
                               inits=ini, check=FALSE, dimensions=dims)
  if (!is.finite(model$calculate())) stop("Non-finite initial model density.")
  # Preserve NIMBLE's identified conjugate Dirichlet update, replace all others.
  monitors <- c("beta","v_sq","u_sq","phi","psi","lambda")
  if (J > 1L) monitors <- c("z","p",monitors)
  conf <- nimble::configureMCMC(model, print=FALSE,monitors=monitors)
  if (J > 1L) {
    samplers <- conf$getSamplers()
    p_sampler <- which(vapply(samplers, function(s) any(grepl("^p\\[", s$target)), logical(1)))
    if (length(p_sampler) != 1L || !grepl("conjugate", samplers[[p_sampler]]$name,
                                        ignore.case=TRUE))
      stop("NIMBLE did not identify the conjugate Dirichlet update.")
    p_target <- samplers[[p_sampler]]$target
  }
  conf$removeSamplers()
  if (J > 1L) {
    for (i in seq_len(n)) conf$addSampler(sprintf("z[%d]",i), type="categorical")
    conf$addSampler(p_target, type="conjugate")
  }
  for (j in seq_len(J)) conf$addSampler(sprintf("beta[%d, 1:%d]",j,P),
    type="RW_block", control=list(adaptive=TRUE, scale=0.05, adaptInterval=200))
  for (nm in c("v_sq", "u_sq")) for (t in seq_len(TT))
    conf$addSampler(sprintf("%s[%d]",nm,t), type="slice")
  for (nm in c("phi","psi","lambda")) for (k in 1:2)
    conf$addSampler(sprintf("%s[%d]",nm,k), type="slice")
  # Exact Method2 sequence: z -> p -> beta -> variances -> dependence.
  sampler_info <- lapply(conf$getSamplers(), function(s) list(target=s$target, type=s$name))
  cmodel <- nimble::compileNimble(model, showCompilerOutput=FALSE)
  mcmc <- nimble::buildMCMC(conf)
  cmcmc <- nimble::compileNimble(mcmc, project=model, showCompilerOutput=FALSE)
  # A single continuous chain preserves adaptation and RNG state. No restarts.
  post <- as.matrix(nimble::runMCMC(cmcmc, niter=niter, nburnin=nburnin,
    thin=thin, progressBar=interactive(), samplesAsCodaMCMC=FALSE))
  if (any(!is.finite(post))) stop("Non-finite posterior sample.")
  out <- .bpfc_summarize(post, J, P, Z0, n)
  pm <- out$posterior_mean
  out$posterior_mean$covariance <- sad_covariance(dat$times_single,
    pm$phi, pm$psi, pm$lambda, pm$v_sq, pm$u_sq)$R
  out$model_info <- list(J=J,n=n,d_single=TT,d_binary=d,order=4L,P=P,
    niter=niter,thin=thin,burnin_frac=burnin_frac,nburnin=nburnin,
    retained_iterations=nburnin+seq_len(nrow(post))*thin,
    covariance="BIVARIATE_SAD1_GOmegaGt",method="Method2_20261120",
    priors=pr,samplers=sampler_info,seed=seed,
    init_method=if(is.null(init)) "kmeans" else "user",
    package_version=as.character(utils::packageVersion("BPFC")))
  feature_ids <- rownames(y) %||% as.character(seq_len(n))
  names(out$clustering) <- names(out$cluster_confidence) <-
    names(out$cluster_uncertainty) <- feature_ids
  rownames(out$cluster_prob) <- feature_ids
  colnames(out$posterior_samples$z) <- feature_ids
  out$data_info <- dat[setdiff(names(dat), "y")]
  out$data_info$feature_ids <- feature_ids
  out$data_info$col_order <- "interleaved: [state1_t1,state2_t1,...,state1_tT,state2_tT]"
  out$diagnostics <- .bpfc_single_diagnostics(out$posterior_samples$all_params)
  class(out) <- "mcmcgraph_result"
  out
}

.bpfc_priors <- function(priors, J, TT) {
  required <- c("alpha","beta_mean","beta_cov","a_v","b_v","a_u","b_u",
                "mu_phi","var_phi","mu_psi","var_psi","a_lambda","b_lambda")
  if (!is.list(priors) || !length(priors) || is.null(names(priors)) ||
      any(!nzchar(names(priors))) || anyDuplicated(names(priors)))
    stop("Supply a uniquely named, complete priors list; Method2 gives no numerical defaults.")
  unknown <- setdiff(names(priors),required)
  if (length(unknown)) stop("Unknown priors: ",paste(unknown,collapse=", "))
  missing <- setdiff(required,names(priors))
  if (length(missing)) stop("Method2 requires explicit hyperparameters: ",paste(missing,collapse=", "))
  p <- priors
  expand <- function(v,len,nm,positive=TRUE) {
    if (!is.numeric(v) || !length(v) || !length(v) %in% c(1L,len) ||
        any(!is.finite(v)) || (positive && any(v <= 0))) stop("Invalid prior ",nm)
    rep(v,length.out=len)
  }
  p$alpha <- expand(p$alpha,J,"alpha")
  for (nm in c("a_v","b_v","a_u","b_u")) p[[nm]] <- expand(p[[nm]],TT,nm)
  for (nm in c("mu_phi","var_phi","mu_psi","var_psi"))
    p[[nm]] <- rep(expand(p[[nm]],1L,nm,!startsWith(nm,"mu_")),2L)
  for (nm in c("a_lambda","b_lambda")) p[[nm]] <- expand(p[[nm]],2L,nm)
  if (is.numeric(p$beta_mean) && is.null(dim(p$beta_mean)) && length(p$beta_mean)==10L)
    p$beta_mean <- matrix(rep(p$beta_mean,each=J),J,10L)
  if (!is.numeric(p$beta_mean) || !identical(dim(p$beta_mean),c(J,10L)) ||
      any(!is.finite(p$beta_mean))) stop("beta_mean must be length 10 or J x 10.")
  if (is.matrix(p$beta_cov) && identical(dim(p$beta_cov),c(10L,10L)))
    p$beta_cov <- aperm(array(rep(p$beta_cov,J),c(10L,10L,J)),c(3,1,2))
  if (!is.numeric(p$beta_cov) || !identical(dim(p$beta_cov),c(J,10L,10L)) ||
      any(!is.finite(p$beta_cov))) stop("beta_cov must be 10 x 10 or J x 10 x 10.")
  for(j in seq_len(J)) {
    b <- p$beta_cov[j,,]
    if (!isSymmetric(b) || inherits(try(chol(b),silent=TRUE),"try-error"))
      stop("beta_cov must be symmetric positive definite for each cluster.")
  }
  p
}

.bpfc_init <- function(y,J,Z0,times,init) {
  if (!is.null(init) && !is.list(init)) stop("init must be a list.")
  if (!is.null(init) && any(!names(init) %in% c("z","beta","p","v_sq","u_sq","phi","psi","lambda")))
    stop("Unknown initialization field.")
  n <- nrow(y); TT <- length(times)
  z <- init$z
  if (is.null(z)) {
    if (J==1L) z <- rep(1L,n) else {
      scaled <- scale(y); scaled[!is.finite(scaled)] <- 0
      if (nrow(unique(scaled)) < J) stop("Too few distinct profiles for kmeans; provide init$z.")
      z <- stats::kmeans(scaled,centers=J,nstart=10,iter.max=100)$cluster
    }
  }
  if (length(z)!=n || any(!is.finite(z)) || any(z != floor(z)) || any(z<1|z>J))
    stop("init$z requires n integer labels in 1:J.")
  beta <- matrix(0,J,10L)
  global <- as.numeric(solve(crossprod(Z0)+diag(1e-6,10),crossprod(Z0,colMeans(y))))
  for(j in seq_len(J)) {
    ii <- which(z==j)
    beta[j,] <- if(length(ii)) as.numeric(solve(crossprod(Z0)+diag(1e-6,10),
      crossprod(Z0,colMeans(y[ii,,drop=FALSE])))) else global
  }
  residual <- y-beta[z,,drop=FALSE] %*% t(Z0)
  variances <- pmax(colMeans(residual^2),1e-4)
  p <- tabulate(z,nbins=J)+0.1
  ini <- list(z=as.integer(z), beta=beta, p=p/sum(p),
    v_sq=variances[seq.int(1L,2L*TT,2L)],u_sq=variances[seq.int(2L,2L*TT,2L)],
    phi=c(0,0),psi=c(0,0),lambda=c(1,1)/max(times))
  if (!is.null(init)) ini <- utils::modifyList(ini,init)
  if (!is.matrix(ini$beta) || !identical(dim(ini$beta),c(J,10L)) ||
      any(!is.finite(ini$beta))) stop("init$beta requires a finite J x 10 matrix.")
  if (length(ini$p)!=J || any(!is.finite(ini$p)) || any(ini$p<=0) || abs(sum(ini$p)-1)>1e-8)
    stop("init$p requires positive weights summing to one.")
  sad_covariance(times,ini$phi,ini$psi,ini$lambda,ini$v_sq,ini$u_sq)
  ini$z <- as.integer(ini$z)
  ini
}

.bpfc_summarize <- function(post,J,P,Z0,n) {
  S <- nrow(post)
  extract <- function(nm) {
    cols <- grep(paste0("^",nm,"\\["),colnames(post),value=TRUE)
    idx <- as.integer(sub(".*\\[([0-9]+)\\]$","\\1",cols))
    post[,cols[order(idx)],drop=FALSE]
  }
  z <- if (J>1L) extract("z") else matrix(1L,S,n)
  storage.mode(z) <- "integer"
  p <- if (J>1L) extract("p") else matrix(1,S,1L)
  beta <- array(0,c(S,J,P))
  bcols <- grep("^beta\\[",colnames(post),value=TRUE)
  bix <- regmatches(bcols,regexec("^beta\\[([0-9]+), *([0-9]+)\\]$",bcols))
  for(k in seq_along(bcols)) beta[,as.integer(bix[[k]][2]),as.integer(bix[[k]][3])] <- post[,bcols[k]]
  rl <- .relabel_ecr(z,J)
  z <- rl$z; p <- .apply_perm_vec(p,rl$perm); beta <- .apply_perm_beta(beta,rl$perm)
  aligned <- post[,!grepl("^z\\[",colnames(post)),drop=FALSE]
  if (J==1L) aligned <- cbind(aligned, "p[1]"=rep(1,S))
  for(j in seq_len(J)) aligned[,sprintf("p[%d]",j)] <- p[,j]
  for(k in seq_along(bcols)) aligned[,bcols[k]] <- beta[,as.integer(bix[[k]][2]),as.integer(bix[[k]][3])]
  prob <- t(matrix(vapply(seq_len(ncol(z)),function(i) tabulate(z[,i],nbins=J)/S,
                          numeric(J)),nrow=J,ncol=ncol(z)))
  colnames(prob) <- paste0("cluster",seq_len(J))
  clustering <- max.col(prob,ties.method="first")
  confidence <- apply(prob,1L,max)
  beta_mean <- apply(beta,c(2,3),mean)
  draws <- lapply(c("phi","psi","lambda","v_sq","u_sq"),extract)
  names(draws) <- c("phi","psi","lambda","v_sq","u_sq")
  pm <- lapply(draws,colMeans)
  pm <- c(list(beta=beta_mean,p=colMeans(p),mu=beta_mean %*% t(Z0)),pm)
  list(clustering=clustering,cluster_counts=table(factor(clustering,levels=seq_len(J))),
    cluster_prob=prob,cluster_confidence=confidence,cluster_uncertainty=1-confidence,
    posterior_mean=pm,posterior_samples=c(list(z=z,p=p,beta=beta,all_params=aligned),draws),
    relabeling=list(method="iterative ECR",pivot=rl$pivot,permutation=rl$perm))
}
