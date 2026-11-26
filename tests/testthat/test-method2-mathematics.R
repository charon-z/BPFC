# Dense reference uses the inverse innovation operator rather than powers of A.
reference_cov <- function(times,phi,psi,lambda,v,u) {
  T <- length(times); A <- rbind(c(phi[1],psi[1]),c(psi[2],phi[2]))
  H <- diag(2*T); O <- matrix(0,2*T,2*T)
  for(t in seq_len(T)) {
    ii <- c(2*t-1,2*t)
    if(t>1) H[ii,ii-2] <- -A
    r <- exp(-lambda[1]*times[t])-exp(-lambda[2]*times[t])
    O[ii,ii] <- matrix(c(v[t],r*sqrt(v[t]*u[t]),r*sqrt(v[t]*u[t]),u[t]),2)
  }
  L <- solve(H)
  L %*% O %*% t(L)
}
dense_log <- function(x,mu,R) {
  r <- x-mu; U <- chol(R)
  -length(r)/2*log(2*pi)-sum(log(diag(U)))-sum(backsolve(U,r,transpose=TRUE)^2)/2
}

test_that("Eq 8 and innovation density agree with an independent dense MVN", {
  set.seed(412)
  for(times in list(c(0,0.5,2,5,10),c(1,3,7,10,20))) {
    phi <- c(0.4,-0.3); psi <- c(0.23,-0.17); lambda <- c(0.07,0.4)
    v <- c(.5,1,1.8,.9,2); u <- rev(v)+.2
    ref <- reference_cov(times,phi,psi,lambda,v,u)
    got <- sad_covariance(times,phi,psi,lambda,v,u)
    expect_equal(got$R,ref,tolerance=1e-12)
    expect_true(all(eigen(got$R,symmetric=TRUE)$values>0))
    x <- rnorm(10); mu <- rnorm(10)
    ll <- BPFC:::dBivSAD(x,mu,phi,psi,lambda,v,u,times,log=1L)
    expect_equal(ll,dense_log(x,mu,ref),tolerance=1e-11)
    # Eq. 9, both directions and different observation times.
    A <- got$A
    pow <- function(n) { ans<-diag(2); if(n>0) for(k in seq_len(n)) ans<-ans%*%A;ans }
    for(t in 1:5) for(s in 1:5) {
      C <- matrix(0,2,2)
      for(h in seq_len(min(t,s))) C <- C+pow(t-h)%*%got$Omega[c(2*h-1,2*h),c(2*h-1,2*h)]%*%t(pow(s-h))
      expect_equal(got$R[2*t-1,2*s],C[1,2],tolerance=1e-12)
    }
  }
})

test_that("original times and unrestricted finite-horizon A follow Method2", {
  times <- c(0,2,6,12)
  cv <- sad_covariance(times,c(.8,.8),c(.6,.6),c(.1,.5),rep(1,4),rep(1,4))
  expect_gt(max(Mod(eigen(cv$A)$values)),1)
  expect_true(all(eigen(cv$R,symmetric=TRUE)$values>0))
  expect_equal(cv$rho[1],0)
  expect_equal(cv$rho,exp(-.1*times)-exp(-.5*times))
  expect_equal(sad_covariance(times,c(.8,.8),c(.6,.6),c(.5,.1),rep(1,4),rep(1,4))$rho,-cv$rho)
})

test_that("basis, data layouts and long-format pairing implement Eqs 1-4", {
  times <- c(0,1,4,7,10)
  concat <- matrix(1:20,2)
  wanted <- concat[,as.vector(rbind(1:5,6:10))]
  got <- as_binary_data(concat,times=times,layout="concat")
  expect_equal(got$y,wanted,ignore_attr=TRUE)
  expect_equal(as_binary_data(wanted,times=times)$y,got$y)
  B <- BPFC:::make_Z0_binary(times)
  expect_equal(B[1,1:5],c(1,-1,1,-1,1))
  expect_equal(B[10,6:10],rep(1,5))
  expect_equal(B[seq(1,10,2),6:10],matrix(0,5,5))
  long <- expand.grid(snp=c("a","b"),ind=c("control","treated"),time=times)
  long$y <- wanted[cbind(match(long$snp,c("a","b")),2*(match(long$time,times)-1)+match(long$ind,c("control","treated")))]
  expect_equal(as_binary_data(long,format="long")$y,wanted,ignore_attr=TRUE)
  expect_error(as_binary_data(rbind(long,long[1,])),"Duplicate")
  expect_error(as_binary_data(long[-1,]),"Complete")
  expect_error(as_binary_data(wanted),"original")
  expect_error(as_binary_data(wanted,times=c(-1,1,4,7,10)),"nonnegative")
})

test_that("beta priors retain full covariance and reject silent unknown options", {
  C <- diag(2,10); C[1,2] <- C[2,1] <- .3
  pr <- BPFC:::.bpfc_priors(utils::modifyList(test_priors(),list(beta_cov=C,alpha=c(.2,.4))),2L,5L)
  expect_equal(pr$beta_cov[2,,],C)
  expect_equal(pr$alpha,c(.2,.4))
  expect_error(BPFC:::.bpfc_priors(list(alpha=1),2L,5L),"explicit hyperparameters")
  expect_error(BPFC:::.bpfc_priors(c(test_priors(),list(two_phi=TRUE)),2L,5L),"Unknown")
  expect_error(BPFC:::.bpfc_priors(utils::modifyList(test_priors(),list(beta_cov=-diag(10))),2L,5L),"positive definite")
})

test_that("ECR aligns allocations, weights and means with one permutation", {
  z <- rbind(c(1,1,2,2),c(2,2,1,1))
  rl <- BPFC:::.relabel_ecr(z,2L)
  expect_equal(rl$z[1,],rl$z[2,])
  weights <- BPFC:::.apply_perm_vec(rbind(c(.7,.3),c(.3,.7)),rl$perm)
  expect_equal(weights[1,],weights[2,])
  b <- array(0,c(2,2,10));b[1,1,]<-1;b[1,2,]<-2;b[2,1,]<-2;b[2,2,]<-1
  aligned <- BPFC:::.apply_perm_beta(b,rl$perm)
  expect_equal(aligned[1,,],aligned[2,,])
})

test_that("metrics honor strict rare-cluster threshold and rectangular matching", {
  truth <- c(rep("a",85),rep("b",10),rep("c",5))
  pred <- c(rep("x",85),rep("y",10),rep("z",5))
  expect_equal(unname(clustering_metrics(truth,pred)),rep(1,6))
  expect_true(is.na(clustering_metrics(rep(1:2,each=10),rep(2:1,each=10))["SmallCluster_Recall"]))
  expect_true(is.na(clustering_metrics(c(rep(1,90),rep(2,10)),c(rep(1,90),rep(2,10)))["SmallCluster_Recall"]))
  # An extra predicted label cannot collide with a true label after matching.
  split <- c(rep(1,40),rep(2,45),rep(3,10),rep(4,5))
  vals <- clustering_metrics(truth,split)
  expect_equal(vals["ACC"],c(ACC=.60))
  expect_equal(vals["SmallCluster_Recall"],c(SmallCluster_Recall=1))
  expect_equal(clustering_metrics(rep(1,5),rep(2,5))[c("ARI","NMI")],c(ARI=1,NMI=1))
})

test_that("NIMBLE inverse-gamma scale matches the declared prior density", {
  a <- 3.2; b <- 2.1; x <- c(.3,1,3.7)
  expected <- a*log(b)-lgamma(a)-(a+1)*log(x)-b/x
  expect_equal(nimble::dinvgamma(x,shape=a,scale=b,log=TRUE),expected,tolerance=1e-12)
})
