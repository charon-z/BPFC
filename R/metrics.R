#' Evaluate clustering recovery as in Method2 Section 2.8
#'
#' @param truth,prediction Nonmissing true and predicted labels of equal length.
#' @return Named vector ACC, ARI, NMI, Macro_F1, Balanced_ACC,
#'   SmallCluster_Recall. Matching maximizes ACC using a padded one-to-one
#'   assignment. Unmatched predicted clusters remain unmatched, avoiding label
#'   collisions. Macro-F1 averages over TRUE clusters; NMI uses arithmetic
#'   entropy normalization 2*MI/(Htruth+Hprediction). These conventions are
#'   explicit because Method2 does not specify them. Small-cluster recall
#'   averages true clusters with proportion strictly below 0.10, and is NA
#'   when none exists. ARI and NMI use original partitions.
#' @export
clustering_metrics <- function(truth,prediction) {
  if (!length(truth) || length(truth)!=length(prediction) || anyNA(truth) || anyNA(prediction))
    stop("Provide equal-length nonmissing clustering labels.")
  a <- as.integer(factor(truth)); b <- as.integer(factor(prediction)); n <- length(a)
  ka <- max(a); kb <- max(b); K <- max(ka,kb)
  tab <- table(factor(a,levels=seq_len(ka)),factor(b,levels=seq_len(kb)))
  padded <- matrix(0,K,K); padded[seq_len(ka),seq_len(kb)] <- tab
  # Rows are predicted labels; columns are true labels (plus dummy labels).
  mapping <- as.integer(clue::solve_LSAP(t(padded),maximum=TRUE))
  aligned <- mapping[b]
  tp <- vapply(seq_len(ka),function(k) sum(a==k & aligned==k),numeric(1))
  fp <- vapply(seq_len(ka),function(k) sum(a!=k & aligned==k),numeric(1))
  size <- tabulate(a,nbins=ka); recall <- tp/size
  choose2 <- function(v) v*(v-1)/2
  N <- choose2(n); aa <- sum(choose2(rowSums(tab))); bb <- sum(choose2(colSums(tab)))
  cc <- sum(choose2(tab)); expected <- if(N>0) aa*bb/N else 0
  den <- (aa+bb)/2-expected
  ari <- if(abs(den)<1e-12) 1 else (cc-expected)/den
  pij <- tab/n; pi <- rowSums(pij); pj <- colSums(pij)
  independent <- outer(pi,pj); nz <- pij>0
  mi <- sum(pij[nz]*log(pij[nz]/independent[nz]))
  entropy <- -sum(pi*log(pi))-sum(pj*log(pj))
  nmi <- if(entropy<1e-12) 1 else 2*mi/entropy
  small <- size/n < 0.10
  c(ACC=sum(tp)/n,ARI=ari,NMI=nmi,Macro_F1=mean(2*tp/(size+tp+fp)),
    Balanced_ACC=mean(recall),SmallCluster_Recall=if(any(small)) mean(recall[small]) else NA_real_)
}
