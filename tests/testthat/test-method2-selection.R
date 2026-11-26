test_that("overfit mass selection preserves all Eq 20 probabilities", {
  z <- rbind(c(1,1,2,3),c(1,1,2,3),c(1,2,2,3),c(1,1,2,3))
  p <- matrix(rep(c(.65,.30,.05),each=4),4,3)
  prob <- t(vapply(1:4,function(i) tabulate(z[,i],nbins=3)/4,numeric(3)))
  fit <- list(clustering=max.col(prob),cluster_prob=prob,
    posterior_samples=list(z=z,p=p),model_info=list(J=3))
  class(fit) <- "mcmcgraph_result"
  a <- BPFC:::.summarize_overfitted_fit(fit,mass_threshold=.9)
  expect_equal(a$selected_K,2L)
  expect_equal(a$occupied_count,rep(3L,4))
  expect_equal(a$cluster_prob,prob)
  expect_equal(a$clustering[4],3L)
  expect_equal(a$retained_probability[4],0)
  expect_equal(a$cluster_confidence,apply(prob,1,max))
  expect_equal(a$component_spectrum$cumulative_mass,c(.65,.95,1))
})
