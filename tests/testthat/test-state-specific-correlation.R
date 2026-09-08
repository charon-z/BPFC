test_that("the manuscript covariance model is the default", {
  expect_identical(formals(run_mcmc_binary)$two_phi, TRUE)
})

test_that("binary initialization keeps state-specific correlation estimates", {
  data(example_binary)
  init <- BPFC:::initialization_kmeans_binary(
    example_binary$y, J = 3,
    Z0_binary = BPFC:::make_Z0_binary(example_binary$times)
  )

  expect_true(all(c("phi", "psi") %in% names(init)))
  expect_length(init$phi, 1L)
  expect_length(init$psi, 1L)
  expect_true(is.finite(init$phi))
  expect_true(is.finite(init$psi))
})

test_that("BIC uses psi and counts both state correlations", {
  data(example_binary)
  y <- example_binary$y[1:8, , drop = FALSE]
  J <- 1L
  P <- 10L
  base_fit <- list(
    posterior_mean = list(
      beta = matrix(0, J, P), p = 1,
      v_sq = rep(1, ncol(y)), phi = 0.25, psi = 0.25
    ),
    model_info = list(J = J, P = P, two_phi = FALSE),
    data_info = list(times_single = example_binary$times)
  )
  separate_fit <- base_fit
  separate_fit$model_info$two_phi <- TRUE

  shared <- eval_bic(base_fit, y, format = "wide")
  separate_equal <- eval_bic(separate_fit, y, format = "wide")
  expect_equal(separate_equal$loglik, shared$loglik, tolerance = 1e-12)
  expect_equal(separate_equal$BIC - shared$BIC, log(nrow(y)), tolerance = 1e-12)

  separate_fit$posterior_mean$psi <- -0.55
  separate_different <- eval_bic(separate_fit, y, format = "wide")
  expect_false(isTRUE(all.equal(separate_different$loglik, shared$loglik)))
})
