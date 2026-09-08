make_fake_overfit_fit <- function() {
  z <- rbind(
    c(1, 1, 2, 2, 3),
    c(1, 1, 2, 2, 2),
    c(1, 1, 2, 2, 4),
    c(1, 1, 2, 2, 2)
  )
  p <- matrix(rep(c(0.60, 0.30, 0.08, 0.02), each = 4L), nrow = 4L)
  cluster_prob <- t(vapply(
    seq_len(ncol(z)),
    function(i) tabulate(z[, i], nbins = 4L) / nrow(z),
    numeric(4L)
  ))
  colnames(cluster_prob) <- paste0("cluster", 1:4)
  fit <- list(
    clustering = apply(z, 2, function(x) which.max(tabulate(x, nbins = 4L))),
    cluster_prob = cluster_prob,
    posterior_samples = list(z = z, p = p),
    model_info = list(J = 4L, e0 = 0.02)
  )
  class(fit) <- "mcmcgraph_result"
  fit
}

test_that("overfitted summary implements occupied-count and mass equations", {
  fit <- make_fake_overfit_fit()
  out <- BPFC:::.summarize_overfitted_fit(fit, mass_threshold = 0.90)

  expect_equal(out$occupied_count, c(3L, 2L, 3L, 2L))
  expect_equal(out$occupied_posterior$probability, c(0, 0.5, 0.5, 0))
  expect_identical(out$K_occ_mode, 2L)
  expect_equal(unname(out$K_occ_interval), c(2L, 3L))
  expect_identical(out$selected_K, 2L)
  expect_equal(out$component_spectrum$posterior_mean_weight,
               c(0.60, 0.30, 0.08, 0.02))
  expect_equal(out$component_spectrum$cumulative_mass,
               c(0.60, 0.90, 0.98, 1.00))
  expect_equal(ncol(out$cluster_prob), 2L)
  expect_true(all(out$clustering %in% 1:2))
})

test_that("mass threshold changes the working resolution without hard-coding K", {
  fit <- make_fake_overfit_fit()
  out <- BPFC:::.summarize_overfitted_fit(fit, mass_threshold = 0.91)
  expect_identical(out$selected_K, 3L)
})

test_that("overfitted interfaces validate selection settings", {
  expect_equal(formals(run_overfitted_binary)$e0, 0.02)
  expect_equal(formals(run_overfitted_binary)$mass_threshold, 0.90)
  expect_error(run_overfitted_binary(matrix(1, 4, 4), J_max = 1), "J_max")
  expect_error(run_overfitted_binary(matrix(1, 4, 4), J_max = 3, e0 = 0), "e0")
  expect_error(
    BPFC:::.summarize_overfitted_fit(make_fake_overfit_fit(), mass_threshold = 0),
    "mass_threshold"
  )
})

test_that("overfitted diagnostic plot accepts a summarized result", {
  summary_fields <- BPFC:::.summarize_overfitted_fit(make_fake_overfit_fit())
  result <- c(
    list(J_max = 4L, e0 = 0.02, mass_threshold = 0.90,
         credible_prob = 0.95),
    summary_fields
  )
  class(result) <- "bpfc_overfit_result"
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  plotted <- plot_overfitted_selection(result)
  expect_equal(plotted, result$component_spectrum)
})

test_that("fixed-J prior keeps e0 at one", {
  default_priors <- eval(formals(run_mcmc_binary)$priors)
  expect_identical(default_priors$e0, 1)
})
