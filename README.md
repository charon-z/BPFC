# BPFC 0.3.0

Bayesian Paired Functional Clustering fits a finite mixture model to features
measured as two paired longitudinal trajectories. Each feature has one shared
cluster label and two state-specific means represented by fourth-degree
Legendre polynomials. The bivariate first-order structured antedependence
(SAD(1)) covariance includes within-state and cross-state lags, time-specific
innovation variances, and time-varying correlation between simultaneous
innovations.

## Installation

```r
install.packages("remotes")
remotes::install_github("charon-z/BPFC")
```

A C++ toolchain is required for NIMBLE model compilation. The source package
and the example below were tested with R 4.6.0 and NIMBLE 1.4.2. Install the
optional `coda` package for chain diagnostics.

## Run a small example

The prior values below are those of the reference simulation settings. The
short chain demonstrates the API.

```r
library(BPFC)

times <- 0:5
beta <- rbind(
  c(-1, .5, 0, 0, 0,  1, -.5, 0, 0, 0),
  c( 1, -.5, 0, 0, 0, -1, .5, 0, 0, 0)
)
sim <- simulate_bpfc(
  n = 24, times = times, beta = beta, p = c(.5, .5),
  phi = c(.45, .25), psi = c(.20, -.15), lambda = c(.08, .30),
  v_sq = rep(.12, length(times)), u_sq = rep(.14, length(times)),
  seed = 621124
)

priors <- list(
  alpha = 1, beta_mean = rep(0, 10), beta_cov = diag(25, 10),
  a_v = 2, b_v = 1, a_u = 2, b_u = 1,
  mu_phi = 0, var_phi = 1, mu_psi = 0, var_psi = 1,
  a_lambda = 2, b_lambda = 10
)

fit <- run_mcmc_binary(
  sim$y, J = 2, times = times, priors = priors,
  niter = 64, thin = 2, seed = 7
)
head(fit$cluster_prob)
eval_bic(fit, sim$y)
```

For wide input, default columns alternate states at each time
(`state1_t1, state2_t1, ..., state1_tT, state2_tT`). Use
`layout = "concat"` if all State 1 columns precede all State 2 columns.
For long input, supply the columns `snp`, `ind`, `time`, and `y`, with
complete paired observations. Supply original nonnegative, increasing times
for wide input.

## Select a working cluster count

```r
auto_priors <- priors
auto_priors$alpha <- NULL
auto <- fit_bpfc(
  sim$y, selection = "automatic", J_max = 3, times = times,
  auto_control = list(e0 = .02, mass_threshold = .90),
  priors = auto_priors, niter = 64, thin = 2, seed = 7
)
auto$selected_K
auto$selection_details$occupied_posterior
```

The reported count is the smallest number of leading components reaching the
specified cumulative posterior weight. Allocation probabilities retain the
full set of fitted components. For a separate comparison of fixed-`J`
models, use `fit_bpfc(selection = "bic", J_grid = ...)`.

## Reproduce the article analyses

The package implements the full bivariate model. The companion reproduction
bundle sets the initialization and proposal controls used in the formal
simulations. It also provides the reduced constant-correlation mouse model
(M1) and its analysis scripts.

Run `citation("BPFC")` for the software citation. The package is available
under the MIT license.
