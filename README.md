# BPFC

<!-- badges: start -->
[![R-CMD-check](https://github.com/charon-z/BPFC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/charon-z/BPFC/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE.md)
<!-- badges: end -->

Bayesian MCMC functional clustering for paired (K = 2) longitudinal data with a
structured antedependence (SAD) covariance, built on the `nimble` backend.

Each feature (e.g. a SNP effect) carries the longitudinal trajectories of **two
individuals**. Instead of clustering each individual separately, BPFC
clusters the two trajectories **jointly**:

- a single cluster label `z_i` is shared by both individuals of a feature;
- the cluster mean is a **paired** Legendre-orthogonal-polynomial (order 4)
  curve `mu_j = X beta_j` spanning both individuals, so a cluster can express
  synergistic or antagonistic joint dynamics;
- within-state residual correlation is modelled by a block-diagonal SAD(1)
  covariance with separate correlation-decay parameters `phi` (State 1) and
  `psi` (State 2), each truncated to `(-1, 1)`, and time-point-specific
  innovation variances `v_sq`.

Inference uses a hybrid MCMC scheme combining categorical/conjugate updates,
slice sampling, and adaptive block random-walk Metropolis updates via `NIMBLE`.
Label switching is resolved with an ECR-style relabeling step before computing
MAP cluster labels and posterior cluster probabilities. A sparse overfitted
finite mixture can infer a working module resolution from one MCMC fit;
repeated fixed-J fits and BIC remain available as a supplementary cross-check.

## Installation

```r
install.packages("remotes")
remotes::install_github("charon-z/BPFC")
```

### System requirements

- R >= 3.5.
- **A working C++ toolchain.** `nimble` compiles each model to C++, so a
  compiler must be available *before* installation. This is the most common
  point of failure.
  - macOS: install the Xcode command-line tools (`xcode-select --install`) and,
    for the Fortran/OpenMP bits, the macOS R toolchain from
    <https://mac.r-project.org/tools/>.
  - Linux: `r-base-dev` plus `gcc`/`gfortran`.
  - Windows: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
    matching your R version.
  - Verify with `nimble::nimbleOptions()` and by compiling a trivial model; see
    the [nimble installation guide](https://r-nimble.org/download).
- Imported R packages: `nimble`, plus base `stats`, `graphics`, `utils`. The
  reproduce scripts additionally use `mclust`, `clue`, `e1071`, `xgboost`,
  `coda`, `ggplot2` and `igraph`.

## Quick start

```r
library(BPFC)

# Example data: a list with $y (n x 2d wide matrix), $times, $d_single, $truth
data(example_binary)

fit <- run_mcmc_binary(
  example_binary$y,
  J      = 3,
  times  = example_binary$times,
  niter  = 3000,
  thin   = 1
)

# The paper model (default) estimates phi and psi separately.
fit$posterior_mean[c("phi", "psi")]

# Optional restricted special case, if scientifically justified:
# fit_shared <- run_mcmc_binary(..., two_phi = FALSE)

fit$clustering          # MAP cluster label per feature
fit$cluster_prob        # posterior cluster probabilities (n x J)
fit$cluster_uncertainty # 1 - max posterior probability per feature

# Trace / density diagnostics
plot_trace_density(fit, params = c("^phi", "^psi", "^p\\["))

# Readable mixing-proportion traces for all components
plot_mixing_trace(fit)
```

## Choosing the number of clusters in one run

`run_overfitted_binary()` makes the manuscript's overfitted-mixture procedure
explicit. It fits exactly one model at `J_max`, computes the posterior of the
occupied-component count, sorts the posterior-mean component weights, and
selects the smallest number reaching the stated cumulative-mass threshold.
The default 90% rule is explicit and can be changed; no value of K is hard-coded.

```r
overfit <- run_overfitted_binary(
  example_binary$y,
  J_max = 6,
  e0 = 0.02,
  mass_threshold = 0.90,
  times = example_binary$times,
  niter = 3000
)

overfit$selected_K        # working resolution from sorted posterior mass
overfit$K_occ_mode        # posterior mode of occupied components
overfit$occupied_posterior
overfit$component_spectrum
plot_overfitted_selection(overfit)
```

The full `J_max` fit is retained as `overfit$fit`. The top-level clustering is
computed from that same run after restricting posterior assignment
probabilities to the selected dominant components.

### Supplementary fixed-J BIC cross-check

`fit_many_J()` is deliberately separate: it refits each candidate J and must
not be described as the one-run overfitted-mixture procedure.

```r
res <- fit_many_J(example_binary$y, J_grid = 2:6,
                  times = example_binary$times, niter = 3000)
res$eval            # data.frame: J, loglik, BIC
plot_bic(res$eval)  # BIC vs J
```

## Input formats

`run_mcmc_binary()` accepts either layout (auto-detected):

- **wide**: an `n x (2d)` matrix ordered as
  `[ind1_t1..ind1_td, ind2_t1..ind2_td]`;
- **long**: a data frame with columns `snp` (feature id), `ind` (exactly two
  individuals), `time`, `y`.

## Method summary

| Component | Specification |
|-----------|---------------|
| Functional basis | Legendre orthogonal polynomials, order 4 (q = 5; paired P = 10) |
| Mixture | finite mixture, `z ~ Categorical(p)`, `p ~ Dirichlet` |
| Cluster mean | `mu_j = X beta_j`, `beta_j ~ MVN` |
| Covariance | block-diagonal SAD(1); State 1 uses `phi`, State 2 uses `psi`; both are truncated-normal on `(-1, 1)`; `v_sq ~ InvGamma` |
| Sampler | `nimble`: categorical (`z`), Dirichlet-conjugate (`p`), slice (`phi`, `psi`, `v_sq`), RW-block (`beta`) |
| Post-processing | ECR relabeling, MAP labels, posterior cluster probabilities, BIC |

## Main exported functions

| Function | Purpose |
|----------|---------|
| `run_mcmc_binary()` | Fit the Bayesian SAD-MCMC functional clustering at a fixed J |
| `run_overfitted_binary()` | Infer a working module resolution from one sparse overfitted-mixture MCMC run |
| `plot_overfitted_selection()` | Plot occupied-count, sorted-weight and cumulative-mass diagnostics |
| `as_binary_data()` | Coerce wide (`n x 2d`) or long data into the required format |
| `eval_bic()` | Posterior-mean plug-in BIC for a fitted object |
| `fit_many_J()` | Sweep a grid of J and return fits + a BIC table |
| `plot_bic()` | Plot BIC versus J |
| `plot_mixing_trace()` | Plot mixing-proportion traces with explicit axis labels and mathematical component notation |
| `plot_trace_density()` | Trace + density diagnostics for selected parameters |
| `as_eval_df()` | Normalise evaluation results into a `J`/`BIC` data frame |

## Reproducing the paper simulations and figures

Scripts live in [`inst/reproduce/`](inst/reproduce). They cover the
multi-replicate stability study, the joint/concatenation/independent ablation
(plus the two mirror scenarios), the secondary ablations (covariance form, LOP
order, priors) and the convergence diagnostics.

```bash
# from the project root (contains benchmark/data/sim_truth_K*.rds)
export MCG_ROOT=$(pwd)
# full study, 6 workers (checkpointed + resumable; many hours):
bash BPFC/inst/reproduce/run_full.sh 6
# or a quick validation pass:
MCG_R=5 MCG_NITER=3000 bash BPFC/inst/reproduce/run_full.sh 6
# then build tables and figures:
Rscript BPFC/inst/reproduce/02_make_tables.R
Rscript BPFC/inst/reproduce/03_make_figures.R
# fixed-J/BIC real-data workflow (supplementary cross-check):
Rscript BPFC/inst/reproduce/04_real_data_analysis.R   # bundled toy example
MCG_DATASET=mouse Rscript BPFC/inst/reproduce/04_real_data_analysis.R
# one-run overfitted-mixture selection (J_max=25, e0=.02 for mouse):
MCG_DATASET=mouse Rscript BPFC/inst/reproduce/05_overfitted_selection.R
```

Simulation data are regenerated deterministically from the stored scenario
parameters and fixed seeds (`mcg_seed()`), so no large data files are needed to
reproduce the simulation tables.

## Computational cost

Measured on an Apple-silicon 10-core machine (16 GB RAM), single thread:

| Setting | Time |
|---------|------|
| n = 3000, 32 dims, J = 3, niter = 2000 | ~90 s |
| n = 3000, 32 dims, J = 3, niter = 8000 | ~6 min |
| n = 3000, 32 dims, J = 8, niter = 8000 | ~10 min |

Cost scales roughly linearly in both `n` (the per-feature categorical update
dominates) and `niter`, and grows with J through the number of cluster-mean
blocks and categorical levels. The EM comparator (no compilation) is well under
a minute even with 10 random starts. `run_full.sh` parallelises across
replicates; the full R = 50 study is on the order of 100 core-hours.

A resumable single-thread benchmark for BPFC-MCMC, BPFC-EM, Gaussian mixtures,
K-means and agglomerative clustering is provided in
[`inst/benchmark/`](inst/benchmark). It records end-to-end elapsed time and peak
resident memory in separate fresh R processes.

## License

MIT. See [LICENSE.md](LICENSE.md).
