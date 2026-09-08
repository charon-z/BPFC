# BPFC 0.1.3

- Added `run_overfitted_binary()`, which performs one sparse finite-mixture fit
  at `J_max` and implements the manuscript's occupied-component posterior,
  sorted posterior-mean weights, cumulative-mass resolution, and dominant
  clustering without a fixed-J grid search.
- Added `plot_overfitted_selection()` and a full mouse/example reproduction
  driver for the one-run selection analysis.
- Made the symmetric Dirichlet concentration `e0` configurable in
  `run_mcmc_binary()`; fixed-J fits retain `e0 = 1`, while the overfitted helper
  defaults to the manuscript sparse setting `e0 = 0.02`.
- Clarified throughout the documentation that `fit_many_J()` is repeated
  fixed-J fitting for a supplementary BIC cross-check, not the single-run
  overfitted-mixture method.

# BPFC 0.1.2

- Restored the manuscript model as the default: State 1 and State 2 now use
  separate SAD(1) correlation-decay parameters `phi` and `psi`. The previous
  shared-parameter model remains available explicitly with `two_phi = FALSE`.
- Updated initialization and plug-in BIC evaluation for the two-parameter
  covariance model.
- Added the analysis-ready Gray/Gough mouse SNP-effect matrix and a dedicated
  `MCG_DATASET=mouse` reproduction path using J = 15 and 30,000 iterations.

# BPFC 0.1.1

- Added `plot_mixing_trace()` for publication-ready mixing-proportion MCMC
  diagnostics with explicit axis labels, mathematical component notation,
  running means and posterior means.
- Added a resumable computational-performance benchmark covering BPFC-MCMC,
  BPFC-EM, Gaussian mixtures, K-means and agglomerative clustering, with
  end-to-end timing and peak-memory recording.
- Fixed benchmark summarization so methods without an MCMC iteration count are
  retained in the aggregated performance table.
- Declared testthat edition 3.

# BPFC 0.1.0

First public release accompanying the *Briefings in Bioinformatics*
manuscript.

## Features

* `run_mcmc_binary()` fits Bayesian finite-mixture functional clustering of
  paired (K = 2) longitudinal trajectories at a fixed number of clusters `J`.
* Paired cluster means use Legendre orthogonal polynomial bases (order 4).
* Within-state residual correlation uses a block-diagonal SAD(1) covariance
  with separate `phi` and `psi` parameters and time-specific innovation
  variances.
* `fit_many_J()` / `eval_bic()` / `plot_bic()` provide BIC-based model selection
  over a grid of `J`.
* `plot_trace_density()` provides MCMC trace and density diagnostics.
* `as_binary_data()` accepts both wide (`n x 2d`) and long input layouts.
* Label switching is resolved with an ECR-style relabeling step; the fit returns
  MAP labels, posterior cluster-membership probabilities, and per-feature
  assignment uncertainty.

## Reproducibility

* `inst/reproduce/` contains the simulation suite (stability, the
  joint/concatenation/independent ablation with two mirror scenarios, secondary
  ablations, and convergence diagnostics) plus real-data workflows, including
  the Gray/Gough mouse analysis used in the manuscript.
* All simulation inputs are regenerated deterministically from stored scenario
  parameters and fixed seeds.
