# Reproduce suite

Scripts that regenerate the paper's Task-1 simulation study (stability,
ablations, convergence) and the real-data workflow.

## Layout

| File | Role |
|------|------|
| `config.R` | shared paths, scale knobs (`MCG_R`, `MCG_NITER`, ...), per-replicate seeds |
| `lib_metrics.R` | ACC / ARI / NMI / Macro-F1 / Balanced-ACC / small-cluster recall (LSAP-aligned) |
| `lib_sim.R` | scenario + mirror generators, same-model EM, generalised SAD-MCMC kernel, baselines |
| `01_run_simulations.R` | runs 1A / 1B / 1C / 1D, one checkpointed `.rds` per cell |
| `02_make_tables.R` | mean±SD + median[IQR] tables and paired Wilcoxon tests |
| `03_make_figures.R` | boxplots, init-sensitivity, ablation, mirror mechanism, traces |
| `04_real_data_analysis.R` | fixed-J/BIC cross-check → clustering → mean curves → diagnostics → module network; includes the manuscript mouse data |
| `05_overfitted_selection.R` | one sparse `J_max` MCMC fit → occupied-count posterior → sorted component mass → working resolution |
| `run_full.sh` | bounded-parallel, resumable launcher for the whole study |

## Experiments

- **1A** — multi-replicate stability: 3 scenarios (K=3/5/8) × R replicates × 5
  methods (MCMC, EM, LOP-k-means, SVM-RBF, XGBoost), EM run from multiple random
  starts; an init-sensitivity cell runs EM and MCMC from identical diverse starts.
- **1B** — clustering-unit ablation: **joint** vs **concatenation** vs
  **independent**, all sharing the *same* Bayesian SAD-MCMC kernel (only the
  block structure + basis change). Plus two mirror scenarios:
  - *pairing* mirror (Cluster 5/7 motif) — breaks the naive independent pipeline;
  - *orthogonal* mirror — breaks concatenation (cluster differences are
    orthogonal to the single concatenated-curve basis).
- **1C** — secondary ablations: covariance form (SAD vs AR(1) vs independence vs
  unstructured GMM), LOP order (3/4/5), prior sensitivity (InverseGamma `(a,b)`,
  `sigma_beta`).
- **1D** — convergence: ≥3 chains from dispersed starts, Gelman-Rubin R-hat and
  effective sample size, trace plots.

## Running

```bash
export MCG_ROOT=/path/to/project        # contains benchmark/data/sim_truth_K*.rds
bash BPFC/inst/reproduce/run_full.sh 3      # 3 parallel workers
Rscript BPFC/inst/reproduce/02_make_tables.R
Rscript BPFC/inst/reproduce/03_make_figures.R

# Manuscript mouse fit: bundled 11833 x 32 matrix, J=15, 30000 iterations,
# separate State-1 phi and State-2 psi (run inside tmux/nohup on a server).
MCG_DATASET=mouse Rscript BPFC/inst/reproduce/04_real_data_analysis.R

# Manuscript one-run overfitted mixture: J_max=25, e0=.02, 90% mass rule,
# 30000 iterations. Run inside tmux/nohup on a server.
MCG_DATASET=mouse Rscript BPFC/inst/reproduce/05_overfitted_selection.R
```

Scale is set inside `run_full.sh` (and overridable via `MCG_R` / `MCG_NITER`).
The committed defaults are tuned to finish overnight on a 16 GB machine; for the
as-published R = 50 / niter = 8000 setting, raise those literals on a larger
host. All cells are seeded (`mcg_seed()`) and checkpointed, so runs are
deterministic and resumable. Outputs land in `results/reproduce/`.

The analysis-ready mouse matrix is bundled in `inst/extdata/mouse/`, together
with provenance, dimensions, and SHA-256 checksums. The fixed-J mouse command
writes the fitted object, BIC, cluster assignments, posterior mean curves,
parameter diagnostics, and module network to `results/reproduce/real/`. The
overfitted command writes the complete `J_max` fit, occupied-count posterior,
sorted component spectrum, selected dominant clustering, and three-panel
selection diagnostic to `results/reproduce/overfitted/`.
