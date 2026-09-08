# Mouse data used by the manuscript workflow

`gray_mouse_snp_effects_11833x32.tsv.gz` is the whitespace-delimited,
analysis-ready matrix used as
input to the BPFC mouse analysis. It contains 11,833 SNP-effect trajectories in
rows and 32 measurements in columns: 16 weekly values for State 1 (`T1`-`T16`)
followed by 16 weekly values for State 2 (`V1`-`V16`). Row order is preserved
from the original analysis; the source file did not contain SNP identifiers, so
the reproduction script assigns stable row identifiers `SNP_00001` through
`SNP_11833`.

The underlying mouse phenotypes and genotypes are available from the Jackson
Laboratory Mouse Phenome Database, project 539, as cited in the manuscript:

<https://phenome.jax.org/db/q?reqprojid=539&rtn=projects/projdet>

The bundled matrix is a derived, analysis-ready input rather than a replacement
for the primary MPD records. It is included so that reviewers can reproduce the
BPFC fit without relying on a changing external download interface.

Integrity checks:

- uncompressed matrix SHA-256:
  `ced218f31577a30df426cc7a646e11dc8bcafdab6f57d8ab79fd8afef9c7960e`
- compressed file SHA-256:
  `081211bf30765901ed36cf295c97a883a25ef09c958afb38727432106fa4060b`

From the repository root, run the manuscript configuration with:

```bash
MCG_DATASET=mouse Rscript inst/reproduce/04_real_data_analysis.R
```

The mouse defaults are fixed at `J = 15`, 30,000 MCMC iterations, a 25% burn-in,
and the two-parameter covariance model (`phi` for State 1 and `psi` for State
2). Override `MCG_JGRID` or `MCG_NITER` only for sensitivity or quick checks.
