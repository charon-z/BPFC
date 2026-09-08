#!/usr/bin/env Rscript
# 04_real_data_analysis.R -------------------------------------------------
# End-to-end real-data workflow with the published BPFC functions:
#   1. load a paired longitudinal matrix (n x 2d)
#   2. supplementary BIC cross-check over fixed-J fits (fit_many_J / eval_bic)
#   3. final clustering at the selected J
#   4. cluster mean-curve table + posterior cluster probabilities
#   5. trace/density diagnostics for key parameters
#   6. a module-level network graph (clusters linked by mean-curve similarity)
#
# Dataset selection (environment variable MCG_DATASET):
#   "example" (default) : the bundled example_binary toy data (runs in minutes)
#   "mouse"             : manuscript mouse SNP-effect trajectories (11833 x 32)
#   "lincs"             : LINCS L1000 vorinostat MCF7/PC3 (978 x 12), if present
#   "smillie"           : Smillie UC pseudobulk (if present)
# The analysis-ready mouse matrix is bundled under inst/extdata/mouse. Its
# provenance and checksum are documented there. The mouse defaults reproduce
# the manuscript fixed-J fit (J = 15, 30000 iterations, separate phi/psi).
# See 05_overfitted_selection.R for the distinct one-run J_max procedure.

suppressPackageStartupMessages(library(BPFC))
src_dir <- if (nzchar(Sys.getenv("MCG_LIB_DIR"))) Sys.getenv("MCG_LIB_DIR") else
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(src_dir) || !nzchar(src_dir)) src_dir <- "BPFC/inst/reproduce"
source(file.path(src_dir, "config.R"))
OUT <- file.path(MCG_OUT_DIR, "real"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dataset <- tolower(Sys.getenv("MCG_DATASET", unset = "example"))
J_grid_text <- Sys.getenv("MCG_JGRID", unset = "")
if (!nzchar(J_grid_text)) J_grid_text <- if (dataset == "mouse") "15L" else "2:6"
J_grid <- eval(parse(text = J_grid_text))
niter_text <- Sys.getenv("MCG_NITER", unset = "")
if (!nzchar(niter_text)) niter_text <- if (dataset == "mouse") "30000" else "8000"
niter <- as.integer(niter_text)

load_dataset <- function(name) {
  if (name == "mouse") {
    f <- system.file("extdata", "mouse", "gray_mouse_snp_effects_11833x32.tsv.gz",
                     package = "BPFC")
    if (!nzchar(f)) {
      f <- normalizePath(file.path(src_dir, "..", "extdata", "mouse",
                                   "gray_mouse_snp_effects_11833x32.tsv.gz"),
                         mustWork = FALSE)
    }
    if (!file.exists(f)) stop("Bundled mouse SNP-effect matrix not found: ", f)
    Y <- as.matrix(read.table(gzfile(f), header = TRUE, check.names = FALSE))
    storage.mode(Y) <- "double"
    rownames(Y) <- sprintf("SNP_%05d", seq_len(nrow(Y)))
    return(list(Y = Y, times = 1:16,
                label = "Gray/Gough mouse SNP-effect trajectories"))
  }
  if (name == "lincs") {
    f <- file.path(MCG_ROOT, "data/LINCS_L1000_vorinostat/processed",
                   "LINCS_Vorinostat_MCF7_PC3_978x12_expression.tsv")
    if (!file.exists(f)) stop("LINCS file not found: ", f)
    df <- read.delim(f, check.names = FALSE, row.names = 1)
    df$gene_symbol <- NULL
    return(list(Y = as.matrix(df), times = 1:6, label = "LINCS vorinostat (MCF7|PC3)"))
  }
  if (name == "smillie") {
    f <- file.path(MCG_ROOT, "data/SmillieUC/processed", "SmillieUC_5000x20_expression.tsv")
    if (!file.exists(f)) stop("Smillie file not found: ", f)
    Y <- as.matrix(read.delim(f, check.names = FALSE, row.names = 1))
    return(list(Y = Y, times = seq_len(ncol(Y) / 2), label = "Smillie UC"))
  }
  data("example_binary", package = "BPFC", envir = environment())
  list(Y = example_binary$y, times = example_binary$times, label = "example_binary (toy)")
}

ds <- load_dataset(dataset)
message(sprintf("[04] dataset=%s  n=%d  d=%d  J=%s  niter=%d",
                ds$label, nrow(ds$Y), ncol(ds$Y),
                paste(J_grid, collapse = ","), niter))

## 1-2. BIC sweep over J --------------------------------------------------
sweep <- BPFC::fit_many_J(ds$Y, J_grid = J_grid, times = ds$times,
                          niter = niter, two_phi = TRUE)
ev <- sweep$eval
write.csv(ev, file.path(OUT, sprintf("%s_BIC.csv", dataset)), row.names = FALSE)
bestJ <- ev$J[which.min(ev$BIC)]
if (length(J_grid) > 1L) {
  message(sprintf("[04] BIC-selected J = %d", bestJ))
} else {
  message(sprintf("[04] fitted manuscript J = %d", bestJ))
}

pdf(file.path(OUT, sprintf("%s_BIC.pdf", dataset)), width = 6, height = 4.5)
BPFC::plot_bic(ev); dev.off()

## 3. final clustering at best J ------------------------------------------
fit <- sweep$fits[[paste0("J", bestJ)]]
if (is.null(fit)) {
  fit <- BPFC::run_mcmc_binary(ds$Y, J = bestJ, times = ds$times,
                               niter = niter, two_phi = TRUE)
}
saveRDS(fit, file.path(OUT, sprintf("%s_J%d_fit.rds", dataset, bestJ)))

samples <- as.matrix(fit$posterior_samples$all_params)
posterior_summary <- data.frame(
  parameter = colnames(samples),
  mean = colMeans(samples),
  sd = apply(samples, 2, stats::sd),
  q2.5 = apply(samples, 2, stats::quantile, probs = 0.025),
  q50 = apply(samples, 2, stats::median),
  q97.5 = apply(samples, 2, stats::quantile, probs = 0.975),
  row.names = NULL,
  check.names = FALSE
)
write.csv(posterior_summary,
          file.path(OUT, sprintf("%s_J%d_posterior_summary.csv", dataset, bestJ)),
          row.names = FALSE)

feature_ids <- rownames(ds$Y)
if (is.null(feature_ids)) feature_ids <- sprintf("feature_%05d", seq_len(nrow(ds$Y)))
clusters <- data.frame(feature = feature_ids, cluster = fit$clustering,
                       fit$cluster_prob, uncertainty = fit$cluster_uncertainty)
write.csv(clusters, file.path(OUT, sprintf("%s_J%d_clusters.csv", dataset, bestJ)), row.names = FALSE)

## 4. cluster mean curves -------------------------------------------------
d_single <- fit$model_info$d_single
Z0 <- BPFC:::make_Z0_binary(ds$times)
mu <- fit$posterior_mean$beta %*% t(Z0)        # J x 2d
write.csv(data.frame(cluster = seq_len(bestJ), mu),
          file.path(OUT, sprintf("%s_J%d_mean_curves.csv", dataset, bestJ)), row.names = FALSE)

pdf(file.path(OUT, sprintf("%s_J%d_mean_curves.pdf", dataset, bestJ)),
    width = 10, height = max(4, 2.2 * ceiling(bestJ / 3)))
op <- par(mfrow = c(ceiling(bestJ / 3), 3), mar = c(2.4, 2.8, 2.1, 0.8),
          oma = c(2.5, 2.5, 0.5, 0.5))
yr <- range(mu)
for (j in seq_len(bestJ)) {
  plot(ds$times, mu[j, seq_len(d_single)], type = "l", col = "#E69F00",
       lwd = 2, ylim = yr, xlab = "", ylab = "", main = paste0("M", j))
  lines(ds$times, mu[j, d_single + seq_len(d_single)], col = "#56B4E9", lwd = 2)
  abline(h = 0, col = "grey75", lty = 3)
}
mtext("Time", side = 1, outer = TRUE, line = 1)
mtext("Estimated effect", side = 2, outer = TRUE, line = 1)
par(op)
dev.off()

## 5. trace/density diagnostics -------------------------------------------
pdf(file.path(OUT, sprintf("%s_J%d_trace.pdf", dataset, bestJ)), width = 7, height = 7)
BPFC::plot_trace_density(fit, params = c("^phi", "^psi", "^p\\[")); dev.off()

## 6. module-level network graph (clusters linked by mean-curve similarity)-
S <- cor(t(mu))                                 # J x J correlation of mean curves
adj <- S; diag(adj) <- 0; adj[abs(adj) < 0.5] <- 0
pdf(file.path(OUT, sprintf("%s_J%d_module_network.pdf", dataset, bestJ)), width = 6, height = 6)
if (requireNamespace("igraph", quietly = TRUE)) {
  g <- igraph::graph_from_adjacency_matrix(abs(adj), mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  igraph::V(g)$size <- 8 + 30 * as.numeric(table(factor(fit$clustering, levels = 1:bestJ))) / nrow(ds$Y)
  ecol <- ifelse(adj[igraph::as_edgelist(g, names = FALSE)] > 0, "#d95f02", "#1b9e77")
  plot(g, edge.width = 2 * igraph::E(g)$weight, edge.color = ecol,
       vertex.label = paste0("M", 1:bestJ), main = "Module similarity network")
} else {
  image(seq_len(bestJ), seq_len(bestJ), S, xlab = "module", ylab = "module",
        main = "Module mean-curve correlation")
}
dev.off()
message("[04] outputs written to ", OUT)
