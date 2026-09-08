#!/usr/bin/env Rscript
# 05_overfitted_selection.R ----------------------------------------------
# Reproduce the manuscript's single-run overfitted finite-mixture analysis.
# This script makes one call to run_overfitted_binary() at J_max. It does not
# call fit_many_J(); the fixed-J BIC cross-check is kept in script 04.

suppressPackageStartupMessages(library(BPFC))
src_dir <- if (nzchar(Sys.getenv("MCG_LIB_DIR"))) Sys.getenv("MCG_LIB_DIR") else
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(src_dir) || !nzchar(src_dir)) src_dir <- "BPFC/inst/reproduce"
source(file.path(src_dir, "config.R"))
OUT <- file.path(MCG_OUT_DIR, "overfitted")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dataset <- tolower(Sys.getenv("MCG_DATASET", unset = "example"))
if (!dataset %in% c("example", "mouse")) {
  stop("05_overfitted_selection.R currently supports MCG_DATASET=example or mouse.")
}

if (dataset == "mouse") {
  data_file <- system.file(
    "extdata", "mouse", "gray_mouse_snp_effects_11833x32.tsv.gz",
    package = "BPFC"
  )
  if (!nzchar(data_file)) {
    data_file <- normalizePath(
      file.path(src_dir, "..", "extdata", "mouse",
                "gray_mouse_snp_effects_11833x32.tsv.gz"),
      mustWork = FALSE
    )
  }
  if (!file.exists(data_file)) stop("Bundled mouse data not found: ", data_file)
  Y <- as.matrix(read.table(gzfile(data_file), header = TRUE, check.names = FALSE))
  storage.mode(Y) <- "double"
  rownames(Y) <- sprintf("SNP_%05d", seq_len(nrow(Y)))
  times <- 1:16
  default_J_max <- 25L
  default_niter <- 30000L
} else {
  data("example_binary", package = "BPFC", envir = environment())
  Y <- example_binary$y
  times <- example_binary$times
  default_J_max <- 6L
  default_niter <- 3000L
}

J_max <- as.integer(Sys.getenv("MCG_JMAX", unset = as.character(default_J_max)))
niter <- as.integer(Sys.getenv("MCG_NITER", unset = as.character(default_niter)))
e0 <- as.numeric(Sys.getenv("MCG_E0", unset = "0.02"))
mass_threshold <- as.numeric(Sys.getenv("MCG_MASS_THRESHOLD", unset = "0.90"))
seed <- as.integer(Sys.getenv("MCG_SEED", unset = "20260530"))

message(sprintf(
  "[05] dataset=%s n=%d d=%d J_max=%d e0=%g threshold=%g niter=%d",
  dataset, nrow(Y), ncol(Y), J_max, e0, mass_threshold, niter
))

overfit <- BPFC::run_overfitted_binary(
  Y,
  J_max = J_max,
  e0 = e0,
  mass_threshold = mass_threshold,
  times = times,
  niter = niter,
  seed = seed,
  two_phi = TRUE
)

prefix <- sprintf("%s_Jmax%d_e0_%g", dataset, J_max, e0)
saveRDS(overfit, file.path(OUT, paste0(prefix, "_overfit.rds")))
write.csv(overfit$occupied_posterior,
          file.path(OUT, paste0(prefix, "_occupied_posterior.csv")),
          row.names = FALSE)
write.csv(overfit$component_spectrum,
          file.path(OUT, paste0(prefix, "_component_spectrum.csv")),
          row.names = FALSE)

feature_ids <- rownames(Y)
if (is.null(feature_ids)) feature_ids <- sprintf("feature_%05d", seq_len(nrow(Y)))
assignments <- data.frame(
  feature = feature_ids,
  cluster = overfit$clustering,
  overfit$cluster_prob,
  uncertainty = overfit$cluster_uncertainty,
  retained_probability = overfit$retained_probability,
  check.names = FALSE
)
write.csv(assignments,
          file.path(OUT, paste0(prefix, "_dominant_assignments.csv")),
          row.names = FALSE)

pdf(file.path(OUT, paste0(prefix, "_selection_diagnostics.pdf")),
    width = 13, height = 4.5)
BPFC::plot_overfitted_selection(overfit)
dev.off()

summary_row <- data.frame(
  dataset = dataset,
  n = nrow(Y),
  d = ncol(Y),
  J_max = J_max,
  e0 = e0,
  mass_threshold = mass_threshold,
  selected_K = overfit$selected_K,
  K_occ_mode = overfit$K_occ_mode,
  K_occ_lower = overfit$K_occ_interval["lower"],
  K_occ_upper = overfit$K_occ_interval["upper"],
  niter = niter,
  seed = seed
)
write.csv(summary_row, file.path(OUT, paste0(prefix, "_summary.csv")),
          row.names = FALSE)
print(overfit)
message("[05] outputs written to ", OUT)
