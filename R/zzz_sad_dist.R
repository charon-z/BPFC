.bpfc_state <- new.env(parent = emptyenv())
.bpfc_state$registered <- FALSE
.bpfc_register_bivariate <- function() {
  for (nm in c("dBivSAD", "rBivSAD")) {
    assign(nm, get(nm, envir = environment(.bpfc_register_bivariate)), .GlobalEnv)
  }
  if (!.bpfc_state$registered) {
    nimble::registerDistributions(list(dBivSAD = list(
      BUGSdist = "dBivSAD(mean, phi, psi, lambda, v_sq, u_sq, times)",
      types = c("value=double(1)", "mean=double(1)", "phi=double(1)",
                "psi=double(1)", "lambda=double(1)", "v_sq=double(1)",
                "u_sq=double(1)", "times=double(1)"), discrete = FALSE)))
    .bpfc_state$registered <- TRUE
  }
  invisible(TRUE)
}
