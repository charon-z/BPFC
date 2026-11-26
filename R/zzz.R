.onLoad <- function(libname, pkgname) {
  utils::getFromNamespace(".bpfc_register_bivariate", pkgname)()
}
