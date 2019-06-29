# generate roxygen namespace file

# old namespace backup
# useDynLib(temp.spat.pca3, .registration=TRUE)
# importFrom(Rcpp, evalCpp)
# exportPattern("^[[:alpha:]]+")

# new namespace
#' @useDynLib temp.spat.pca3
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL