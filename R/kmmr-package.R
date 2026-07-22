#' @useDynLib kmmr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom methods as
#' @importFrom Matrix nearPD Matrix t tcrossprod crossprod sparse.model.matrix Cholesky chol forceSymmetric diag chol2inv
#' @importFrom stats na.omit terms as.formula
NULL