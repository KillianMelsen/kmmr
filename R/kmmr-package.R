#' @keywords internal
#' @import Matrix
#' @importFrom Rcpp sourceCpp
#' @importFrom methods as
#' @importFrom Matrix nearPD Matrix t tcrossprod crossprod sparse.model.matrix Cholesky chol forceSymmetric diag chol2inv solve rowSums
#' @importFrom stats na.omit terms as.formula
#' @importFrom fad fad
#' @importFrom tictoc tic toc tic.log tic.clearlog
"_PACKAGE"

## usethis namespace: start
#' @useDynLib kmmr, .registration = TRUE
## usethis namespace: end
NULL
