#' exfalkSetup
#'
#' @param K Kinship matrix.
#' @param m Number of factors to be used.
#' @param method Which method should be used to estimate loadings and unique variances?
#'
#' @returns A list with the number of factors, scaled loadings matrix, and scaled inverse unique variances.
#' @export
exfalkSetup <- function(K, m = 50, method = "factanal") {
  if (!is.matrix(K)) {
    K <- as.matrix(K)
  }
  scaling <- sqrt(diag(K))
  if (method == "factanal") {
    fit <- factanal(covmat = K, factors = m)
  } else if (method == "fad") {
    fit <- fad(covmat = K, factors = m)
    rownames(fit$loadings) <- rownames(K)
  } else {
    stop("No valid method!")
  }
  L <- scaling * Matrix(fit$loadings[,, drop = FALSE])
  PSIi <- Diagonal(nrow(L), 1 / (fit$uniquenesses * scaling^2))
  return(list(m = ncol(L),
              L = L,
              PSIi = PSIi))
}
