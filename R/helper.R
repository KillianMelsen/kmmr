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

#' getDeviceProps
#' 
#' The \code{getDeviceProps} function returns some properties of the CUDA
#' device specified by the \code{id} argument.
#'
#' @param id ID of the CUDA device for which properties should be returned.
#' Defaults to \code{0} for the first device.
#'
#' @export
getDeviceProps = function(id = 0) {
  
  rst <- .getDeviceProps(id = id)
  
  txt <- c("Found CUDA capable device.\n",
           "Device ID:\t\t%s\n",
           "Device name:\t\t%s\n",
           "Integraded device:\t%s\n",
           "Compute capability:\t%s\n")
  
  msg <- sprintf(paste(txt, collapse = ""),
                 rst$id,
                 rst$deviceName,
                 as.logical(rst$integr),
                 sprintf("%s.%s", rst$mjr, rst$mnr))
  cat(msg)
}


