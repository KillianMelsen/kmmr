#' predict.kmmrFit
#'
#' @description
#' \code{predict.kmmrFit} returns predictions and associated standard errors for
#' any linear combination of the BLUEs and BLUPs estimated from a linear mixed model
#' fitted using \code{kmmr}. It is based on the idea that any predicted value can
#' be computed given a prediction column vector \eqn{\mathbf{k}} of coefficients associated with each BLUE
#' or BLUP as \eqn{\mathbf{k}^\top \boldsymbol{\tau}} where \eqn{\boldsymbol{\tau} = [\boldsymbol{\beta}^\top \mathbf{u}^\top]^\top}.
#' The standard error of that predicted value is then \eqn{\mathbf{k}^\top \mathbf{C}^{-1} \mathbf{k}}
#' which can be efficiently computed using the (sparse) Cholesky factor of \eqn{\mathbf{C}}
#' by solving the system \eqn{\mathbf{C}\mathbf{x} = \mathbf{k}} and subsequently computing
#' \eqn{\mathbf{k}^\top \mathbf{x}}. The actual inverse \eqn{\mathbf{C}^{-1}} is never formed.
#' 
#' For making many predictions at the same time, individual prediction vectors can be
#' combined into a sparse prediction matrix so that all predictions are obtained using
#' \eqn{\mathbf{K}^\top \boldsymbol{\tau}} and standard errors are obtained as the diagonal
#' elements of \eqn{\mathbf{K}^\top \mathbf{C}^{-1} \mathbf{K}}.
#'
#' @param object kmmrFit object.
#' @param predictionMatrix A sparse vector or matrix of prediction coefficients.
#' @param ... Extra arguments.
#'
#' @keywords internal
#' @exportS3Method predict kmmrFit
predict.kmmrFit <- function(object, predictionMatrix, chunks = FALSE, chunkSize = 500, ...) {
  if (is(predictionMatrix, "sparseVector")) {
    predictionMatrix <- as(predictionMatrix, "sparseMatrix")
  }
  tau <- rbind(object$b, object$u)
  Cmode <- object$infoPEV$Cmode
  pvals <- crossprod(predictionMatrix, tau) # Predicted values are easy.
  if (Cmode == "mixed") {
    # Preparing the data:
    dense.idx <- object$infoPEV$dense.idx
    sparse.idx <- object$infoPEV$sparse.idx
    kDense <- predictionMatrix[dense.idx, , drop = FALSE]
    kSparse <- predictionMatrix[sparse.idx, , drop = FALSE]
    
    # Sparse part:
    kSparsePermuted <- solve(object$infoPEV$Ccholss, kSparse, system = "P")
    
    # Do the entire thing at once (needs more memory), or divide into chunks?
    if (!chunks) {
      gamma <- solve(object$infoPEV$Ccholss, as.matrix(kSparsePermuted), system = "L")
      pevSparse <- colSums(gamma * gamma) # diag(crossprod(gamma))
    } else {
      n_preds <- ncol(kSparsePermuted)
      pevSparse <- numeric(n_preds)
      for (i in seq(1, n_preds, by = chunkSize)) {
        idx <- i:min(i + chunkSize - 1, n_preds)
        k_chunk <- as.matrix(kSparsePermuted[, idx, drop = FALSE])
        gamma_chunk <- solve(object$infoPEV$Ccholss, k_chunk, system = "L")
        pevSparse[idx] <- colSums(gamma_chunk * gamma_chunk)
      }
    }
    
    # Dense part:
    kDenseAdjusted <- kDense - (object$infoPEV$Qt %*% kSparse)
    z <- object$infoPEV$Cidd %*% kDenseAdjusted
    pevDense <- colSums(kDenseAdjusted * z)
    pev <- pevSparse + pevDense
    
  } else if (Cmode == "sparse") {
    sparse.idx <- object$infoPEV$sparse.idx
    tauSparse <- tau[sparse.idx, 1, drop = FALSE]
    kSparse <- predictionMatrix[sparse.idx, , drop = FALSE]
    kSparsePermuted <- solve(object$infoPEV$Ccholss, kSparse, system = "P")
    
    # Do the entire thing at once (needs more memory), or divide into chunks?
    if (!chunks) {
      gamma <- solve(object$infoPEV$Ccholss, as.matrix(kSparsePermuted), system = "L")
      pevSparse <- colSums(gamma * gamma) # diag(crossprod(gamma))
    } else {
      n_preds <- ncol(kSparsePermuted)
      pevSparse <- numeric(n_preds)
      for (i in seq(1, n_preds, by = chunkSize)) {
        idx <- i:min(i + chunkSize - 1, n_preds)
        k_chunk <- as.matrix(kSparsePermuted[, idx, drop = FALSE])
        gamma_chunk <- solve(object$infoPEV$Ccholss, k_chunk, system = "L")
        pevSparse[idx] <- colSums(gamma_chunk * gamma_chunk)
      }
    }
    pev <- pevSparse
    
  } else if (Cmode == "dense") {
    dense.idx <- object$infoPEV$dense.idx
    tauDense <- tau[dense.idx, 1, drop = FALSE]
    kDense <- predictionMatrix[dense.idx, , drop = FALSE]
    z <- object$infoPEV$Cidd %*% kDense
    pevDense <- colSums(kDense * z)
    pev <- pevDense
    
  }
  
  return(data.frame(Prediction = pvals@x,
                    SE = sqrt(pev)))
}

#' summary.kmmrFit
#'
#' @param object kmmrFit object.
#' @param round Number of decimals to be printed for parameter estimates and SEs.
#' @param ... Extra arguments.
#' 
#' @keywords internal
#' @exportS3Method summary kmmrFit
summary.kmmrFit <- function(object, round = 5, ...) {
  estimates <- object$k[, ncol(object$k)]
  ses <- rep(NA, length(estimates))
  ses[which(object$isFixed == 0)] <- sqrt(diag(solve(object$AI[which(object$isFixed == 0), which(object$isFixed == 0)])))
  fixed <- ifelse(object$isFixed == 0, "No", "Yes")
  
  outputData <- list(logLik = object$LL[length(object$LL)],
                     varcomp = data.frame(Parameter = object$params,
                                          Estimate = estimates,
                                          SE = ses,
                                          Fixed = fixed))
  attr(outputData, "round") <- round
  class(outputData) <- c("summary.kmmrFit", "list")
  return(outputData)
  # output <- apply(outputData, 1, function(row) {
  #   sprintf("%-20s %-12.8f %-12.8f %-8s",
  #           row["parameter"],
  #           as.numeric(row["estimate"]),
  #           as.numeric(row["SE"]),
  #           row["isFixed"])
  # })
  # 
  # header <- sprintf("%-20s %-12s %-12s %-8s", "Parameter", "Estimate", "SE", "Fixed?")
  # ansi <- interactive() && !is.null(Sys.getenv("TERM"))
  # if (ansi) {
  #   cat("\033[1m", header, "\033[0m", "\n", sep = "")
  # } else {
  #   cat(header, "\n")
  # }
  # cat(paste(output, collapse = "\n"))
}

#' print.summary.kmmrFit
#'
#' @param x summary.kmmrFit object.
#' @param ... Extra arguments.
#' 
#' @keywords internal
#' @exportS3Method print summary.kmmrFit
print.summary.kmmrFit <- function(x, ...) {
  # Is bold text possible?
  ansi <- interactive() && !is.null(Sys.getenv("TERM"))
  round <- attr(x, "round")
  
  # Formatting parameter estimates, SEs, etc.
  max_len <- max(nchar(x$varcomp$Parameter))
  output <- apply(x$varcomp, 1, function(row) {
    sprintf(gettextf("%%-%ds %%-12.%df %%-12.%df %%-8s", max_len + 5, round, round),
            row["Parameter"],
            as.numeric(row["Estimate"]),
            as.numeric(row["SE"]),
            row["Fixed"])
  })
  
  # Formatting parameter table header:
  header <- sprintf(gettextf("%%-%ds %%-12s %%-12s %%-8s", max_len + 5), "Parameter", "Estimate", "SE", "Fixed?")

  # Printing:
  if (ansi) {
    cat("\033[1m", "Log-likelihood: ", "\033[0m", sprintf("%.3f", x$logLik), "\n\n", sep = "")
    cat("\033[1m", header, "\033[0m", "\n", sep = "")
  } else {
    cat("Log-likelihood: ", sprintf("%.3f", x$logLik), "\n\n", sep = "")
    cat(header, "\n")
  }
  cat(paste(output, collapse = "\n"), "\n")
}




