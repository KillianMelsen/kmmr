timer <- function(command, phase, iter) {
  if (command == "start") {
    tic(sprintf("%s:iter_%d", phase, iter))
  } else if (command == "stop") {
    toc(log = TRUE, quiet = TRUE)
  } else if (command == "report") {
    log <- tic.log()
    timings <- data.frame(phase = gsub("(.*):iter_(.*):.*", "\\1", unlist(tic.log())),
                          iter = gsub("(.*):iter_(.*):.*", "\\2", unlist(tic.log())),
                          duration = as.numeric(unlist(lapply(tic.log(format = FALSE), function(x) x$toc - x$tic))))
    timings$phase <- factor(timings$phase, unique(timings$phase))
    tic.clearlog()
    return(timings)
  }
}

formGi <- function(sigmai, nvarGi) {
  # sigmai is a 1x1 matrix:
  if (identical(dim(sigmai), c(1L, 1L))) {
    return(sigmai@x * nvarGi)
  }
  
  # nvarGi is a 1x1 matrix, don't think this should ever happen:
  if (identical(dim(nvarGi), c(1L, 1L))) {
    stop("One of the nvarGi matrices is a 1x1 matrix!")
    # return(sigmai * nvarGi@x)
  }
  
  # General case for p > 1 (dense or sparse)
  res <- kronecker(sigmai, nvarGi)
  
  # Making sure the result is a ds.Matrix:
  if (!is(res, "symmetricMatrix")) {
    res <- forceSymmetric(res)
  }
  
  return(res)
}

# A function for checking and constraining updates to the loadings of fa and exfa structures:
checkLoadings <- function(currentParams, delta, loadingIndices, isFixed, p) {
  kUpdate <- currentParams + delta
  signflip <- sign(currentParams[loadingIndices, 1]) * sign(kUpdate[loadingIndices, 1]) == -1 # Detecting sign flips.
  mags <- abs(kUpdate / currentParams[, 1])[loadingIndices, 1] # Computing the magnitudes of the changes.
  meanMags <- mean(mags[!isFixed[loadingIndices]]) # Should be smaller than 10.
  con.loadings <- c()
  bound.loadings <- c()
  constrained <- FALSE
  checkConstraints <- FALSE
  constraintsUsed <- FALSE
  hitBoundary <- FALSE
  if (meanMags > 10) { # Shrink delta for all loadings.
    con.loadings <- loadingIndices[!isFixed[loadingIndices]] # Loadings that are fixed should not be considered as constrained.
    scaling <- 9 / meanMags # Close enough
    delta[con.loadings, 1] <- scaling * delta[con.loadings, 1]
    constrained <- TRUE
    checkConstraints <- TRUE
    constraintsUsed <- TRUE
    
    # We have changed delta, so we need to recompute the magnitudes (mean magnitude should now be lower than 10).
    # We also need to recheck the sign-flips, because by shrinking delta some of them could have disappeared.
    kUpdate <- currentParams + delta # Easiest way is to start with recomputing the proposed update
    signflip <- sign(currentParams[loadingIndices, 1]) * sign(kUpdate[loadingIndices, 1]) == -1 # Detecting sign flips.
    mags <- abs(kUpdate / currentParams[, 1])[loadingIndices, 1] # Computing the magnitudes of the changes.
  }
  if (any(mags[signflip] > 2)) { # Should all be smaller than 2.
    # We need the sort(unique()) because if the mean magnitude was larger than 10, all loadings were already constrained.
    # But some of them might need to be constrained further (sign-flipped and mag > 2).
    # Subset loading indices and mags to sign-flipped ones, then subset to those for which the magnitude > 2.
    magOver2 <- loadingIndices[signflip][mags[signflip] > 2] # These can never be fixed loadings, because fixed loadings never change.
    con.loadings <- sort(unique(c(con.loadings, magOver2)))
    # We need the -1e-9 because if we use exactly 2 we can run into numerical issues where it keeps saying that mags[signflip] > 2 == TRUE even after constraining...
    delta[magOver2, 1] <- -(2 - 1e-9) * currentParams[magOver2, 1] - currentParams[magOver2, 1]
    constrained <- TRUE
    checkConstraints <- TRUE
    constraintsUsed <- TRUE
    
    # We have changed delta so we need to recompute the proposed update.
    # We don't need to recompute the magnitudes, because the mean magnitude can only
    # go down if we further constrain updates for sign-flipped loadings.
    kUpdate <- currentParams + delta
  }
  # Now the final part where we check whether all communalities are smaller than 1 - 1e-3:
  L <- Matrix(kUpdate[loadingIndices, 1], nrow = p, sparse = FALSE, doDiag = FALSE)
  comms <- rowSums(L^2)
  m <- ncol(L) # Number of factors.
  if (any(comms > 1 - 1e-3)) {
    commOver1 <- which(comms > 1 - 1e-3) # These are the rows of L for which comm > 1 - 1e-3.
    # We need the indices of all loadings involved in these rows:
    loadingIndicesMat <- Matrix(loadingIndices, p, m)
    # Add them to the full vector of constrained, but remove any fixed loadings using setdiff().
    bound.loadings <- sort(unique(c(bound.loadings, setdiff(loadingIndicesMat[commOver1,, drop = FALSE]@x, loadingIndices[isFixed[loadingIndices]]))))
    for (row in commOver1) {
      # We need to make sure that the scaling results in communalities that are slightly below the threshold of 1 - 1e-3.
      # This is purely due to machine precision issues. If we compute the scaling factor exactly as sqrt((1 - 1e-3) / comms[row])
      # the communalities will still be slightly over 0.999. They will be something like 1e-15 too large...
      scaling <- sqrt((1 - 1e-3 - 1e-9) / comms[row])
      delta[loadingIndicesMat[row,], 1] <- scaling * kUpdate[loadingIndicesMat[row,], 1] - currentParams[loadingIndicesMat[row,], 1]
    }
    hitBoundary <- TRUE
    checkConstraints <- TRUE
    constraintsUsed <- TRUE
  }
  return(list(delta = delta,
              constrained = constrained,
              checkConstraints = checkConstraints,
              constraintsUsed = constraintsUsed,
              hitBoundary = hitBoundary,
              con.loadings = con.loadings,
              bound.loadings = bound.loadings,
              isFixed = isFixed))
}

getCormat_corg <- function(p, params, doCormat = TRUE, doCormati = TRUE, doLogDet = TRUE) {
  if (p != (1 + sqrt(1 + 8 * length(params))) / 2) {
    stop("Number of factor levels and parameters for corg model do not match.")
  }
  cormat <- logDet <- cormati <- NULL
  if (doCormat) {
    cormat <- Matrix(diag(p), sparse = FALSE, doDiag = FALSE)
    cormat[lower.tri(cormat)] <- params
    cormat <- forceSymmetric(cormat, uplo = "L")
  }
  if (doLogDet) {
    cormatChol <- chol(cormat)
    logDet <- 2 * sum(log(diag(cormatChol)))
  }
  if (doCormati) {
    cormati <- chol2inv(cormatChol)
  }
  return(list(cormat, logDet, cormati))
}

getCormat_cors <- function(p, params, doCormat = TRUE, doCormati = TRUE, doLogDet = TRUE) {
  cormat <- logDet <- cormati <- NULL
  if (doCormat) {
    cormat <- Matrix(diag(p), sparse = FALSE, doDiag = FALSE)
    cormat[lower.tri(cormat)] <- params[1]
    cormat <- forceSymmetric(cormat, uplo = "L")
  }
  if (doChol) {
    cormatChol <- chol(cormat)
    logDet <- 2 * sum(log(diag(cormatChol)))
  }
  if (doCormati) {
    cormati <- chol2inv(cormatChol)
  }
  return(list(cormat, logDet, cormati))
}

getCormat_gk <- function(params, base, doCormat = TRUE, doCormati = TRUE, doLogDet = TRUE) {
  cormat <- logDet <- cormati <- NULL
  if (doCormat) {
    cormat <- exp(-params[1] * base)
  }
  if (doLogDet) {
    cormatChol <- chol(cormat)
    logDet <- 2 * sum(log(diag(cormatChol)))
  }
  if (doCormati) {
    cormati <- chol2inv(cormatChol)
  }
  return(list(cormat, logDet, cormati))
}

getCormat_ar1 <- function(params, base, doCormat = TRUE, doCormati = TRUE, doLogDet = TRUE) {
  cormat <- logDet <- cormati <- NULL
  n <- nrow(base)
  if (doCormat) {
    cormat <- params[1]^base
  }
  if (doLogDet) {
    superDiag <- base@x[(1:(n - 1)) * (n + 1)]
    logDet <- sum(log(1 - params[1]^(2 * superDiag)))
  }
  if (doCormati) {
    superDiag <- params[1]^base@x[(1:(n - 1)) * (n + 1)]
    denom <- 1 - superDiag^2
    a <- superDiag / denom
    b <- superDiag^2 / denom
    c <- 1 / denom
    d <- c(c[1], c[-1] + b[-length(b)], c[length(c)])
    cormati <- Matrix(0, n, n, sparse = TRUE, doDiag = FALSE)
    diag(cormati) <- d
    cormati[cbind(1:(n - 1), 2:n)] <- -a
    cormati <- forceSymmetric(cormati)
  }
  return(list(cormat, logDet, cormati))
}

getCormat_lk <- function(base, doCormat = TRUE, doCormati = TRUE, doLogDet = TRUE) {
  cormat <- logDet <- cormati <- NULL
  if (doCormat) {
    cormat <- base
  }
  if (doLogDet) {
    cormatChol <- chol(cormat)
    logDet <- 2 * sum(log(diag(cormatChol)))
  }
  if (doCormati) {
    cormati <- chol2inv(cormatChol)
  }
  return(list(cormat, logDet, cormati))
}

getCormat_id <- function(p, doCormat = TRUE, doCormati = TRUE, doLogDet = TRUE) {
  cormat <- logDet <- cormati <- NULL
  if (doCormat) {
    cormat <- Matrix(diag(p), sparse = TRUE, doDiag = TRUE)
  }
  if (doLogDet) {
    cormatChol <- Matrix(diag(p), sparse = TRUE, doDiag = TRUE)
    logDet <- 0
  }
  if (doCormati) {
    cormati <- Matrix(diag(p), sparse = TRUE, doDiag = TRUE)
  }
  return(list(cormat, logDet, cormati))
}

getCormat_fa <- function(p, params, doCormat = TRUE, doCormati = TRUE, doLogDet = TRUE) {
  cormat <- logDet <- cormati <- NULL
  if (doCormat) {
    L <- Matrix(params, nrow = p, sparse = FALSE, doDiag = FALSE)
    # tcrossprod(L) + (diag(p) - diag(p) * tcrossprod(L)), second term is equal to PSI
    cormat <- tcrossprod(L)
    diag(cormat) <- 1
    cormat <- forceSymmetric(cormat)
  }
  if (doLogDet) {
    cormatChol <- chol(cormat)
    logDet <- 2 * sum(log(diag(cormatChol)))
  }
  if (doCormati) {
    cormati <- chol2inv(cormatChol)
  }
  return(list(cormat, logDet, cormati))
}

getCormat_exfa <- function(p, params, doCormat = TRUE, doCormati = TRUE, doLogDet = TRUE) {
  cormat <- logDet <- cormati <- NULL
  if (doCormat) {
    L <- Matrix(params, nrow = p, sparse = FALSE, doDiag = FALSE)
    m <- ncol(L)
    LLt <- tcrossprod(L)
    cormat <- rbind(cbind(LLt, -L),
                    cbind(-t(L), Matrix(diag(m), sparse = FALSE, doDiag = FALSE)))
    diag(cormat) <- 1
    cormat <- forceSymmetric(cormat)
  }
  if (doLogDet) {
    cormatChol <- chol(cormat)
    logDet <- 2 * sum(log(diag(cormatChol)))
  }
  if (doCormati) {
    PSIi <- 1 / (1 - diag(LLt)) # diag(solve(diag(p) - diag(p) * Matrix::tcrossprod(L)))
    PSIi_L <- PSIi * L
    cormati <- rbind(cbind(diag(PSIi), PSIi_L),
                     cbind(t(PSIi_L), diag(m) + crossprod(L, PSIi_L)))
    cormati <- forceSymmetric(as(cormati, "CsparseMatrix"))
  }
  return(list(cormat, logDet, cormati))
}

#' getCormat forms, for a given structure and set of parameters, the correlation matrix, inverse of the correlation matrix, and determinant of the correlation matrix.
#'
#' @param struc Character indicating which correlation structure must be formed from the parameters.
#' @param p Numeric indicating the number of levels of this factor.
#' @param params Numeric vector containing the parameters needed to form the correlation matrix.
#' @param base Optional argument containing a numeric matrix to be used in forming the correlation matrix (distance matrix for Gaussian kernel, correlation matrix or kinship for linear kernel).
#' @param doCormat Logical indicating whether the correlation matrix should be formed.
#' @param doCormati Logical indicating whether the inverse of the correlation matrix should be formed.
#' @param doLogDet Logical indicating whether the log determinant should be computed.
#'
#' @returns A list containing the formed correlation matrix as well as its upper triangular Cholesky decomposition.
#' @keywords internal
getCormat <- function(struc, p = NULL, params = NULL, base = NULL, doCormat = TRUE, doCormati = TRUE, doLogDet = TRUE) {
  if (struc == "corg") {
    ret <- getCormat_corg(p = p, params = params, doCormat = doCormat, doCormati = doCormati, doLogDet = doLogDet)
  } else if (struc == "cors") {
    ret <- getCormat_cors(p = p, params = params, doCormat = doCormat, doCormati = doCormati, doLogDet = doLogDet)
  } else if (struc == "gk") {
    ret <- getCormat_gk(params = params, base = base, doCormat = doCormat, doCormati = doCormati, doLogDet = doLogDet)
  } else if (struc == "lk") {
    ret <- getCormat_lk(base = base, doCormat = doCormat, doCormati = doCormati, doLogDet = doLogDet)
  } else if (struc == "id") {
    ret <- getCormat_id(p = p, doCormat = doCormat, doCormati = doCormati, doLogDet = doLogDet)
  } else if (struc == "fa") {
    ret <- getCormat_fa(p = p, params = params, doCormat = doCormat, doCormati = doCormati, doLogDet = doLogDet)
  } else if (struc == "exfa") {
    ret <- getCormat_exfa(p = p, params = params, doCormat = doCormat, doCormati = doCormati, doLogDet = doLogDet)
  } else if (struc == "ar1") {
    ret <- getCormat_ar1(params = params, base = base, doCormat = doCormat, doCormati = doCormati, doLogDet = doLogDet)
  }
  return(ret)
}

getDeriv_cor <- function(cormatIndex, p) {
  deriv <- Matrix(0, p, p, sparse = FALSE, doDiag = FALSE)
  deriv[lower.tri(deriv)][cormatIndex] <- 1
  deriv <- forceSymmetric(deriv, uplo = "L")
  return(as(deriv, "generalMatrix"))
}

getDeriv_gk <- function(params, base) {
  deriv <- -base * exp(-params[1] * base)
  return(deriv)
}

getDeriv_ar1 <- function(params, base) {
  if (params[1] == 0) {
    deriv <- 1 * (base == 1)
  } else {
    deriv <- base * params[1]^(base - 1)
  }
  return(deriv)
}

getDeriv_fa <- function(cormatIndex, params, p) {
  L <- Matrix(params, nrow = p, sparse = FALSE, doDiag = FALSE)
  m <- ncol(L)
  dL <- Matrix(0, nrow = p, ncol = m, sparse = TRUE, doDiag = FALSE)
  dL[cormatIndex] <- 1
  # This is equal to dLL^T/dl + dPSI/dl, where dLL^T/dl = L %*% (dL/dl)^T + dL/dl %*% L^T and dPSI/dl = -I * dLL^t/dl, easier just to set the diagonal to 0...
  deriv <- tcrossprod(L, dL) + tcrossprod(dL, L) # L %*% t(dL) + dL %*% t(L)
  diag(deriv) <- 0
  # deriv <- forceSymmetric(as(deriv, "CsparseMatrix"))
  return(deriv)
}

getDeriv_exfa <- function(cormatIndex, params, p, inverse) {
  L <- Matrix(params, nrow = p, sparse = FALSE, doDiag = FALSE)
  m <- ncol(L)
  dL <- Matrix(0, nrow = p, ncol = m, sparse = TRUE, doDiag = FALSE)
  dL[cormatIndex] <- 1
  if (!inverse) {
    # See getDeriv_fa, ignore dPSI/dl and just set the diagonal elements to 0:
    deriv <- rbind(cbind(tcrossprod(L, dL) + tcrossprod(dL, L), -dL),
                   cbind(-t(dL), Matrix(0, m, m, sparse = FALSE, doDiag = FALSE)))
    diag(deriv) <- 0
    # deriv <- forceSymmetric(as(deriv, "CsparseMatrix"))
  } else { # If the derivative of the inverse extended factor analytic correlation matrix is needed:
    PSIi <- Matrix(diag(1 / (1 - diag(tcrossprod(L)))), sparse = TRUE, doDiag = FALSE)
    dPSI <- Matrix(diag(-diag(tcrossprod(L, dL) + tcrossprod(dL, L))), sparse = TRUE, doDiag = FALSE)
    dPSIi <- forceSymmetric(-PSIi %*% dPSI %*% PSIi)
    # dPSIi <- diag(PSIi)^2 * dPSI
    deriv <- rbind(cbind(dPSIi, PSIi %*% dL + dPSIi %*% L),
                   cbind(t(PSIi %*% dL + dPSIi %*% L), t(dL) %*% PSIi %*% L + t(L) %*% dPSIi %*% L + t(L) %*% PSIi %*% dL))
    deriv <- forceSymmetric(deriv)
  }
  return(deriv)
}

getDeriv <- function(struc, type, cormatIndex = NULL, params = NULL, p = NULL, base = NULL, inverse) {
  if (struc == "fa") {
    ret <- getDeriv_fa(cormatIndex = cormatIndex, params = params, p = p)
  } else if (struc == "exfa") {
    ret <- getDeriv_exfa(cormatIndex = cormatIndex, params = params, p = p, inverse = inverse)
  } else if (struc == "gk") {
    ret <- getDeriv_gk(params = params, base = base)
  } else if (struc %in% c("corg", "cors")) {
    ret <- getDeriv_cor(cormatIndex = cormatIndex, p = p)
  } else if (struc == "ar1") {
    ret <- getDeriv_ar1(params = params, base = base)
  }
  return(ret)
}

#' solveLMM
#' 
#' Implementation of the C-based average information REML algorithm for fitting linear mixed models.
#'
#' @param y Numeric vector containing the response data.
#' @param X Numeric fixed effect design matrix.
#' @param Z List containing numeric random effect design matrices. One for each random effect.
#' @param init Numeric vector of initial values for the parameters.
#' @param weights Numeric matrix of weights. This can be based on the full covariance matrix $\Omega^{-1}$.
#' Alternatively, for a diagonal approximation, the weights can be based on the reciprocals of known variances, i.e., the diagonal entries of $\Omega^{-1}$ or plot numbers.
#' Note that 'weights' is an upper triangular Cholesky factor if the full $\Omega^{-1}$ is used, or a diagonal matrix containing the square roots of the plot numbers or diagonal entries of $\Omega^{-1}$.
#' @param useWeights Logical indicating whether to use weights or not.
#' @param maxit Numeric value indicating the maximum number of iterations for the AI algorithm.
#' @param stepSize Numeric vector indicating the scaling factor for parameter updates for each iteration.
#' @param avInfReg Numeric vector indicating the penalty to be applied to the average information matrix. Using \code{1.1} inflates the diagonal
#' elements of the AI matrix by 10%.
#' @param verbose  Logical indicating verbosity.
#' @param tolConvLL Numeric value indicating the tolerance for deciding on convergence based on the change in log-likelihood.
#' @param indices Indices used for selecting the correct parameters for the covariance matrices of random effects.
#' @param u.applyvars To which factors should variances be applied?
#' @param ur.types Parameter types.
#' @param u.nlevels Numbers of levels for all factors.
#' @param u.forFactor For which factor is each parameter.
#' @param u.indexCor Indices used to count correlations.
#' @param u.sparse Provides info on whether a random effect is sparse or not.
#' @param bases Base matrices
#' @param struc Structure for each factor.
#' @param i.var Which factors have variable covariance structures?
#' @param i.nvar Which factors have non-variable covariance structures?
#' @param isFixed Which parameters are fixed and should not be estimated?
#' @param u.isStack Are factors within random effects made up of a stack of structures?
#' @param u.nparams.list Number of parameters per factor/random effect in list format.
#' @param timesToBoundary How many times is a parameter or set of loadings allowed to cross the boundary before it is fixed?
#' @param super Should a supernodal (\code{super = TRUE}) or simplicial Cholesky factor be computed?
#' @param useGPU Should the GPU be used for dense models?
#' @param denseThreshold Threshold for determining whether a fixed effect is dense based on the proportion of non-zeros in each column of \code{X}.
#' @param debugging If \code{TRUE}, additional information such as derivatives, log determinants, etc., is returned for debugging purposes.
#' @param timing Can be set to \code{TRUE} to report runtimes.
#'
#' @returns A list containing the log-likelihood, parameter estimates, average information matrix,
#' scores, BLUEs for fixed effects, and BLUPs for random effects.
#' @keywords internal
solveLMM <- function(y, X, Z, init, weights, useWeights, maxit, stepSize, avInfReg, verbose, tolConvLL,
                     indices, u.applyvars, ur.types, u.nlevels, u.nparams, u.forFactor, u.indexCor, u.sparse,
                     bases, struc, i.var, i.nvar, isFixed, u.isStack, u.nparams.list,
                     timesToBoundary, super, useGPU, denseThreshold, debugging, timing) {
  ## Prelims ====
  start <- proc.time()[3]
  k <- Matrix(0, length(init), maxit + 1, sparse = FALSE, doDiag = FALSE) # Matrix of parameter estimates.
  k[, 1] <- init # Filling with initial values.
  boundary <- numeric(length(init)) # Vector to hold how many times each parameter has crossed the boundary.
  n <- nrow(y)
  
  nU <- length(Z) - 1 # Number of random effects
  nR <- 1 # Number of residual effects
  nUR <- nU + nR # Number of random effects + residual
  nX <- ncol(X) # Number of fixed effects
  
  # Determining the number of dense and sparse fixed effects:
  colDensitiesX <- diff(X@p) / nrow(X) # Proportion of nonzeros for each column of X
  nXdense <- sum(colDensitiesX > denseThreshold) # Number of dense fixed effects
  nXsparse <- nX - nXdense # Number of sparse fixed effects
  Xdense.idx <- which(colDensitiesX > denseThreshold)
  Xsparse.idx <- setdiff(1:nX, Xdense.idx)
  
  # Doing something similar for the random effects:
  Z.idx <- lapply(Z[1:nU], function(Z) return(1:ncol(Z)))
  Z.idx[[1]] <- Z.idx[[1]] + nX
  if (length(Z.idx) > 1) {
    for (i in 2:length(Z.idx)) {
      Z.idx[[i]] <- Z.idx[[i]] + max(Z.idx[[i - 1]])
    }
  }
  nZdense <- sum(!(u.sparse[1:nU])) # Ignore the residual structure.
  nZsparse <- nU - nZdense
  Zdense.idx <- which(!(u.sparse[1:nU]))
  Zsparse.idx <- setdiff(1:nU,Zdense.idx)
  
  ## Computing determinants of matrices associated with non-variable terms before iterating ====
  n.var <- numeric(nUR)
  n.nvar <- numeric(nUR)
  dim.var <- numeric(nUR)
  dim.nvar <- numeric(nUR)
  u.nlevels.var <- vector("list", nUR) 
  u.nlevels.nvar <- vector("list", nUR) 
  for (re in 1:nUR) { # Start loop over random effects and residual effect
    n.var[re] <- length(i.var[[re]]) # How many variable factors are there?
    n.nvar[re] <- length(i.nvar[[re]]) # How many non-variable factors are there?
    u.nlevels.var[[re]] <- u.nlevels[[re]][i.var[[re]]] # Vector of levels for each variable factor for each random effect.
    u.nlevels.nvar[[re]] <- u.nlevels[[re]][i.nvar[[re]]] # Vector of levels for each non-variable factor for each random effect.
    dim.var[re] <- prod(u.nlevels.var[[re]]) # Total dimension of the variable part of each random effect.
    dim.nvar[re] <- prod(u.nlevels.nvar[[re]]) # Total dimension of the non-variable part of each random effect.
  } # End loop over random effects
  
  # Forming non-variable structures for random effects:
  nvarGi <- vector("list", nU)
  logDetnvarGtotal <- numeric(nU)
  for (re in 1:nU) { # Start loop over random effects
    if (n.nvar[re] > 0) { # All of this can only be done if there are non-variable factors (such as a kinship/linear kernel for genotypes).
      logDetnvarG <- numeric(n.nvar[re]) # A log determinant for each non-variable base matrix (linear kernel).
      for (i in i.nvar[[re]]) {
        if (!is.null(attr(bases[[re]][[i]], "logDet"))) {
          logDetnvarG[i - n.var[re]] <- attr(bases[[re]][[i]], "logDet")
        } else {
          if (is(bases[[re]][[i]], "sparseMatrix")) {
            logDetnvarG[i - n.var[re]] <- -1 * determinant(Cholesky(bases[[re]][[i]], LDL = FALSE, perm = TRUE, super = TRUE), sqrt = FALSE, logarithm = TRUE)$modulus
          } else {
            logDetnvarG[i - n.var[re]] <- -1 * sum(log(eigen(bases[[re]][[i]], only.values = TRUE)$values)) # Times -1 because the non-variable base matrices (linear kernels) are in inverse form. # Maybe on the GPU for large matrices? 
          }
        }
      }
      
      # If the non-variable part of G consists of more than 1 factor we need to compute the total log determinant:
      if (n.nvar[re] > 1) {
        # u.levels.nvar <- u.nlevels[[re]][i.nvar[[re]]]
        logDetnvarGtotal[re] <- 0
        for (i in 1:n.nvar[re]) {
          logDetnvarGtotal[re] <- logDetnvarGtotal[re] + logDetnvarG[i] * (dim.nvar[re] / u.nlevels.nvar[[re]][i])
        }
      } else {
        logDetnvarGtotal[re] <- logDetnvarG
      }
      
      # Form the non-variable part of Gi using kronecker products before iterating if it exists:
      nvarGi[[re]] <- Reduce(kronecker, bases[[re]][i.nvar[[re]]])
    }
  } # End loop over random effects
  
  # Forming non-variable structure for the residual effect:
  re <- nUR
  if (n.nvar[re] > 0) { # This should always be true for the residual effect due to the id(units) component.
    # Form the non-variable part of R using kronecker products before iterating if it exists:
    S <- Reduce(kronecker, bases[[re]][i.nvar[[re]]])
  }
  
  ## Forming W = [X Z] for the dense and sparse equations ====
  Wdense <- cbind(X[, Xdense.idx, drop = FALSE], Reduce(cbind, Z[which(!(u.sparse[1:nU]))]))
  Wsparse <- cbind(X[, Xsparse.idx, drop = FALSE], Reduce(cbind, Z[which(u.sparse[1:nU])]))
  dense.idx <- c(Xdense.idx, unlist(Z.idx[Zdense.idx])) # This will come in handy for putting the solutions back into the order specified by the user.
  sparse.idx <- c(Xsparse.idx, unlist(Z.idx[Zsparse.idx])) # This will come in handy for putting the solutions back into the order specified by the user.
  dimC <- ncol(Wsparse) + ncol(Wdense) # Total number of equations.
  if (ncol(Wdense) == dimC) {
    Cmode <- "dense"
  } else if (ncol(Wsparse) == dimC) {
    Cmode <- "sparse"
  } else {
    Cmode <- "mixed"
  }
  
  stop <- proc.time()[3]
  duration <- as.numeric(stop - start)
  duration <- sprintf("%02d:%02d:%06.3f", duration %/% 3600, (duration %% 3600) %/% 60, duration %% 60)
  if (verbose) {
    ansi <- interactive() && !is.null(Sys.getenv("TERM"))
    if (Cmode == "mixed") {
      cat(sprintf("Solving %d dense equations and %d sparse equations %s GPU acceleration.\n", ncol(Wdense), ncol(Wsparse), ifelse(useGPU, "with", "without")))
    } else if (Cmode == "sparse") {
      cat(sprintf("Solving %d sparse equations.\n", ncol(Wsparse)))
    } else {
      cat(sprintf("Solving %d dense equations %s GPU acceleration.\n", ncol(Wdense), ifelse(useGPU, "with", "without")))
    }
    if (ansi) {
      cat("\033[1m", sprintf("%-11s %-15s %-20s", "Iteration", "LogLik", "Duration (hh:mm:ss)"), "\033[0m", "\n", sep = "")
    } else {
      cat(sprintf("%-11s %-15s %-20s \n", "Iteration", "LogLik", "Duration (hh:mm:ss)"))
    }
    cat(sprintf("%-11s %-15.3f %-20s\n", "Setup", NA, duration))
  }
  
  ## Setup the GPU CUDA context ====
  ##tic("Initializing CUDA")
  if (useGPU && Cmode %in% c("dense", "mixed")) {
    gpu_ctx <- cuSolverInit(N = ncol(Wdense), rhs_cols_1 = 1, rhs_cols_2 = nrow(k)) # N = dim of the dense block of C (Cdd) matrix, rhs_cols_2 is the number of parameters.
  } else if (useGPU && Cmode == "sparse") {
    stop("GPU cannot be used for a fully sparse model...")
  }
  ##toc()
  
  ## Start iterating ====
  LL <- numeric(maxit + 1)
  converged <- FALSE
  iter <- 1
  while (!(converged) && iter <= maxit) {
    start <- proc.time()[3]
    ### Computing Ri ====
    if (timing) {timer("start", "Computing Ri", iter)}
    # These lists will also be used later for the random effects:
    cormats <- vector("list", nUR)
    sigmai <- vector("list", nUR)
    logDetvarGtotal <- numeric(nUR)
    vars <- vector("list", nUR)
    re <- nUR # First do the very last structure, which is the residual effect.
    vars[[re]] <- tcrossprod(sqrt(k[indices[[re]][[1]], iter, drop = FALSE]))
    if (n.var[re] > 0) { # Only form the correlation matrices for variable factors if there are variable factors.
      cormats[[re]] <- vector("list", n.var[re]) # One correlation matrix for every variable factor.
      for (i in i.var[[re]]) { # Loop to form the correlations matrices
        # For the residual effect we are not going to compute the log determinants and inverses
        # of component matrices because we need to select columns/rows so it matches the data
        # anyway.
        result <- getCormat(struc = struc[[re]][i], # Which structure should be formed?
                            p = u.nlevels[[re]][i], # How many levels are there in this factor?
                            params = k[indices[[re]][[1 + i]], iter], # Which parameters should be used to form the correlation matrix?
                            base = bases[[re]][[i]], # Which, if any, base matrix should be used?
                            doLogDet = FALSE,
                            doCormati = FALSE) 
        cormats[[re]][[i]] <- result[[1]] # Correlation matrix
      } # End of loop to form the correlation matrices
      if (u.applyvars[[re]][1] == 0) { # If there is just a single variance.
        sigma <- vars[[re]]@x * Reduce(kronecker, cormats[[re]]) # Form the covariance matrix.
      } else if (all(u.applyvars[[re]] == i.var[[re]])) { # If variances are applied to all variable factors.
        sigma <- vars[[re]] * Reduce(kronecker, cormats[[re]]) # Form the inverse covariance matrix.
      } else { # If variances are applied to 1 or more factors, but not all factors.
        sigma <- kronecker(vars[[re]] * Reduce(kronecker, cormats[[re]][u.applyvars[[re]]]), Reduce(kronecker, cormats[[re]][-u.applyvars[[re]]])) # Form the covariance matrix.
      }
    } else { # if there are no variable factors.
      sigma <- vars[[re]]
    }
    # This should always be true for the residual effect as we always have an
    # identity matrix for the units:
    if (n.nvar[re] > 0) {
      R <- kronecker(sigma, S) # This is a bit slow as well.
      # This is where we multiple on the left and right with the residual effect
      # "design" matrices. These are actually just selection matrices that select
      # the appropriate rows/cols of R. This is necessary because the data may be
      # unbalanced: for a simple 3-trait GBLUP model, not all genotypes may have
      # data on all three traits. Some may only have data for a single trait or 2
      # of the 3 traits. This means that unlike the G matrices, the R matrix may
      # not be a clean Kronecker product. That means we have to explicitly form the
      # potentially messy R matrix before computing it's log determinant and inverse.
      # That's fine though, because the messy R will be incredibly sparse.
      # It's going to banded or maybe block-diagonal with potentially different sized
      # blocks.
      R <- forceSymmetric(tcrossprod(Z[[re]] %*% R, Z[[re]]))
      if (useWeights) {
        R <- weights %*% R %*% weights
      }
      if (isDiagonal(R)) {
        logDetR <- sum(log(R@x))
        Ri <- .symDiagonal(nrow(R), 1 / R@x)
      } else {
        Rchol <- Cholesky(R, LDL = FALSE, perm = TRUE, super = FALSE)
        logDetR <- determinant(Rchol, sqrt = FALSE, logarithm = TRUE)$modulus
        Ri <- solve(Rchol)
      }
    } else {
      R <- sigmai[[re]]
      stop("Something went wrong, the residual effect should always have a nonvariable part corresponding to id(units)...")
    }
    if (timing) {timer("stop")}
    
    ### Computing C = t(W) %*% Ri %*% W ====
    if (timing) {timer("start", "Forming base C", iter)}
    if (Cmode == "mixed") { # Both dense and sparse equations.
      Ri_Wdense <- Ri %*% Wdense
      Ri_Wsparse <- Ri %*% Wsparse
      Cdd <- as(forceSymmetric(crossprod(Wdense, Ri_Wdense)), "unpackedMatrix") # Top-left dense block of C.
      Css <- as(forceSymmetric(crossprod(Wsparse, Ri_Wsparse)), "TsparseMatrix") # Bottom-right sparse block of C.
      Cds <- crossprod(Wdense, Ri_Wsparse) # Top-right block of C.
      dimCss <- dim(Css)
      dimCdd <- dim(Cdd)
    } else if (Cmode == "dense") { # Only dense equations.
      Cdd <- as(forceSymmetric(crossprod(Wdense, Ri %*% Wdense)), "unpackedMatrix") # Full C!
    } else { # Only sparse equations.
      Css <- as(forceSymmetric(crossprod(Wsparse, Ri %*% Wsparse)), "TsparseMatrix") # Full C!
      dimCss <- dim(Css)
    }
    if (timing) {timer("stop")}
    
    ### Adding Gi to C ====
    if (timing) {timer("start", "Computing Gi and adding to base C", iter)}
    offsetDense <- nXdense + 1
    offsetSparse <- nXsparse + 1
    # if we have a sparse block of C, or if C is fully sparse.
    if (Cmode %in% c("mixed", "sparse")) {
      iList <- vector("list", nZsparse + 1)
      jList <- vector("list", nZsparse + 1)
      xList <- vector("list", nZsparse + 1)
      iList[[1]] <- Css@i
      jList[[1]] <- Css@j
      xList[[1]] <- Css@x
    }
    
    for (re in 1:nU) { # Start loop over random effects
      #### Computing Gi ====
      tic("Forming Sigmai")
      if (any(u.isStack[[re]])) { # If there is a stack in the random effect
        # Determining which factor is a stack and the size of the stack:
        i.stack <- which(u.isStack[[re]])
        sos <- length(u.nparams.list[[re]][[i.stack + 1]]) # + 1 due to the variances
        # Always the outer products of standard deviations for each stack first:
        vars[[re]] <- vector("list", sos)
        nvars <- length(indices[[re]][[1]]) / sos # Number of variances per stack layer
        for (s in 1:sos) {
          vars[[re]][[s]] <- tcrossprod(sqrt(k[indices[[re]][[1]][(s * nvars - nvars + 1):(s * nvars)], iter, drop = FALSE]))
        }
        if (n.var[re] > 0) { # This should always be the case for a random effect with a stack
          cormatis <- vector("list", n.var[re]) # One inverse correlation matrix for every variable factor.
          cormats[[re]] <- vector("list", n.var[re]) # One correlation matrix for every variable factor.
          logDetvarG <- numeric(n.var[re])
          for (i in i.var[[re]]) { # Loop to form the covariance stack and correlation matrices
            # The first variable factor has to be the one with a stack:
            if (i == i.stack) {
              # We need to calculate the complete weighed sum of the structures using the variances:
              covmat.stack <- Matrix(0, u.nlevels.var[[re]][i.var[[re]][i]], u.nlevels.var[[re]][i.var[[re]][i]], sparse = FALSE, doDiag = FALSE)
              cormats[[re]][[i]] <- vector("list", sos)
              indices.start <- 1
              for (s in 1:sos) {
                indices.stack <- (indices.start):(indices.start - 1 + u.nparams.list[[re]][[1 + i]][s])
                result <- getCormat(struc = struc[[re]][[i]][s], # Which structure should be formed?
                                    p = u.nlevels[[re]][i], # How many levels are there in this factor?
                                    params = k[indices[[re]][[1 + i]][indices.stack], iter], # Which parameters should be used to form the correlation matrix?
                                    base = bases[[re]][[i]][[s]], # Which, if any, base matrix should be used?
                                    doLogDet = FALSE, doCormati = FALSE) # Log determinant and inverses should not be computed for the individual layers of the stack.
                indices.start <- indices.start + u.nparams.list[[re]][[1 + i]][s]
                cormats[[re]][[i]][[s]] <- result[[1]] # Correlation matrix
                if (nvars == 1) {
                  covmat.stack <- covmat.stack + vars[[re]][[s]]@x * result[[1]] # We first construct sigma by summing all stack layers, then invert
                } else {
                  covmat.stack <- covmat.stack + vars[[re]][[s]] * result[[1]] # We first construct sigma by summing all stack layers, then invert
                }
              }
              cormatis[[i]] <- chol(covmat.stack, pivot = FALSE) # Now we compute the log determinant and inverse after completing the stack.
              logDetvarG[i] <- 2 * sum(log(diag(cormatis[[i]]))) # Log determinant.
              cormatis[[i]] <- chol2inv(cormatis[[i]]) # Inverse.
            } else {
              # Do the remaining variable structures that are not part of the stack (there can only be a single stack in a random effect):
              result <- getCormat(struc = struc[[re]][[i]], # Which structure should be formed?
                                  p = u.nlevels[[re]][i], # How many levels are there in this factor?
                                  params = k[indices[[re]][[1 + i]], iter], # Which parameters should be used to form the correlation matrix?
                                  base = bases[[re]][[i]]) # Which, if any, base matrix should be used?
              cormats[[re]][[i]] <- result[[1]] # Correlation matrix
              logDetvarG[i] <- result[[2]] # Log determinant of the correlation matrix
              cormatis[[i]] <- result[[3]] # Inverse of the correlation matrix
            }
          } # End of loop to form the covariance stack and variable correlation matrices
          # This is where we construct the full sigmai and logDetvarGtotal, which is easier for random effects with a stack:
          sigmai[[re]] <- Reduce(kronecker, cormatis) # The first entry of cormatis is actually the inverted stacked covariance matrix
          if (n.var[re] > 1) {
            logDetvarGtotal[re] <- 0
            for (i in i.var[[re]]) {
              logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * (dim.var[re] / u.nlevels.var[[re]][i])
            }
          } else {
            logDetvarGtotal[re] <- logDetvarG
          }
        } else { # End of if (n.var[re] > 0), which should always be TRUE for random effects with a stack
          stop("Something is wrong, n.var[re] > 0 is FALSE for a random effect with a stack...")
        }
      } else { # If there is no stack in the random effect
        # Always the outer product of standard deviations first:
        # For now, exfa is guaranteed/assumed to be the first factor in a random
        # effect, with variances applied to that factor, e.g., environments.
        if (any(struc[[re]] == "exfa")) {
          # 1 + 1, because the first entry are the variances, and exfa is guaranteed to be the second entry.
          # 1 for the u.nlevels as exfa is guaranteed to be the first structure.
          m <- u.nparams.list[[re]][[1 + 1]] / u.nlevels[[re]][1]
          vars[[re]] <- tcrossprod(rbind(sqrt(k[indices[[re]][[1]], iter, drop = FALSE]), Matrix(1, m, 1)))
        } else {
          m <- 0
          vars[[re]] <- tcrossprod(sqrt(k[indices[[re]][[1]], iter, drop = FALSE]))
        }
        if (n.var[re] > 0) { # Only form the correlation matrices for variable factors if there are variable factors.
          cormatis <- vector("list", n.var[re]) # One inverse correlation matrix for every variable factor.
          cormats[[re]] <- vector("list", n.var[re]) # One correlation matrix for every variable factor.
          logDetvarG <- numeric(n.var[re])
          for (i in i.var[[re]]) { # Loop to form the correlations matrices
            result <- getCormat(struc = struc[[re]][i], # Which structure should be formed?
                                p = u.nlevels[[re]][i], # How many levels are there in this factor?
                                params = k[indices[[re]][[1 + i]], iter], # Which parameters should be used to form the correlation matrix?
                                base = bases[[re]][[i]]) # Which, if any, base matrix should be used?
            cormats[[re]][[i]] <- result[[1]] # Correlation matrix
            logDetvarG[i] <- result[[2]] # Log determinant of the correlation matrix
            cormatis[[i]] <- result[[3]] # Inverse of the correlation matrix
          } # End of loop to form the correlation matrices
          if (u.applyvars[[re]][1] == 0) { # If there is just a single variance.
            sigmai[[re]] <- (1 / vars[[re]]@x) * Reduce(kronecker, cormatis) # Form the inverse covariance matrix.
            # Form the log determinant of sigma:
            if (n.var[re] > 1) {
              logDetvarGtotal[re] <- 0
              for (i in i.var[[re]]) {
                logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * (dim.var[re] / u.nlevels.var[[re]][i])
              }
            } else {
              logDetvarGtotal[re] <- logDetvarG
            }
            logDetvarGtotal[re] <- dim.var[re] * log(as.numeric(vars[[re]])) + logDetvarGtotal[re]
          } else if (all(u.applyvars[[re]] == i.var[[re]])) { # If variances are applied to all variable factors.
            sigmai[[re]] <- (1 / vars[[re]]) * Reduce(kronecker, cormatis) # Form the inverse covariance matrix.
            # Form the log determinant of sigma:
            logDetvarGtotal[re] <- sum(log(diag(vars[[re]])))
            for (i in i.var[[re]]) {
              logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * (dim.var[re] / u.nlevels.var[[re]][i])
            }
          } else { # If variances are applied to 1 or more factors, but not all factors.
            sigmai[[re]] <- kronecker((1 / vars[[re]]) * Reduce(kronecker, cormatis[u.applyvars[[re]]]), Reduce(kronecker, cormatis[-u.applyvars[[re]]])) # Form the inverse covariance matrix.
            # Form the log determinant of sigma:
            # First form the log determinant of the part with variances:
            # If there is an exfa structure, m is the number of levels and u.applyvars[[re]] is guaranteed to be 1.
            # If there is no exfa structure, m is 0 so nothing changes.
            u.nlevels.applyvar <- u.nlevels[[re]][u.applyvars[[re]]] + m # Numbers of levels for (variable, obviously) factors that have variances applied to them
            u.nlevels.napplyvar <- u.nlevels.var[[re]][-u.applyvars[[re]]] # Number of levels for variable factors that do not have variances applied to them
            logDetvarGtotal[re] <- sum(log(diag(vars[[re]])))
            for (i in u.applyvars[[re]]) { # for (i in 1) in the case of exfa.
              logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * prod(u.nlevels.applyvar[-i]) # prod(numeric(0)) = 1
            }
            logDetvarGtotal[re] <- logDetvarGtotal[re] * prod(u.nlevels.napplyvar)
            # Now we can form the full log determinant:
            if (any(struc[[re]] == "exfa")) {
              dim.var.corrected <- prod(c(u.nlevels.var[[re]][1] + m, u.nlevels.var[[re]][-1]))
            } else {
              dim.var.corrected <- dim.var[re]
            }
            for (i in i.var[[re]][-u.applyvars[[re]]]) {
              logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * (dim.var.corrected / u.nlevels.var[[re]][i])
            }
          }
        } else { # if there are no variable factors.
          sigmai[[re]] <- 1 / vars[[re]]
          logDetvarGtotal[re] <- -1 * log(sigmai[[re]])@x
        }
      }
      toc()
      ### Adding it to C ====
      if (n.nvar[re] > 0) {
        tic("Kronecker product to form Gi")
        Gi <- formGi(sigmai[[re]], nvarGi[[re]])
        toc()
      } else {
        Gi <- sigmai[[re]]
      }
      
      # Now add to C:
      tic("Adding to Ci")
      dimGi <- nrow(Gi)
      if (u.sparse[re]) { # Add to Css.
        if (!is(Gi, "TsparseMatrix")) { # Convert to triplet format to add to Css.
          Gi <- as(Gi, "TsparseMatrix")
        }
        # Collecting the values of the previous Css and current Gi.
        # Padding Gi with zeros by adding the offset to the coords to make it conformable with the sparse base Css.
        iList[[re + 1]] <- Gi@i + offsetSparse - 1 # - 1 is needed because i and j use 0-based indexing.
        jList[[re + 1]] <- Gi@j + offsetSparse - 1 # - 1 is needed because i and j use 0-based indexing.
        xList[[re + 1]] <- Gi@x
        offsetSparse <- offsetSparse + dimGi # Moving offset to next random effect.
        
      } else { # Add to Cdd.
        if (!is(Gi, "dsyMatrix")) { # Convert to dsyMatrix in case for some reason a sparse matrix is forced into the dense equations.
          Gi <- as(Gi, "unpackedMatrix")
        }
        addBlockInplace(Cdd, Gi, offset = offsetDense)
        offsetDense <- offsetDense + dimGi # Moving offset to next random effect.
      }
      toc()
    } # End loop over random effects.
    
    # If C is either fully or partly sparse, we can actually re-construct Css only once we have collected all
    # i's, j's, and x's for all sparse random effects.
    tic("Reconstructing Css from accumulated data and row/col indices")
    if (Cmode %in% c("mixed", "sparse")) {
      Css <- sparseMatrix(i = unlist(iList) + 1, j = unlist(jList) + 1, x = unlist(xList), dims = dimCss, symmetric = TRUE, repr = "C")
    }
    toc()
    if (timing) {timer("stop")}
    
    ### Cholesky decomposition of C and computing BLUEs/BLUPs ====
    if (timing) {timer("start", "CF of C and solves BLUEs/BLUPs", iter)}
    Riy <- Ri %*% y
    if (Cmode == "mixed") { # Mix of dense and sparse equations.
      Ccholss <- Cholesky(Css, LDL = ifelse(super, FALSE, TRUE), perm = TRUE, super = super)
      logDetCss <- as.numeric(determinant(Ccholss, sqrt = FALSE, logarithm = TRUE)$modulus)
      denseRHS <- crossprod(Wdense, Riy)
      sparseRHS <- crossprod(Wsparse, Riy)
      ##tic("Solving for BLUEs and BLUPs")
      Q <- solve(Ccholss, t(as(Cds, "unpackedMatrix")))
      SC <- Cdd - forceSymmetric(Cds %*% Q)
      sparseTemp <- solve(Ccholss, sparseRHS)
      adjDenseRHS <- denseRHS - (Cds %*% sparseTemp)
      if (!useGPU) {
        SCchol <- Cholesky(SC, perm = FALSE)
        # Looks strange because the 2 * is missing, but is correct as diag(Cholesky(SC, perm = FALSE))
        # returns the squared diagonal elements of the Cholesky factor, see ?Matrix::`Cholesky-class`
        logDetSC <- sum(log(diag(SCchol)))
        buDense <- solve(SCchol, adjDenseRHS)
        buSparse <- sparseTemp - (Q %*% buDense)
      } else {
        gpuResult <- cuSolverPhase1(ctx_ptr = gpu_ctx, C_in = SC, RHS_in = adjDenseRHS, returnCchol = FALSE)
        logDetSC <- gpuResult$logDetC
        buDense <- gpuResult$bu
        rownames(buDense) <- colnames(SC)
        buSparse <- sparseTemp - (Q %*% buDense)
      }
      ##toc()
      # logDetC <- as.numeric(determinant(Ccholss, sqrt = FALSE, logarithm = TRUE)$modulus + determinant(SC, logarithm = TRUE)$modulus)
      logDetC <- logDetCss + logDetSC
      yPy <- as.numeric(crossprod(y, Riy) - (crossprod(sparseRHS, buSparse) + crossprod(denseRHS, buDense)))
      
    } else if (Cmode == "sparse") { # Fully sparse equations.
      Ccholss <- Cholesky(Css, LDL = ifelse(super, FALSE, TRUE), perm = TRUE, super = super)
      RHS <- crossprod(Wsparse, Riy)
      ##tic("Solving for BLUEs and BLUPs")
      buSparse <- solve(Ccholss, RHS)
      ##toc()
      logDetC <- determinant(Ccholss, sqrt = FALSE, logarithm = TRUE)$modulus
      yPy <- as.numeric(crossprod(y, Riy) - crossprod(RHS, buSparse))
      
    } else { # Fully dense equations.
      RHS <- crossprod(Wdense, Riy)
      if (!useGPU) {
        Ccholdd <- Cholesky(Cdd, perm = FALSE)
        ##tic("Solving for BLUEs and BLUPs")
        buDense <- solve(Ccholdd, RHS)
        ##toc()
        # Looks strange because the 2 * is missing, but is correct as diag(Cholesky(Cdd, perm = FALSE))
        # returns the squared diagonal elements of the Cholesky factor, see ?Matrix::`Cholesky-class`
        logDetC <- sum(log(diag(Ccholdd)))
      } else {
        gpuResult <- cuSolverPhase1(ctx_ptr = gpu_ctx, C_in = Cdd, RHS_in = RHS, returnCchol = FALSE)
        buDense <- gpuResult$bu
        rownames(buDense) <- colnames(Cdd)
        logDetC <- gpuResult$logDetC
      }
      yPy <- as.numeric(crossprod(y, Riy) - crossprod(RHS, buDense))
    }
    
    # Splitting BLUEs and BLUPs:
    if (Cmode == "mixed") {
      if (nXdense > 0) { bDense <- buDense[1:nXdense, 1, drop = FALSE] }
      if (nZdense > 0) { uDense <- buDense[(nXdense + 1):nrow(buDense), 1, drop = FALSE] }
      if (nXsparse > 0) { bSparse <- buSparse[1:nXsparse, 1, drop = FALSE] }
      if (nZsparse > 0) { uSparse <- buSparse[(nXsparse + 1):nrow(buSparse), 1, drop = FALSE] }
    } else if (Cmode == "sparse") {
      bSparse <- buSparse[1:nXsparse, 1, drop = FALSE]
      uSparse <- buSparse[(nXsparse + 1):nrow(buSparse), 1, drop = FALSE]
    } else {
      bDense <- buDense[1:nXdense, 1, drop = FALSE]
      uDense <- buDense[(nXdense + 1):nrow(buDense), 1, drop = FALSE]
    }
    if (timing) {timer("stop")}
    
    ### Computing the log-likelihood ====
    if (timing) {timer("start", "Log-likelihood", iter)}
    # logDetG:
    logDetG <- 0
    for (re in 1:nU) { # Start loop over random effects
      if (any(struc[[re]] == "exfa")) {
        m <- u.nparams.list[[re]][[1 + 1]] / u.nlevels[[re]][1]
        dim.var.corrected <- prod(c(u.nlevels.var[[re]][1] + m, u.nlevels.var[[re]][-1]))
      } else {
        dim.var.corrected <- dim.var[re]
      }
      if (n.nvar[re] > 0) {
        logDetG <- logDetG + dim.nvar[re] * logDetvarGtotal[re] + dim.var.corrected * logDetnvarGtotal[re]
      } else {
        logDetG <- logDetG + dim.nvar[re] * logDetvarGtotal[re]
      }
    } # End loop over random effects
    
    # Full LL:
    LL[iter] <- -0.5 * (logDetC + logDetR + logDetG + yPy)
    if (timing) {timer("stop")}
    
    ### Working variates ====
    #### Genetic working variates ====
    if (timing) {timer("start", "Random effect working variates", iter)}
    Wu <- Matrix(0, n, nrow(k), sparse = FALSE, doDiag = FALSE) # Matrix for the working variates.
    umat_sigmai <- vector("list", nU) # We will store these matrices for computing the scores later.
    derivs <- vector("list", nUR) # Same thing for all derivatives, except the last entry is reserved for the residual effect.
    uDense.start <- 1
    uSparse.start <- 1
    param <- 1 # Keeping track of which column/parameter of Wu we are working on.
    for (re in 1:nU) { # Start of loop over random effects
      if (any(struc[[re]] == "exfa")) {
        m <- u.nparams.list[[re]][[1 + 1]] / u.nlevels[[re]][1]
        dim.var.corrected <- prod(c(u.nlevels.var[[re]][1] + m, u.nlevels.var[[re]][-1]))
      } else {
        m <- 0
        dim.var.corrected <- dim.var[re]
      }
      dim.total <- dim.var.corrected * dim.nvar[re]
      if (u.sparse[re]) {
        u <- uSparse[uSparse.start:(uSparse.start - 1 + dim.total), 1]
        uSparse.start <- uSparse.start + dim.total
      } else {
        u <- uDense[uDense.start:(uDense.start - 1 + dim.total), 1]
        uDense.start <- uDense.start + dim.total
      }
      umat <- Matrix(u, nrow = ifelse(n.nvar[re] > 0, dim.nvar[re], 1), ncol = dim.var.corrected)
      umat_sigmai[[re]] <- umat %*% sigmai[[re]]
      derivs[[re]] <- Matrix(0, dim.var.corrected, dim.var.corrected * u.nparams[[re]], sparse = FALSE, doDiag = FALSE) # Allocate big matrix for batched derivatives.
      if (any(u.isStack[[re]])) { # If there is a stack in this random effect
        # Determining which factor is a stack and the size of the stack:
        i.stack <- which(u.isStack[[re]])
        sos <- length(u.nparams.list[[re]][[i.stack + 1]]) # + 1 due to the variances
        nvars <- length(indices[[re]][[1]]) / sos # Number of variances per stack layer
        for (i in 1:u.nparams[[re]]) {
          col.start <- (i - 1) * dim.var.corrected + 1 # Start col index of the matrix containing batched derivatives.
          col.end <- i * dim.var.corrected # End col index of the matrix containing batched derivatives.
          idx.start <- (col.start - 1) * dim.var.corrected + 1
          idx.end <- col.end * dim.var.corrected
          if (isFixed[param]) {
            param <- param + 1
          } else {
            if (ur.types[param] == "V") {
              # We need to figure out for which stack layer the current variance is:
              s <- ceiling(i / nvars)
              if (u.applyvars[[re]][1] == 0) { # If there is only a single variance.
                derivs[[re]]@x[idx.start:idx.end] <- Reduce(kronecker, c(cormats[[re]][[i.stack]][s], cormats[[re]][-i.stack]))@x # The stack is always the first factor.
              } else { # If variances are applied to the factor in the stack.
                # In this case u.applyvars[[re]] is always equal to 1.
                A <- Matrix(0, prod(u.nlevels[[re]][u.applyvars[[re]]]), prod(u.nlevels[[re]][u.applyvars[[re]]]), sparse = FALSE, doDiag = FALSE)
                A[i - (s - 1) * nvars,] <- 0.5 # we can't use i directly
                if (n.var[re] > 1) { # If there are more variable factors besides the stack.
                  derivs[[re]]@x[idx.start:idx.end] <- kronecker((A + t(A)) * vars[[re]][[s]] * cormats[[re]][[u.applyvars[[re]]]][[s]] / k[param, iter], Reduce(kronecker, cormats[[re]][-u.applyvars[[re]]]))@x
                } else { # If the stack is the only variable factor.
                  derivs[[re]]@x[idx.start:idx.end] <- ((A + t(A)) * vars[[re]][[s]] * cormats[[re]][[u.applyvars[[re]]]][[s]] / k[param, iter])@x
                }
              }
            } else {
              j <- u.forFactor[[re]][i] # Which factor does this parameter belong to?
              cormats1 <- cormats[[re]] # Making a copy of the correlation matrices so we can replace the matrix of factor j with the derivative.
              if (j == i.stack) { # If we are looking at a parameter that is part of a structure in the stack
                # We need to determine which layer of the stack this parameter appears in:
                s <- findInterval(i - nvars * sos, c(0, cumsum(u.nparams.list[[re]][[j + 1]])), left.open = TRUE)
                indices.stack.start <- 1 + c(0, cumsum(u.nparams.list[[re]][[j + 1]]))[s]
                indices.stack.end <- indices.stack.start - 1 + u.nparams.list[[re]][[j + 1]][s]
                indices.stack <- indices.stack.start:indices.stack.end # Indices of the parameters that belong to this layer of the stack.
                cormats1[[j]] <- getDeriv(struc = struc[[re]][[j]][s],
                                          type = ur.types[param],
                                          cormatIndex = u.indexCor[[re]][[i]],
                                          params = k[indices[[re]][[1 + j]][indices.stack], iter],
                                          p = u.nlevels[[re]][j],
                                          base = bases[[re]][[j]][[s]],
                                          inverse = FALSE)
                # Now add the relevant variances to compute the final derivative.
                # We're dealing with a stack so there can either be 1 variance per layer or a variance for each level of the stacked factor:
                if (u.applyvars[[re]][1] == 0) { # If there is a single variance per layer.
                  derivs[[re]]@x[idx.start:idx.end] <- (vars[[re]][[s]]@x * Reduce(kronecker, cormats1))@x
                } else { # If there is a variance for each level of the stacked factor.
                  derivs[[re]]@x[idx.start:idx.end] <- Reduce(kronecker, c(list(vars[[re]][[s]] * cormats1[[i.stack]]), cormats1[-i.stack]))@x
                }
              } else { # If we are looking at a parameter for a structure that is not part of the stack
                cormats1[[j]] <- getDeriv(struc = struc[[re]][[j]],
                                          type = ur.types[param],
                                          cormatIndex = u.indexCor[[re]][[i]],
                                          params = k[indices[[re]][[1 + j]], iter],
                                          p = u.nlevels[[re]][j],
                                          base = bases[[re]][[j]],
                                          inverse = FALSE)
                # Now add the relevant variances to compute the final derivative.
                # We're dealing with a stack so there can either be 1 variance per layer or a variance for each level of the stacked factor:
                if (u.applyvars[[re]][1] == 0) { # If there is a single variance per layer.
                  covmat.stack <- Reduce(`+`, lapply(1:sos, function(s) {vars[[re]][[s]]@x * cormats1[[i.stack]][[s]]}))
                } else { # If there is a variance for each level of the stacked factor.
                  covmat.stack <- Reduce(`+`, lapply(1:sos, function(s) {vars[[re]][[s]] * cormats1[[i.stack]][[s]]}))
                }
                derivs[[re]]@x[idx.start:idx.end] <- kronecker(covmat.stack, Reduce(kronecker, cormats1[-i.stack]))@x
              }
            }
            param <- param + 1
          }
        } # End of loop over the parameters for this random effect.
      } else { # If there is no stack
        for (i in 1:u.nparams[[re]]) {
          col.start <- (i - 1) * dim.var.corrected + 1
          col.end <- i * dim.var.corrected
          idx.start <- (col.start - 1) * dim.var.corrected + 1
          idx.end <- col.end * dim.var.corrected
          if (isFixed[param]) {
            param <- param + 1
          } else {
            if (ur.types[param] == "V") { # If the parameter is a variance.
              if (u.applyvars[[re]][1] == 0) { # If there is only a single variance.
                if (n.var[re] > 0) {
                  derivs[[re]]@x[idx.start:idx.end] <- Reduce(kronecker, cormats[[re]])@x
                } else {
                  derivs[[re]]@x[idx.start:idx.end] <- Matrix(1, sparse = FALSE, doDiag = FALSE)@x
                }
              } else { # If variances are applied to 1 or more factors.
                # If there is no exfa structure, m = 0 and nothing changes.
                # If there is an exfa structure, u.applyvars[[re]] = 1 by definition and m is added to the number of levels of the factor modeled by the exfa structure.
                A <- Matrix(0, prod(u.nlevels[[re]][u.applyvars[[re]]]) + m, prod(u.nlevels[[re]][u.applyvars[[re]]]) + m, sparse = FALSE, doDiag = FALSE)
                A[i,] <- 0.5
                if (all(u.applyvars[[re]] == i.var[[re]])) { # If variances are applied to all variable factors.
                  derivs[[re]]@x[idx.start:idx.end] <- ((A + t(A)) * vars[[re]] * Reduce(kronecker, cormats[[re]]) / k[param, iter])@x
                } else { # If variances are applied to 1 or more factors, but not all factors.
                  derivs[[re]]@x[idx.start:idx.end] <- kronecker((A + t(A)) * vars[[re]] * Reduce(kronecker, cormats[[re]][u.applyvars[[re]]]) / k[param, iter], Reduce(kronecker, cormats[[re]][-u.applyvars[[re]]]))@x
                }
              }
            } else { # If the parameter is not a variance.
              j <- u.forFactor[[re]][i] # Which factor does this parameter belong to?
              cormats1 <- cormats[[re]] # Making a copy of the correlation matrices so we can replace the matrix of factor j with the derivative.
              cormats1[[j]] <- getDeriv(struc = struc[[re]][j],
                                        type = ur.types[param],
                                        cormatIndex = u.indexCor[[re]][[i]],
                                        params = k[indices[[re]][[1 + j]], iter],
                                        p = u.nlevels[[re]][j],
                                        base = bases[[re]][[j]],
                                        inverse = FALSE)
              if (u.applyvars[[re]][1] == 0) {
                derivs[[re]]@x[idx.start:idx.end] <- (vars[[re]]@x * Reduce(kronecker, cormats1))@x
              } else if (all(u.applyvars[[re]] == i.var[[re]])) {
                derivs[[re]]@x[idx.start:idx.end] <- (vars[[re]] * Reduce(kronecker, cormats1))@x
              } else {
                derivs[[re]]@x[idx.start:idx.end] <- kronecker(vars[[re]] * Reduce(kronecker, cormats1[u.applyvars[[re]]]), Reduce(kronecker, cormats1[-u.applyvars[[re]]]))@x
              }
            }
            param <- param + 1
          }
        }
      }
      # We have now collected all derivatives in a big matrix ready for the batched execution.
      # Get all working variates for this random effect at once:
      col.start <- param - u.nparams[[re]]
      col.end <- param - 1
      USD <- batchedWorkingVariates(umat_sigmai[[re]], derivs[[re]]) # Batched execution and reshaping in C++.
      Wu[, col.start:col.end] <- Z[[re]] %*% USD # Single multiplication and dropping into Wu.
    } # End of loop over random effects
    if (timing) {timer("stop")}
    
    #### Residual working variates ====
    if (timing) {timer("start", "Residual working variates", iter)}
    if (Cmode == "mixed") { # If there are both dense and sparse equations.
      e <- y - Wdense %*% buDense - Wsparse %*% buSparse
    } else if (Cmode == "sparse") { # If there are only sparse equations.
      e <- y - Wsparse %*% buSparse
    } else { # If there are only dense equations.
      e <- y - Wdense %*% buDense
    }
    re <- nUR
    derivs[[re]] <- vector("list", u.nparams[[re]])
    for (i in 1:u.nparams[[re]]) {
      if (isFixed[param]) {
        param <- param + 1
      } else {
        if (ur.types[param] == "V") { # If the parameter is a variance.
          if (u.applyvars[[re]][1] == 0) { # If there is only a single variance.
            # Don't forget to multiply with the non-variable structure and selection matrices:
            derivs[[re]][[i]] <- forceSymmetric(tcrossprod(Z[[re]] %*% Reduce(kronecker, c(cormats[[re]], S)), Z[[re]]))
          } else { # If variances are applied to 1 or more factors.
            A <- Matrix(0, prod(u.nlevels[[re]][u.applyvars[[re]]]), prod(u.nlevels[[re]][u.applyvars[[re]]]), sparse = FALSE, doDiag = FALSE)
            A[i,] <- 0.5
            if (all(u.applyvars[[re]] == i.var[[re]])) { # If variances are applied to all variable factors.
              derivs[[re]][[i]] <- forceSymmetric(tcrossprod(Z[[re]] %*% kronecker((A + t(A)) * vars[[re]] * Reduce(kronecker, cormats[[re]]) / k[param, iter], S), Z[[re]]))
            } else { # If variances are applied to 1 or more factors, but not all factors.
              derivs[[re]][[i]] <- forceSymmetric(tcrossprod(Z[[re]] %*% kronecker(kronecker((A + t(A)) * vars[[re]] * Reduce(kronecker, cormats[[re]][u.applyvars[[re]]]) / k[param, iter], Reduce(kronecker, cormats[[re]][-u.applyvars[[re]]])), S), Z[[re]]))
            }
          }
        } else { # If the parameter is not a variance.
          j <- u.forFactor[[re]][i] # Which factor does this parameter belong to?
          cormats1 <- cormats[[re]] # Making a copy of the correlation matrices so we can replace the matrix of factor j with the derivative.
          cormats1[[j]] <- getDeriv(struc = struc[[re]][j],
                                    type = ur.types[param],
                                    cormatIndex = u.indexCor[[re]][[i]],
                                    params = k[indices[[re]][[1 + j]], iter],
                                    p = u.nlevels[[re]][j],
                                    base = bases[[re]][[j]],
                                    inverse = FALSE)
          if (u.applyvars[[re]][1] == 0) {
            derivs[[re]][[i]] <- forceSymmetric(tcrossprod(Z[[re]] %*% (vars[[re]]@x * Reduce(kronecker, c(cormats1, S))), Z[[re]]))
          } else if (all(u.applyvars[[re]] == i.var[[re]])) {
            derivs[[re]][[i]] <- forceSymmetric(tcrossprod(Z[[re]] %*% kronecker(vars[[re]] * Reduce(kronecker, cormats1), S), Z[[re]]))
          } else {
            derivs[[re]][[i]] <- forceSymmetric(tcrossprod(Z[[re]] %*% kronecker(kronecker(vars[[re]] * Reduce(kronecker, cormats1[u.applyvars[[re]]]), Reduce(kronecker, cormats1[-u.applyvars[[re]]])), S), Z[[re]]))
          }
        }
        if (useWeights) {
          derivs[[re]][[i]] <- weights %*% derivs[[re]][[i]] %*% weights
        }
        Wu[, param] <- derivs[[re]][[i]] %*% Ri %*% e
        param <- param + 1
      }
    }
    if (timing) {timer("stop")}
    
    ### Computing the AI matrix ====
    if (timing) {timer("start", "AI matrix and Cidd if on GPU", iter)}
    RiWu <- Ri %*% Wu
    WiWj <- forceSymmetric(crossprod(Wu, RiWu)) # Similar to t(y) %*% Ri %*% y
    if (Cmode == "mixed") { # If we have both dense and sparse equations.
      denseWuRHS <- crossprod(Wdense, RiWu)
      sparseWuRHS <- crossprod(Wsparse, RiWu)
      sparseWuTemp <- solve(Ccholss, sparseWuRHS)
      adjDenseWuRHS <- denseWuRHS - (Cds %*% sparseWuTemp)
      if (!useGPU) {
        buWuDense <- solve(SCchol, adjDenseWuRHS)
      } else {
        gpuResult2 <- cuSolverPhase2(ctx_ptr = gpu_ctx, XWjxZWj_in = adjDenseWuRHS)
        buWuDense <- gpuResult2$buWu
      }
      buWuSparse <- sparseWuTemp - (Q %*% buWuDense)
      buWuXWjxZWj <- forceSymmetric(crossprod(sparseWuRHS, buWuSparse) + crossprod(denseWuRHS, buWuDense))
      
    } else if (Cmode == "sparse") { # If we only have sparse equations.
      XWjxZWj <- solve(Ccholss, crossprod(Wsparse, RiWu), system = "P") # Similar to permRHS
      zWu <- solve(Ccholss, XWjxZWj, system = "L") # Similar to z
      buWu <- solve(Ccholss, zWu, system = "Lt") # Similar to obtaining bu using a backsolve
      buWuXWjxZWj <- forceSymmetric(crossprod(buWu, XWjxZWj))
      
    } else { # If we only have dense equations.
      XWjxZWj <- crossprod(Wdense, RiWu) # t(W) %*% Ri %*% Wu
      if (!useGPU) {
        buWu <- solve(Ccholdd, XWjxZWj)
      } else {
        gpuResult2 <- cuSolverPhase2(ctx_ptr = gpu_ctx, XWjxZWj_in = XWjxZWj)
        buWu <- gpuResult2$buWu
      }
      buWuXWjxZWj <- forceSymmetric(crossprod(buWu, XWjxZWj))
    }
    AI <- (WiWj - buWuXWjxZWj) / 2
    
    if (avInfReg[iter] > 1) { # Potential inflation of diagonal elements of the AI matrix
      # diag(AI) <- diag(AI) * avInfReg[iter]
      diag(AI)[ur.types == "L"] <- diag(AI)[ur.types == "L"] * avInfReg[iter] # Alternative that only inflates loadings.
    }
    if (any(isFixed)) { # If there are any fixed parameters, remove them from the AI matrix before inverting.
      AIi <- Matrix(0, nrow(AI), ncol(AI), sparse = FALSE, doDiag = FALSE)
      notFixed <- which(!isFixed)
      AIi[notFixed, notFixed] <- solve(AI[notFixed, notFixed])
    } else {
      AIi <- solve(AI)
    }
    if (timing) {timer("stop")}
    
    ### Computing the inverse of the coefficient matrix ====
    if (timing) {timer("start", "Inversion of C", iter)}
    if (Cmode == "mixed") { # If we have both dense and sparse equations.
      #tic("Base SIS of Ciss")
      Ciss <- SIS(Ccholss) # SIS() returns a lower triangular dsCMatrix!!!
      #toc()
      # Inverse of the dense Schur complement:
      if (!useGPU) {
        #tic("Inversion of SC")
        Cidd <- chol2inv(as(SCchol, "dtrMatrix"))
        #toc()
      } else {
        Cidd <- gpuResult2$Ci
      }
      # Updated sparse block of the inverse of C:
      #tic("Updating the SIS")
      inplaceUpdateSIS(Ciss, tcrossprod(Cidd, Q), t(Q))
      #toc()
      # Off-diagonal (bottom left) block of the inverse of C:
      #tic("Off-diagonal block of Ci")
      Cisd <- -Q %*% Cidd
      #toc()
      
    } else if (Cmode == "sparse") { # If we only have sparse equations, so Ciss is the full inverse.
      Ciss <- SIS(Ccholss) # SIS() returns a lower triangular dsCMatrix!!!
      
    } else { # if we only have dense equations, so Cidd is the full inverse.
      if (!useGPU) {
        Cidd <- chol2inv(as(Ccholdd, "dtrMatrix"))
      } else {
        Cidd <- gpuResult2$Ci
      }
    }
    if (timing) {timer("stop")}
    
    ### Computing scores ====
    if (timing) {timer("start", "Random effect scores", iter)}
    # See equations (18) and (15) from Jensen, Mantysaari, Madsen, & Thompson 1997:
    # "Residual maximum likelihood estimation of (co) variance components in multivariate
    # mixed linear models using average information"
    #### Genetic scores ====
    if (debugging) {
      traces.debugging <- vector("list", nU)
    }
    scores <- Matrix(0, nrow(k), 1, sparse = FALSE)
    offsetDense <- nXdense + 1
    offsetSparse <- nXsparse + 1
    param <- 1
    for (re in 1:nU) { # Start of loop over random effects
      if (any(struc[[re]] == "exfa")) {
        m <- u.nparams.list[[re]][[1 + 1]] / u.nlevels[[re]][1]
        dimVar <- prod(c(u.nlevels.var[[re]][1] + m, u.nlevels.var[[re]][-1]))
      } else {
        dimVar <- dim.var[re]
      }
      dimNvar <- dim.nvar[re]
      
      if (n.nvar[re] > 0) { # If there is a non-variable structure like a GRM.
        if (u.sparse[re]) { # If the current random effect is part of the sparse equations.
          if (is(nvarGi[[re]], "sparseMatrix")) { # If the non-variable structure is sparse.
            traces <- getBlockTracesSparseSparse(Ciss@p, Ciss@i, Ciss@x, nvarGi[[re]]@p, nvarGi[[re]]@i, nvarGi[[re]]@x, dimVar, dimNvar, offsetSparse - 1)
          } else { # If the non-variable structure is dense.
            traces <- getBlockTracesSparseDense(Ciss@p, Ciss@i, Ciss@x, nvarGi[[re]]@x, dimVar, dimNvar, offsetSparse - 1)
          }
          offsetSparse <- offsetSparse + dimVar * dimNvar
          
        } else { # If the current random effect is part of the dense equations.
          if (is(nvarGi[[re]], "sparseMatrix")) { # If the non-variable structure is sparse, I don't think this will typically ever be the case cause if true, we'd probably be in the sparse equations.
            traces <- getBlockTracesDenseSparse(Cidd@x, nrow(Cidd), nvarGi[[re]]@p, nvarGi[[re]]@i, nvarGi[[re]]@x, dimVar, dimNvar, offsetDense - 1)
          } else { # If the non-variable structure is dense.
            traces <- getBlockTracesDenseDense(Cidd@x, nrow(Cidd), nvarGi[[re]]@x, dimVar, dimNvar, offsetDense - 1)
          }
          offsetDense <- offsetDense + dimVar * dimNvar
        }
        
        if (nrow(traces) > 1) {
          traces <- forceSymmetric(traces)
        }
        if (debugging) {
          traces.debugging[[re]] <- traces
        }
        
      } else { # If there is no non-variable structure like a GRM.
        if (u.sparse[re]) { # If part of the sparse equations.
          traces <- Ciss[offsetSparse:(offsetSparse - 1 + dimVar), offsetSparse:(offsetSparse - 1 + dimVar)]
          offsetSparse <- offsetSparse + dimVar
        } else { # If part of the dense equations.
          traces <- Cidd[offsetSparse:(offsetSparse - 1 + dimVar), offsetSparse:(offsetSparse - 1 + dimVar)]
          offsetDense <- offsetDense + dimVar
        }
      }
      
      if (n.nvar[re] > 0) {
        basisScores <- -0.5 * (dimNvar * sigmai[[re]] - forceSymmetric(sigmai[[re]] %*% traces %*% sigmai[[re]]) - forceSymmetric(crossprod(umat_sigmai[[re]], nvarGi[[re]] %*% umat_sigmai[[re]]))) # Last term maybe on the GPU?
      } else {
        basisScores <- -0.5 * (sigmai[[re]] - forceSymmetric(sigmai[[re]] %*% traces %*% sigmai[[re]]) - crossprod(umat_sigmai[[re]])) # Last term maybe on the GPU?
      }
      
      # Old unbatched version:
      # for (i in 1:length(derivs[[re]])) {
      #   if (!isFixed[param]) {
      #     scores[param, 1] <- sum(derivs[[re]][[i]] * basisScores)
      #   }
      #   param <- param + 1
      # }
      
      # New batched version:
      row.start <- param
      row.end <- param + u.nparams[[re]] - 1
      scores[row.start:row.end, 1] <- batchedScores(basisScores, derivs[[re]])
      param <- param + u.nparams[[re]] # Move to the start of the set of parameters for the next random effect.
    } # End of loop over random effects.
    if (timing) {timer("stop")}
    
    #### Residual scores ====
    if (timing) {timer("start", "Residual scores", iter)}
    re <- nUR
    Rie <- Ri %*% e
    if (Cmode == "mixed") {
      RiWdense <- Ri %*% Wdense
      RiWsparse <- Ri %*% Wsparse
    } else if (Cmode == "sparse") {
      RiWsparse <- Ri %*% Wsparse
    } else {
      RiWdense <- Ri %*% Wdense
    }
    
    for (i in 1:length(derivs[[re]])) {
      if (!isFixed[param]) {
        dR <- derivs[[re]][[i]]
        RidR <- Ri %*% dR
        tr1 <- sum(diag(RidR))
        quad <- as.numeric(crossprod(Rie, dR %*% Rie))
        if (Cmode == "mixed") {
          RidRRi_Wdense  <- RidR %*% RiWdense
          RidRRi_Wsparse <- RidR %*% RiWsparse
          Mdd <- forceSymmetric(crossprod(Wdense, RidRRi_Wdense))
          Msd <- crossprod(Wsparse, RidRRi_Wdense)
          Mss <- forceSymmetric(crossprod(Wsparse, RidRRi_Wsparse))
          ddTrace <- sum(Mdd * Cidd)
          # ssTrace <- sum(diag(Mss %*% Ciss))
          if (!is(Mss, "dsCMatrix")) {
            stop("Mss must be a dsCMatrix...")
          }
          ssTrace <- as.numeric(getBlockTracesSparseSparse(Mss@p, Mss@i, Mss@x,
                                                           Ciss@p, Ciss@i, Ciss@x, dimVar = 1, dimNvar = nrow(Ciss), offset = 0))
          sdTrace <- 2 * sum(Msd * Cisd)
          tr2 <- ssTrace + sdTrace + ddTrace
          
        } else if (Cmode == "sparse") {
          Wsparse_RidRRi_Wsparse <- forceSymmetric(crossprod(RiWsparse, dR %*% RiWsparse))
          # tr2 <- sum(diag(Wsparse_RidRRi_Wsparse %*% Ciss)) # As both Wsparse_RidRRi_Wsparse and Ciss are sparse, this is probably the fastest.
          if (!is(Wsparse_RidRRi_Wsparse, "dsCMatrix")) {
            stop("Wsparse_RidRRi_Wsparse must be a dsCMatrix...")
          }
          tr2 <- as.numeric(getBlockTracesSparseSparse(Wsparse_RidRRi_Wsparse@p, Wsparse_RidRRi_Wsparse@i, Wsparse_RidRRi_Wsparse@x,
                                                       Ciss@p, Ciss@i, Ciss@x, dimVar = 1, dimNvar = nrow(Ciss), offset = 0))
          
        } else {
          Wdense_RidRRi_Wdense <- forceSymmetric(crossprod(RiWdense, dR %*% RiWdense))
          tr2 <- getTraceSparseDense(p = Wdense_RidRRi_Wdense@p, i = Wdense_RidRRi_Wdense@i, x = Wdense_RidRRi_Wdense@x, dense_x = Cidd@x, n_rows = nrow(Cidd))
          
        }
        scores[param, 1] <- -0.5 * (tr1 - tr2 - quad)
      }
      param <- param + 1
    }
    if (timing) {timer("stop")}
    
    ### Updating parameters ====
    if (timing) {timer("start", "Updating parameters and constraining", iter)}
    delta <- AIi %*% scores * stepSize[iter]
    #### Checking constraints ====
    # - We check all parameters and constrain as needed (adding them to con).
    # - We then recompute delta for any unconstrained parameters.
    # - We then check again to see if any of the previously unconstrained parameters
    #   now need to be constrained.
    # - If so, we constrain them and add them to con.
    # - Then we recompute delta again for any parameters that remained unconstrained.
    # - ...
    # - We exit this loop when:
    #     (1) none of the unconstrained parameters have to be constrained (this can happen in the first cycle).
    #     Or
    #     (2) all of the parameters are constrained (this can also happen in the first cycle).
    if (debugging) {
      vars.debugging <- vars
    }
    con <- c() # Vector to keep the total set of constrained parameters across cycles. We need this for recomputing delta for unconstrained params.
    vars <- which(ur.types == "V") # Which of the parameters in k are variances?
    cors <- which(ur.types == "R") # Which of the parameters in k are correlations?
    autoCors <- which(ur.types == "AR") # Which of the parameters in k are autocorrelations?
    bws <- which(ur.types == "B") # Which of the parameters in k are bandwidths?
    checkConstraints <- TRUE # Control variable to break the loop after satisfying (1) or (2).
    constraintsUsed <- FALSE # Variable that lets us print to the console if any params were constrained.
    q <- 1 # If we can't get a good update after 10 cycles we stop.
    while (checkConstraints) {
      if (q > 10) {
        stop("Cannot constrain to get a good update...")
      }
      checkConstraints <- FALSE # Assume everything is fine and we don't need to check again.
      
      constrained <- FALSE # Were any parameters constrained because the deltas were too large?
      hitBoundary <- FALSE # Were any parameters constrained because the delta pushed them across a boundary?
      nonPD <- FALSE # Were any correlations constrained because they resulted in a non-PD corg matrix?
      con.vars <- c() # Vector with indices of variances that were constrained in the current cycle due to the deltas being too large.
      con.cors <- c() # Vector with indices of correlations that were constrained in the current cycle due to the deltas being too large.
      con.autoCors <- c() # Vector with indices of autocorrelations that were constrained in the current cycle due to the deltas being too large.
      con.bws <- c() # Vector with indices of bandwidths that were constrained in the current cycle due to the deltas being too large.
      con.loadings <- c() # Vector with indices of loadings that were constrained in the current cycle due to the deltas being too large.
      bound.vars <- c() # Vector with indices of variances that were constrained in the current cycle due to the deltas pushing them across boundaries.
      bound.cors <- c() # Vector with indices of correlations that were constrained in the current cycle due to the deltas pushing them across boundaries.
      bound.autoCors <- c() # Vector with indices of correlations that were constrained in the current cycle due to the deltas pushing them across boundaries.
      bound.bws <- c() # Vector with indices of bandwidths that were constrained in the current cycle due to the deltas pushing them across boundaries.
      bound.loadings <- c() # Vector with indices of loadings that were constrained in the current cycle due to the deltas pushing the communalities over 1 - 1e-3.
      npd.cors <- c() # Vector with indices of correlations that were constrained in the current cycle due to the presence of non-PD corg matrices.
      
      # Computing the proposed update so we can check for constraints:
      kUpdate <- k[, iter, drop = FALSE] + delta
      
      # Checking variances for changes that result in an increase larger than a factor 15 or decrease larger than a factor 1/15:
      # In other words, check that k.new >= 1/15 * k.old & k.new <= 15 * k.old
      # deltaFactor <- delta[vars, 1] / k[vars, iter]
      deltaFactor <- kUpdate[vars, 1] / k[vars, iter]
      if (any(deltaFactor > 15 | deltaFactor < 1 / 15)) {
        con.vars.up <- vars[deltaFactor > 15] # Indices of variances that go up with more than a factor 15
        con.vars.down <- vars[deltaFactor < 1 / 15] # Indices of variances that go down with more than a factor 1 / 15
        con.vars <- sort(c(con.vars.up, con.vars.down)) # All indices of constrained variances
        # con.vars <- vars[deltaFactor > 15 | deltaFactor < 1 / 15]
        delta[con.vars.up, 1] <- (14 - 1e-5) * k[con.vars.up, iter] # (k + d) / k = 15
        delta[con.vars.down, 1] <- -(14 - 1e-5) / 15 * k[con.vars.down, iter] # (k + d) / k = 1 / 15
        # delta[con.vars, 1] <- sign(delta[con.vars, 1]) * 15 * k[con.vars, iter]
        constrained <- TRUE
        checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
        constraintsUsed <- TRUE
      }
      
      # Checking correlations for absolute changes larger than 0.3:
      deltaChange <- delta[cors, 1]
      if (any(abs(deltaChange) > 0.3)) {
        con.cors <- cors[abs(deltaChange) > 0.3]
        delta[con.cors, 1] <- sign(delta[con.cors, 1]) * 0.3
        constrained <- TRUE
        checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
        constraintsUsed <- TRUE
      }
      
      # Checking autocorrelations for absolute changes larger than 0.3:
      deltaChange <- delta[autoCors, 1]
      if (any(abs(deltaChange) > 0.3)) {
        con.autoCors <- autoCors[abs(deltaChange) > 0.3]
        delta[con.autoCors, 1] <- sign(delta[con.autoCors, 1]) * 0.3
        constrained <- TRUE
        checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
        constraintsUsed <- TRUE
      }
      
      # Checking bandwidths for absolute changes larger than 0.3:
      deltaChange <- delta[bws, 1]
      if (any(abs(deltaChange) > 0.3)) {
        con.bws <- bws[abs(deltaChange) > 0.3]
        delta[con.bws, 1] <- sign(delta[con.bws, 1]) * 0.3
        constrained <- TRUE
        checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
        constraintsUsed <- TRUE
      }
      
      # Computing the proposed update so we can check for boundaries:
      kUpdate <- k[, iter, drop = FALSE] + delta
      
      # Checking variances for boundaries:
      if (any(kUpdate[vars, 1] < 1e-6)) {
        bound.vars <- vars[kUpdate[vars, 1] < 1e-6]
        boundary[bound.vars] <- boundary[bound.vars] + 1
        delta[bound.vars, 1] <- 1e-5 - k[bound.vars, iter]
        hitBoundary <- TRUE
        checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
        constraintsUsed <- TRUE
      }
      
      # Checking correlations for boundaries:
      if (any(abs(kUpdate[cors, 1]) > 0.99)) {
        bound.cors <- cors[abs(kUpdate[cors, 1]) > 0.99]
        boundary[bound.cors] <- boundary[bound.cors] + 1
        signs <- sign(kUpdate[bound.cors, 1])
        delta[bound.cors, 1] <- signs * 0.9333333 - k[bound.cors, iter]
        hitBoundary <- TRUE
        checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
        constraintsUsed <- TRUE
      }
      
      # Checking autocorrelations for boundaries:
      if (any(kUpdate[autoCors, 1] > 0.99)) {
        bound.autoCors <- autoCors[kUpdate[autoCors, 1] > 0.99]
        boundary[bound.autoCors] <- boundary[bound.autoCors] + 1
        delta[bound.autoCors, 1] <- 0.98 - k[bound.autoCors, iter]
        hitBoundary <- TRUE
        checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
        constraintsUsed <- TRUE
      }
      if (any(kUpdate[autoCors, 1] < 0)) {
        bound.autoCors <- autoCors[kUpdate[autoCors, 1] < 0]
        boundary[bound.autoCors] <- boundary[bound.autoCors] + 1
        delta[bound.autoCors, 1] <- 0 - k[bound.autoCors, iter]
        hitBoundary <- TRUE
        checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
        constraintsUsed <- TRUE
      }
      
      # Checking bandwidths for boundaries:
      if (any(kUpdate[bws, 1] < 1e-3)) {
        bound.bws <- bws[kUpdate[bws, 1] < 1e-3]
        boundary[bound.bws] <- boundary[bound.bws] + 1
        delta[bound.bws, 1] <- 1e-2 - k[bound.bws, iter]
        hitBoundary <- TRUE
        checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
        constraintsUsed <- TRUE
      }
      
      # Recomputing the proposed update so we can check the loadings and whether all corg matrices are PD:
      kUpdate <- k[, iter, drop = FALSE] + delta
      
      # Checking loadings ====
      # Checking the loadings is a little bit more complicated. We check three things:
      # (1) The mean magnitude of the changes cannot be larger than 10. If it is, deltas for loadings are shrunken so the mean magnitude is 10.
      # (2) The magnitude of the change for any sign-flipped loading cannot be larger than 2. If any is, the delta for that loading is shrunken so the magnitude is 2.
      # (3) The communalities for all features (environments) must be smaller than 1 - 1e-3. If any communality is larger than that threshold, all loadings for that feature i
      #     are scaled back so that the communality becomes 1 - 1e-3, but the direction of the vector L[i,] in m-dimensional space does not change.
      # See the checkLoadings() function for what checks are done in more detail.
      for (re in 1:nU) {
        if (any(u.isStack[[re]])) { # If there is a stack in the random effect.
          # First we check any potential fa/exfa structures in the stack:
          i.stack <- which(u.isStack[[re]])
          for (i in which(struc[[re]][[i.stack]] %in% c("fa", "exfa"))) { # Check fa/exfa structures in the stack.
            indices.stack.start <- 1 + c(0, cumsum(u.nparams.list[[re]][[1 + i.stack]]))[i]
            indices.stack.end <- indices.stack.start - 1 + u.nparams.list[[re]][[1 + i.stack]][i]
            indices.stack <- indices.stack.start:indices.stack.end # Indices of the parameters that belong to this layer of the stack.
            loadingIndices <- indices[[re]][[1 + i.stack]][indices.stack] # Getting the indices for the loadings for this particular structure.
            check <- checkLoadings(currentParams = k[, iter, drop = FALSE], delta = delta, loadingIndices = loadingIndices, isFixed = isFixed, p = u.nlevels[[re]][i.stack])
            delta <- check$delta
            constrained <- (constrained || check$constrained)
            checkConstraints <- (checkConstraints || check$checkConstraints)
            constraintsUsed <- (constraintsUsed || check$constraintsUsed)
            hitBoundary <- (hitBoundary || check$hitBoundary)
            con.loadings <- sort(unique(c(con.loadings, check$con.loadings)))
            bound.loadings <- sort(unique(c(bound.loadings, check$bound.loadings)))
          }
        }
        # Then we check potential fa or exfa structures in the remaining factors.
        # Note that that block above never runs for non-stack random effects.
        for (i in which(struc[[re]] %in% c("fa", "exfa"))) {
          loadingIndices <- indices[[re]][[1 + i]] # Getting the indices for the loadings for this particular structure.
          check <- checkLoadings(currentParams = k[, iter, drop = FALSE], delta = delta, loadingIndices = loadingIndices, isFixed = isFixed, p = u.nlevels[[re]][i])
          delta <- check$delta
          constrained <- (constrained || check$constrained)
          checkConstraints <- (checkConstraints || check$checkConstraints)
          constraintsUsed <- (constraintsUsed || check$constraintsUsed)
          hitBoundary <- (hitBoundary || check$hitBoundary)
          con.loadings <- sort(unique(c(con.loadings, check$con.loadings)))
          bound.loadings <- sort(unique(c(bound.loadings, check$bound.loadings)))
        } # End of loop over fa/exfa structures for non-stack random effects.
      } # End of loop over random effects.
      boundary[bound.loadings] <- boundary[bound.loadings] + 1
      
      # Recomputing update:
      kUpdate <- k[, iter, drop = FALSE] + delta
      
      # Checking whether all unstructured correlation matrices are PD:
      for (re in 1:nUR) {
        if (any(u.isStack[[re]])) {
          i.stack <- which(u.isStack[[re]])
          sos <- length(u.nparams.list[[re]][[i.stack + 1]]) # + 1 due to the variances
          
          # First we check the stack:
          if (any(struc[[re]][[i.stack]] == "corg")) {
            s <- which(struc[[re]][[i.stack]] == "corg") # The stack can obviously only have a single unstructured correlation matrix, we just need to figure out where in the stack it is.
            cormats[[re]][[i.stack]][[s]] <- Matrix(diag(nrow(cormats[[re]][[i.stack]][[s]])), sparse = FALSE, doDiag = FALSE)
            if (s == 1) {
              indices.start <- 1
            } else {
              indices.start <- sum(u.nparams.list[[re]][[1 + i.stack]][1:(s - 1)]) + 1
            }
            indices.stack <- (indices.start):(indices.start - 1 + u.nparams.list[[re]][[1 + i.stack]][s]) # We need to get the correct indices.
            cormats[[re]][[i.stack]][[s]][lower.tri(cormats[[re]][[i.stack]][[s]])] <- kUpdate[indices[[re]][[1 + i.stack]][indices.stack], 1] # Now construct the proposed corg matrix.
            cormats[[re]][[i.stack]][[s]] <- forceSymmetric(cormats[[re]][[i.stack]][[s]], uplo = "L") # Symmetrising.
            ev <- eigen(cormats[[re]][[i.stack]][[s]], only.values = TRUE)$values
            if (min(ev) <= 0) {
              nonPD <- TRUE
              checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
              constraintsUsed <- TRUE
              npd <- indices[[re]][[1 + i.stack]][indices.stack]
              npd.cors <- c(npd.cors, npd)
              cormats[[re]][[i.stack]][[s]] <- nearPD(cormats[[re]][[i.stack]][[s]], corr = TRUE, posd.tol = 1e-3)$mat
              delta[npd, 1] <- cormats[[re]][[i.stack]][[s]][lower.tri(cormats[[re]][[i.stack]][[s]])] - k[npd, iter]
            }
          }
          
          # Then we check the remaining structures:
          if (any(struc[[re]] == "corg")) { # Note that the entry for the stack will always be FALSE: list(c("gk", "corg"), "cors", "lk") == "corg" -> FALSE FALSE FALSE
            for (i in which(struc[[re]] == "corg")) { # Same thing here, i will never be the index of the stack.
              cormats[[re]][[i]] <- Matrix(diag(nrow(cormats[[re]][[i]])), sparse = FALSE, doDiag = FALSE)
              cormats[[re]][[i]][lower.tri(cormats[[re]][[i]])] <- kUpdate[indices[[re]][[1 + i]], 1]
              cormats[[re]][[i]] <- forceSymmetric(cormats[[re]][[i]], uplo = "L")
              ev <- eigen(cormats[[re]][[i]], only.values = TRUE)$values
              if (min(ev) <= 0) {
                nonPD <- TRUE
                checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
                constraintsUsed <- TRUE
                npd <- indices[[re]][[1 + i]]
                npd.cors <- c(npd.cors, npd)
                cormats[[re]][[i]] <- nearPD(cormats[[re]][[i]], corr = TRUE, posd.tol = 1e-3)$mat
                delta[npd, 1] <- cormats[[re]][[i]][lower.tri(cormats[[re]][[i]])] - k[npd, iter]
              }
            }
          }
        } else {
          if (any(struc[[re]] == "corg")) {
            for (i in which(struc[[re]] == "corg")) {
              cormats[[re]][[i]] <- Matrix(diag(nrow(cormats[[re]][[i]])), sparse = FALSE, doDiag = FALSE)
              cormats[[re]][[i]][lower.tri(cormats[[re]][[i]])] <- kUpdate[indices[[re]][[1 + i]], 1]
              cormats[[re]][[i]] <- forceSymmetric(cormats[[re]][[i]], uplo = "L")
              ev <- eigen(cormats[[re]][[i]], only.values = TRUE)$values
              if (min(ev) <= 0) {
                nonPD <- TRUE
                checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
                constraintsUsed <- TRUE
                npd <- indices[[re]][[1 + i]]
                npd.cors <- c(npd.cors, npd)
                cormats[[re]][[i]] <- nearPD(cormats[[re]][[i]], corr = TRUE, posd.tol = 1e-3)$mat
                delta[npd, 1] <- cormats[[re]][[i]][lower.tri(cormats[[re]][[i]])] - k[npd, iter]
              }
            }
          }
        }
      }
      
      # Now that we have (1) constrained for excessive changes, (2) constrained for boundary
      # issues, and (3) constrained non-PD correlation matrices, we can recompute delta for any
      # unconstrained parameters:
      if (constrained || hitBoundary || nonPD) {
        # con consists of indices of params that were constrained in the previous cycles, with the new constrained params added.
        con <- sort(unique(c(con, con.vars, bound.vars, con.cors, bound.cors, npd.cors, con.autoCors, bound.autoCors, con.bws, bound.bws, con.loadings, bound.loadings)))
        noncon <- setdiff(1:nrow(delta), con)
        nConstrained <- length(con)
        if (length(setdiff(noncon, which(isFixed))) > 0) {
          AIuu <- AI[noncon, noncon, drop = FALSE]
          uuisFixed <- isFixed[noncon]
          AIff <- AI[con, con, drop = FALSE]
          AIuf <- AI[noncon, con, drop = FALSE]
          scoresuu <- scores[noncon, 1, drop = FALSE]
          scoresff <- scores[con, 1, drop = FALSE]
          deltaff <- delta[con, 1, drop = FALSE]
          if (any(uuisFixed)) { # If there are any fixed parameters in the unconstrained part, remove them from the AIuu matrix before inverting.
            AIuui <- Matrix(0, nrow(AIuu), ncol(AIuu), sparse = FALSE, doDiag = FALSE)
            uunotFixed <- which(!uuisFixed)
            AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
            deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
          } else {
            deltauu <- solve(AIuu) %*% (scoresuu - AIuf %*% deltaff)
          }
          delta[noncon, 1] <- deltauu
        }
      }
      q <- q + 1
    }
    
    # Fixing any parameters that have hit the boundary five times:
    if (any(boundary == timesToBoundary)) {
      isFixed[which(boundary == timesToBoundary)] <- TRUE
    }
    if (timing) {timer("stop")}
    
    # Printing:
    ##tic("Printing info")
    stop <- proc.time()[3]
    duration <- as.numeric(stop - start)
    duration <- sprintf("%02d:%02d:%06.3f", duration %/% 3600, (duration %% 3600) %/% 60, duration %% 60)
    if (verbose) {
      if (constraintsUsed) {
        msg <- paste0("(", paste(c(
          if (constraintsUsed) sprintf("%d constrained", nConstrained)
        ), collapse = " and "), ")")
        cat(sprintf("%-11d %-15.3f %-20s %-5s %-s\n", iter, LL[iter], duration, "", msg), sep = "")
      } else {
        cat(sprintf("%-11d %-15.3f %-20s\n", iter, LL[iter], duration))
      }
    }
    ##toc()
    
    # Checking for convergence:
    if (iter > 1) {
      deltaLL <- LL[iter] - LL[iter - 1]
      if (deltaLL > 0 && deltaLL <= tolConvLL && !(constrained || hitBoundary)) {
        LL <- LL[1:(iter + 1)]
        k <- k[, 1:(iter + 1)]
        LL[iter + 1] <- LL[iter]
        k[, iter + 1] <- k[, iter, drop = FALSE] + delta
        converged <- TRUE
      }
    }
    
    # Final update to parameters if not converged and not reached maxit:
    if (iter <= maxit && !(converged)) {
      k[, iter + 1] <- k[, iter, drop = FALSE] + delta
      iter <- iter + 1
    }
  }
  
  # Now that when we have reached maxit or convergence, we need to compute Cchol one
  # last time if we are using the GPU. This is because the Cholesky factor of C
  # only ever existed on the GPU, and we don't want to copy it back in every
  # single iteration.
  if (useGPU && Cmode == "dense") {
    gpuResult <- cuSolverPhase1(ctx_ptr = gpu_ctx, C_in = Cdd, RHS_in = RHS, returnCchol = TRUE)
    cuSolverCleanup(gpu_ctx)
    Ccholdd <- gpuResult$Cchol
  }
  
  # Putting the BLUEs and BLUPs back into specified order:
  if (Cmode == "mixed") {
    bu <- rbind(buDense, buSparse)
  } else if (Cmode == "sparse") {
    bu <- buSparse
  } else {
    bu <- buDense
  }
  
  mapping <- match(1:nrow(bu), c(dense.idx, sparse.idx))
  bu <- bu[mapping, 1, drop = FALSE]
  b <- bu[1:nX, 1, drop = FALSE]
  u <- bu[(nX + 1):nrow(bu), 1, drop = FALSE]
  
  if (Cmode == "mixed") {
    infoPEV <- list(Cmode = "mixed",
                    Ccholss = Ccholss,
                    Q = Q,
                    Cidd = Cidd,
                    dense.idx = dense.idx,
                    sparse.idx = sparse.idx)
    
  } else if (Cmode == "sparse") {
    infoPEV <- list(Cmode = "sparse",
                    Ccholss = Ccholss,
                    sparse.idx = sparse.idx)
    
  } else {
    infoPEV <- list(Cmode = "dense",
                    Ccholdd = Ccholdd,
                    dense.idx = dense.idx)
  }
  
  timings <- NULL
  if (timing) {
    timings <- timer("report")
  }
  
  if (debugging) {
    return(list(LL = LL,
                k = k,
                AI = AI,
                scores = scores,
                b = b,
                u = u,
                infoPEV = infoPEV,
                isFixed = isFixed,
                converged = converged,
                derivs = derivs,
                sigmai = sigmai,
                vars = vars.debugging,
                cormats = cormats,
                logDetC = logDetC,
                logDetR = logDetR,
                logDetG = logDetG,
                yPy = yPy,
                traces = traces.debugging,
                logDetvarGtotal = logDetvarGtotal,
                logDetnvarGtotal = logDetnvarGtotal,
                X = X,
                Z = Z,
                y = y,
                Wu = Wu,
                e = e,
                bases = bases,
                timings = timings))
  } else {
    return(list(LL = LL,
                k = k,
                AI = AI,
                scores = scores,
                b = b,
                u = u,
                infoPEV = infoPEV,
                isFixed = isFixed,
                converged = converged,
                timings = timings))
  }
}

