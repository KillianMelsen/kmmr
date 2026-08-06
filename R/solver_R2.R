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
    delta[magOver2, 1] <- -2 * currentParams[magOver2, 1] - currentParams[magOver2, 1]
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

# A simpler wrapper for the Takahashi equations:
# takahashi <- function(L) {
#   # sparseinv_wrapper <- utils::getFromNamespace(".sparseinv_wrapper", "sparseinv")
#   n <- nrow(L)
#   d <- diag(L)
#   L <- tril(L %*% Diagonal(x = 1 / d), -1)
#   d <- d^2
#   D <- Diagonal(x = d)
#   dp <- diff(L@p)
#   jj <- rep(seq_along(dp), dp)
#   Zpattern <- sparseMatrix(c(L@i + 1, jj, 1:n), c(jj, L@i + 1, 1:n))
#   Z <- forceSymmetric(sparseinv:::.sparseinv_wrapper(L, d, L, Zpattern))
#   # Z <- forceSymmetric(sparseinv_wrapper(L, d, L, Zpattern))
#   return(Z)
# }

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
  }
  return(ret)
}

getDeriv_cor <- function(cormatIndex, p) {
  deriv <- Matrix(0, p, p, sparse = TRUE, doDiag = FALSE)
  deriv[lower.tri(deriv)][cormatIndex] <- 1
  deriv <- forceSymmetric(deriv, uplo = "L")
  return(deriv)
}

getDeriv_gk <- function(params, base) {
  deriv <- -base * exp(-params[1] * base)
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
  deriv <- forceSymmetric(as(deriv, "CsparseMatrix"))
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
    deriv <- forceSymmetric(as(deriv, "CsparseMatrix"))
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
  }
  return(ret)
}

#' solveLMM_R
#' 
#' R implementation of the average information algorithm for fitting linear mixed
#' models.
#'
#' @param y Numeric vector containing the response data.
#' @param X Numeric fixed effect design matrix.
#' @param Z List containing numeric random effect design matrices. One for each random effect.
#' @param S List containing numeric residual covariance matrices.
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
#' @param bases Base matrices
#' @param struc Structure for each factor.
#' @param i.var Which factors have variable covariance structures?
#' @param i.nvar Which factors have non-variable covariance structures?
#' @param isFixed Which parameters are fixed and should not be estimated?
#' @param u.isStack Are factors within random effects made up of a stack of structures?
#' @param u.nparams.list Number of parameters per factor/random effect in list format.
#' @param sparse Should the sparse version of the algorithm be used? Useful if C is reasonable sparse.
#' @param timesToBoundary How many times is a parameter or set of loadings allowed to cross the boundary before it is fixed?
#' @param super Should a supernodal (\code{super = TRUE}) or simplicial Cholesky factor be computed?
#'
#' @returns A list containing the log-likelihood, parameter estimates, average information matrix,
#' scores, BLUEs for fixed effects, and BLUPs for random effects.
#' @keywords internal
solveLMM_R <- function(y, X, Z, S, init, weights, useWeights, maxit, stepSize, avInfReg, verbose, tolConvLL,
                       indices, u.applyvars, ur.types, u.nlevels, u.nparams, u.forFactor, u.indexCor,
                       bases, struc, i.var, i.nvar, isFixed, u.isStack, u.nparams.list, sparse,
                       timesToBoundary, super) {
  ## Prelims ====
  k <- Matrix(0, length(init), maxit + 1, sparse = FALSE, doDiag = FALSE) # Matrix of parameter estimates.
  k[, 1] <- init # Filling with initial values.
  boundary <- numeric(length(init)) # Vector to hold how many times each parameter has crossed the boundary.
  n <- nrow(y)
  
  nU <- length(Z) # Number of random effects
  nR <- length(S) # Number of residual effects
  nUR <- nU + nR # Number of random effects + residual
  nX <- ncol(X) # Number of fixed effects
  intercept <- sum(X[, 1]) == nrow(X) # Is there an intercept?
  
  # Assume for now we have just one random and residual effect:
  S <- S[[1]]
  
  ## Computing determinants of matrices associated with non-variable terms before iterating ====
  n.var <- numeric(nU)
  n.nvar <- numeric(nU)
  dim.var <- numeric(nU)
  dim.nvar <- numeric(nU)
  u.nlevels.var <- vector("list", nU) 
  u.nlevels.nvar <- vector("list", nU) 
  for (re in 1:nU) { # Start loop over random effects
    n.var[re] <- length(i.var[[re]]) # How many variable factors are there?
    n.nvar[re] <- length(i.nvar[[re]]) # How many non-variable factors are there?
    u.nlevels.var[[re]] <- u.nlevels[[re]][i.var[[re]]] # Vector of levels for each variable factor for each random effect.
    u.nlevels.nvar[[re]] <- u.nlevels[[re]][i.nvar[[re]]] # Vector of levels for each non-variable factor for each random effect.
    dim.var[re] <- prod(u.nlevels.var[[re]]) # Total dimension of the variable part of each random effect.
    dim.nvar[re] <- prod(u.nlevels.nvar[[re]]) # Total dimension of the non-variable part of each random effect.
  } # End loop over random effects
  
  nvarGi <- vector("list", nU)
  logDetnvarGtotal <- numeric(nU)
  for (re in 1:nU) { # Start loop over random effects
    if (n.nvar[re] > 0) { # All of this can only be done if there are non-variable factors (such as a kinship/linear kernel for genotypes).
      logDetnvarG <- numeric(n.nvar[re]) # A log determinant for each non-variable base matrix (linear kernel).
      for (i in i.nvar[[re]]) {
        if (!is.null(attr(bases[[re]][[i]], "logDet"))) {
          logDetnvarG[i - n.var[re]] <- attr(bases[[re]][[i]], "logDet")
        } else {
          logDetnvarG[i - n.var[re]] <- -1 * sum(log(eigen(bases[[re]][[i]])$values)) # Times -1 because the non-variable base matrices (linear kernels) are in inverse form. # Maybe on the GPU for large matrices? 
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
  
  ## Start iterating ====
  LL <- numeric(maxit + 1)
  converged <- FALSE
  iter <- 1
  while (!(converged) & iter <= maxit) {
    start <- Sys.time()
    ### Compute Ri ====
    #tic("Computing Ri")
    if (useWeights & iter == 1) { # Apply weights in the first iteration if they need to be used.
      S <- t(weights) %*% S %*% weights
    }
    Ri <- 1 / k[nrow(k), iter] * S # The residual variance is always last in kappa for now
    #toc(log = T)
    
    ### Forming W = [X Z] ====
    #tic("Forming W and C = t(W) %*% Ri %*% W")
    if (iter == 1) {
      W <- X
      for (re in 1:nU) {
        W <- cbind(W, Z[[re]])
      }
    }
    
    ### Computing C = t(W) %*% Ri %*% W
    if (sparse) {
      C <- as(forceSymmetric(crossprod(W, Ri %*% W)), "TsparseMatrix") # Converting to triplet form for adding sparse Gi's later
      # C_i <- C@i
      # C_j <- C@j
      # C_x <- C@x
      C_dim <- dim(C)
      # rm(C)
    } else {
      C <- as(forceSymmetric(crossprod(W, Ri %*% W)), "unpackedMatrix") # Store as a dense matrix since C will be dense anyway after adding Gi.
    }
    # if (useWeights) {
    #   # We need to compute the entire thing every time? Hopefully not for diagonal R...
    #   C <- t(W) %*% Ri %*% W
    # } else {
    #   if (iter == 1) {
    #     C0 <- t(W) %*% W # Maybe use something like syrk()?
    #   }
    #   C <- C0 * (1 / k[nrow(k), iter])
    # }
    #toc(log = T)
    
    ### Adding Gi to C ====
    #tic("Computing Gi")
    cormats <- vector("list", nU)
    sigmai <- vector("list", nU)
    logDetvarGtotal <- numeric(nU)
    vars <- vector("list", nU)
    offset <- nX + 1
    for (re in 1:nU) { # Start loop over random effects
      # Computing Gi every iteration:
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
        if (n.nvar[re] > 0) {
          Gi <- as(kronecker(sigmai[[re]], nvarGi[[re]]), "TsparseMatrix") # This is a bit slow as well.
        } else {
          Gi <- as(sigmai[[re]], "TsparseMatrix")
        }

        # Now add to C:
        dimGi <- nrow(Gi)
        if (sparse) {
          # If we want the full C to be sparse but still have a dense random effect, we need to store the Gi for that
          # dense random effect in sparse format.
          if (!is(Gi, "sparseMatrix")) {
            Gi <- as(Gi, "TsparseMatrix")
          }
          # Collecting the values of the previous C and current Gi.
          # Padding Gi with zeros to make it conformable with the sparse base C.
          C_i_combined <- c(C@i, Gi@i + offset - 1)
          C_j_combined <- c(C@j, Gi@j + offset - 1)
          C_x_combined <- c(C@x, Gi@x)
          rm(Gi, C) # We have all we need in the three vectors above, so get rid of Gi and the old C to free up memory.
          
          # Constructing the new C:
          if (re == nU) {
            repr <- "C" # Finally convert C to dsCMatrix to prepare for CF if there are no more random effects.
          } else {
            repr <- "T" # Keep C in triplet form if we have more random effects.
          }
          C <- sparseMatrix(i = C_i_combined + 1,
                            j = C_j_combined + 1,
                            x = C_x_combined,
                            dims = C_dim,
                            symmetric = TRUE,
                            repr = repr)
          
          # Alternative approach, collecting all i's, j's, x's and building C after looping
          # over all random effects.
          # C_i <- c(C_i, Gi@i + offset - 1)
          # C_j <- c(C_j, Gi@j + offset - 1)
          # C_x <- c(C_x, Gi@x)
          
          # Gi <- sparseMatrix(i = Gi@i + 1 + offset - 1,
          #                    j = rep(seq_len(ncol(Gi)), diff(Gi@p)) + offset - 1,
          #                    x = Gi@x,
          #                    dims = dim(C),
          #                    symmetric = TRUE)
          # C <- C + Gi
        } else {
          C[offset:(offset - 1 + dimGi), offset:(offset - 1 + dimGi)] <- C[offset:(offset - 1 + dimGi), offset:(offset - 1 + dimGi)] + Gi
        }
        offset <- offset + dimGi # Moving offset to the next random effect.
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
        if (n.nvar[re] > 0) {
          Gi <- as(kronecker(sigmai[[re]], nvarGi[[re]]), "TsparseMatrix") # This is a bit slow as well.
        } else {
          Gi <- as(sigmai[[re]], "TsparseMatrix")
        }
        
        # Now add to C:
        dimGi <- nrow(Gi)
        if (sparse) {
          # If we want C to be sparse but still have a dense random effect, we need to store the Gi for that
          # dense random effect in sparse format.
          if (!is(Gi, "sparseMatrix")) {
            Gi <- as(Gi, "TsparseMatrix")
          }
          # Collecting the values of the previous C and current Gi.
          # Padding Gi with zeros to make it conformable with the sparse base C.
          C_i_combined <- c(C@i, Gi@i + offset - 1)
          C_j_combined <- c(C@j, Gi@j + offset - 1)
          C_x_combined <- c(C@x, Gi@x)
          rm(Gi, C) # We have all we need in the three vectors above, so get rid of Gi and the old C to free up memory.
          
          # Constructing the new C:
          if (re == nU) {
            repr <- "C" # Finally convert C to dsCMatrix to prepare for CF if there are no more random effects.
          } else {
            repr <- "T" # Keep C in triplet form if we have more random effects.
          }
          C <- sparseMatrix(i = C_i_combined + 1,
                            j = C_j_combined + 1,
                            x = C_x_combined,
                            dims = C_dim,
                            symmetric = TRUE,
                            repr = repr)
          
          # Alternative approach, collecting all i's, j's, x's and building C after looping
          # over all random effects.
          # C_i <- c(C_i, Gi@i + offset - 1)
          # C_j <- c(C_j, Gi@j + offset - 1)
          # C_x <- c(C_x, Gi@x)
          
          # Padding Gi with zeros to make it conformable with the sparse base C.
          # Gi <- sparseMatrix(i = Gi@i + 1 + offset - 1,
          #                    j = rep(seq_len(ncol(Gi)), diff(Gi@p)) + offset - 1,
          #                    x = Gi@x,
          #                    dims = dim(C),
          #                    symmetric = TRUE)
          # C <- C + Gi
        } else {
          C[offset:(offset - 1 + dimGi), offset:(offset - 1 + dimGi)] <- C[offset:(offset - 1 + dimGi), offset:(offset - 1 + dimGi)] + Gi # This one is also slow
        }
        offset <- offset + dimGi # Moving offset to next random effect.
      }
    } # End loop over random effects.
    # rm(Gi) # No more need for Gi, so get rid of it.
    
    # If C is sparse, we need to actually construct it once we have collected all
    # i's, j's, and x's for all random effects:
    # if (sparse) {
    #   C <- sparseMatrix(i = C_i + 1,
    #                     j = C_j + 1,
    #                     x = C_x,
    #                     dims = C_dim,
    #                     symmetric = TRUE,
    #                     repr = "C")
    # }
    #toc(log = T)
    
    ### Cholesky decomposition of C ====
    # Overwrite C with its simplicial/supernodal Cholesky factor.
    # We don't need C anymore.
    #tic(msg = "Cholesky decomp of C")
    if (sparse) {
      if (super) {
        C <- Cholesky(C, LDL = FALSE, perm = TRUE, super = TRUE)
      } else {
        C <- Cholesky(C, LDL = TRUE, perm = TRUE, super = FALSE)
      }
      # if (iter == 1) {
      #   Cchol <- Cholesky(C, LDL = FALSE, perm = TRUE, super = FALSE)
      # } else {
      #   # To Do: only update the cholesky factor instead
      # }
    } else {
      C <- Cholesky(C, perm = FALSE) # Good for the GPU?
    }
    #toc(log = T)
    
    ### BLUEs/BLUPs ====
    #tic("Solving for BLUEs and BLUPs")
    if (sparse) {
      # Apply permutation to the right-hand side of the MME: P %*% t(W) %*% Ri %*% y:
      permRHS <- solve(C, crossprod(W, Ri %*% y), system = "P")
      # Sparse triangular forward solve: L %*% z = permRHS where L = C (remember that C is now actually the CF of C):
      z <- solve(C, permRHS, system = "L")
      # Sparse triangular backward solve: t(L) %*% bu = z where l = C to obtain BLUEs and BLUPs:
      bu <- solve(C, z, system = "Lt")
      # Final reverse permutation t(P) %*% bu to get the BLUEs and BLUPs into model order:
      bu <- solve(C, bu, system = "Pt")
    } else {
      RHS <- crossprod(W, Ri %*% y)
      bu <- solve(C, RHS)
      z <- as(C, "dtrMatrix") %*% bu # t(U) %*% U %*% bu = RHS and t(U) %*% z = RHS, therefore, t(U) %*% U %*% bu = t(U) %*% z, thus U %*% bu = z
    }
    
    # Splitting BLUEs and BLUPs:
    b <- bu[1:nX, 1, drop = FALSE]
    u <- bu[(nX + 1):nrow(bu), 1, drop = FALSE]
    #toc(log = T)
    
    ### Computing the log-likelihood ====
    #tic("Computing log-likelihood")
    # yPy:
    yPy <- as.numeric(crossprod(y, Ri %*% y) - crossprod(z))
    
    # logdetC:
    if (sparse) {
      logDetC <- determinant(C, sqrt = FALSE, logarithm = TRUE)$modulus
    } else {
      # Looks strange because the 2 * is missing, but is correct as diag(Cholesky(C, perm = FAlSE))
      # returns the squared diagonal elements of the Cholesky factor, see ?Matrix::`Cholesky-class`
      logDetC <- sum(log(diag(C)))
    }
    
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
    
    # logDetR:
    if (useWeights) {
      # For now just assume a diagonal structure:
      logDetR <- sum(log(diag(Ri))) * -1
    } else {
      logDetR <- n * log(k[nrow(k), iter])
    }
    
    # Full LL:
    LL[iter] <- -0.5 * (logDetC + logDetR + logDetG + yPy)
    #toc(log = T)
    
    ### Working variates ====
    #### Genetic working variates ====
    #tic("Genetic working variates")
    Wu <- Matrix(0, n, nrow(k), sparse = FALSE, doDiag = FALSE) # Matrix for the working variates.
    umat_sigmai <- vector("list", nU) # We will store these matrices for computing the scores later.
    derivs <- vector("list", nU) # Same thing for all derivatives.
    u.start <- 1
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
      if (n.nvar[re] > 0) {
        umat <- Matrix(u[u.start:(u.start - 1 + dim.total), 1], nrow = dim.nvar[re], ncol = dim.var.corrected)
      } else {
        umat <- Matrix(u[u.start:(u.start - 1 + dim.total), 1], nrow = 1, ncol = dim.var[re])
      }
      u.start <- u.start + dim.total
      umat_sigmai[[re]] <- umat %*% sigmai[[re]]
      derivs[[re]] <- vector("list", u.nparams[[re]])
      
      if (any(u.isStack[[re]])) { # If there is a stack in this random effect
        # Determining which factor is a stack and the size of the stack:
        i.stack <- which(u.isStack[[re]])
        sos <- length(u.nparams.list[[re]][[i.stack + 1]]) # + 1 due to the variances
        nvars <- length(indices[[re]][[1]]) / sos # Number of variances per stack layer
        for (i in 1:u.nparams[[re]]) {
          if (isFixed[param]) {
            param <- param + 1
          } else {
            if (ur.types[param] == "V") {
              # We need to figure out for which stack layer the current variance is:
              s <- ceiling(i / nvars)
              if (u.applyvars[[re]][1] == 0) { # If there is only a single variance.
                derivs[[re]][[i]] <- Reduce(kronecker, c(cormats[[re]][[i.stack]][s], cormats[[re]][-i.stack])) # The stack is always the first factor.
              } else { # If variances are applied to the factor in the stack.
                # In this case u.applyvars[[re]] is always equal to 1.
                A <- Matrix(0, prod(u.nlevels[[re]][u.applyvars[[re]]]), prod(u.nlevels[[re]][u.applyvars[[re]]]), sparse = FALSE, doDiag = FALSE)
                A[i - (s - 1) * nvars,] <- 0.5 # we can't use i directly
                if (n.var[re] > 1) { # If there are more variable factors besides the stack.
                  derivs[[re]][[i]] <- kronecker((A + t(A)) * vars[[re]][[s]] * cormats[[re]][[u.applyvars[[re]]]][[s]] / k[param, iter], Reduce(kronecker, cormats[[re]][-u.applyvars[[re]]]))
                } else { # If the stack is the only variable factor.
                  derivs[[re]][[i]] <- (A + t(A)) * vars[[re]][[s]] * cormats[[re]][[u.applyvars[[re]]]][[s]] / k[param, iter]
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
                  derivs[[re]][[i]] <- vars[[re]][[s]]@x * Reduce(kronecker, cormats1)
                } else { # If there is a variance for each level of the stacked factor.
                  derivs[[re]][[i]] <- Reduce(kronecker, c(list(vars[[re]][[s]] * cormats1[[i.stack]]), cormats1[-i.stack]))
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
                derivs[[re]][[i]] <- kronecker(covmat.stack, Reduce(kronecker, cormats1[-i.stack]))
              }
            }
            Wu[, param] <- Z[[re]] %*% Matrix((umat_sigmai[[re]] %*% derivs[[re]][[i]])@x)
            param <- param + 1
          }
        } # End of loop over the parameters for this random effect.
      } else { # If there is no stack
        for (i in 1:u.nparams[[re]]) {
          if (isFixed[param]) {
            param <- param + 1
          } else {
            if (ur.types[param] == "V") { # If the parameter is a variance.
              if (u.applyvars[[re]][1] == 0) { # If there is only a single variance.
                if (n.var[re] > 0) {
                  derivs[[re]][[i]] <- Reduce(kronecker, cormats[[re]])
                } else {
                  derivs[[re]][[i]] <- Matrix(1, sparse = FALSE, doDiag = FALSE)
                }
              } else { # If variances are applied to 1 or more factors.
                # If there is no exfa structure, m = 0 and nothing changes.
                # If there is an exfa structure, u.applyvars[[re]] = 1 by definition and m is added to the number of levels of the factor modeled by the exfa structure.
                A <- Matrix(0, prod(u.nlevels[[re]][u.applyvars[[re]]]) + m, prod(u.nlevels[[re]][u.applyvars[[re]]]) + m, sparse = FALSE, doDiag = FALSE)
                A[i,] <- 0.5
                if (all(u.applyvars[[re]] == i.var[[re]])) { # If variances are applied to all variable factors.
                  derivs[[re]][[i]] <- (A + t(A)) * vars[[re]] * Reduce(kronecker, cormats[[re]]) / k[param, iter]
                } else { # If variances are applied to 1 or more factors, but not all factors.
                  derivs[[re]][[i]] <- kronecker((A + t(A)) * vars[[re]] * Reduce(kronecker, cormats[[re]][u.applyvars[[re]]]) / k[param, iter], Reduce(kronecker, cormats[[re]][-u.applyvars[[re]]]))
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
                derivs[[re]][[i]] <- vars[[re]]@x * Reduce(kronecker, cormats1)
              } else if (all(u.applyvars[[re]] == i.var[[re]])) {
                derivs[[re]][[i]] <- vars[[re]] * Reduce(kronecker, cormats1)
              } else {
                derivs[[re]][[i]] <- kronecker(vars[[re]] * Reduce(kronecker, cormats1[u.applyvars[[re]]]), Reduce(kronecker, cormats1[-u.applyvars[[re]]]))
              }
            }
            Wu[, param] <- Z[[re]] %*% Matrix((umat_sigmai[[re]] %*% derivs[[re]][[i]])@x)
            param <- param + 1
          }
        }
      }
    } # End of loop over random effects
    #toc(log = T)
    
    #### Residual working variates ====
    #tic("Residual working variates")
    e <- y - W %*% bu
    if (isFixed[nrow(k)]) {
      # Wu[, nrow(k)] <- 0
    } else {
      Wu[, nrow(k)] <- (1 / k[nrow(k), iter]) * e@x # I think this is always true for diagonal R (also for non-diagonal fixed R, if that would ever be the case), even with weights.
    }
    #toc(log = T)
    
    ### Computing the AI matrix ====
    #tic("Computing the AI matrix and its inverse")
    WiWj <- forceSymmetric(crossprod(Wu, Ri %*% Wu)) # Similar to t(y) %*% Ri %*% y
    if (sparse) {
      XWjxZWj <- solve(C, crossprod(W, Ri %*% Wu), system = "P") # Similar to permRHS
      zWu <- solve(C, XWjxZWj, system = "L") # Similar to z
      buWu <- solve(C, zWu, system = "Lt") # Similar to obtaining bu using a backsolve
    } else {
      XWjxZWj <- crossprod(W, Ri %*% Wu) # t(W) %*% Ri %*% Wu
      buWu <- solve(C, XWjxZWj)
      # zWu <- backsolve(C, XWjxZWj, transpose = TRUE) # Forward solve of lower triangular matrix t(C)
      # buWu <- Matrix(backsolve(C, zWu))
    }
    buWuXWjxZWj <- forceSymmetric(crossprod(buWu, XWjxZWj))
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
    #toc(log = T)
    
    ### Computing the inverse of the coefficient matrix ====
    # We can again overwrite C as we no longer need the CF after this.
    #tic("Inversion of C from cholesky factor")
    if (sparse) {
      C <- SIS(C) # The cholmod_extra code already returns the SIS in original model order.
    } else {
      C <- chol2inv(as(C, "dtrMatrix"))
    }
    #toc(log = TRUE)
    
    ### Computing scores ====
    # tic("Computing genetic scores")
    # See equations (18) and (15) from Jensen, Mantysaari, Madsen, & Thompson 1997:
    # "Residual maximum likelihood estimation of (co) variance components in multivariate
    # mixed linear models using average information"
    #### Genetic scores ====
    scores <- Matrix(0, nrow(k), 1, sparse = FALSE)
    corner <- nX + 1
    param <- 1
    for (re in 1:nU) { # Start of loop over random effects
      if (any(struc[[re]] == "exfa")) {
        m <- u.nparams.list[[re]][[1 + 1]] / u.nlevels[[re]][1]
        dim.var.corrected <- prod(c(u.nlevels.var[[re]][1] + m, u.nlevels.var[[re]][-1]))
      } else {
        dim.var.corrected <- dim.var[re]
      }
      if (n.nvar[re] > 0) {
        if (sparse) {
          C_start <- corner
          C_end <- corner + dim.var.corrected * dim.nvar[re] - 1
          C_sub <- C[C_start:C_end, C_start:C_end, drop = FALSE]
          traces <- forceSymmetric(getTraces(C_sub@p,
                                             C_sub@i,
                                             C_sub@x,
                                             as.matrix(nvarGi[[re]]),
                                             dim.var.corrected,
                                             dim.nvar[re]))
        } else {
          traces <- Matrix(0, nrow(sigmai[[re]]), ncol(sigmai[[re]]), sparse = FALSE, doDiag = FALSE)
          for (i in 1:nrow(sigmai[[re]])) {
            for (j in 1:ncol(sigmai[[re]])) {
              if (j >= i) { # Only do the upper triangle
                rstart <- (corner + (i - 1) * dim.nvar[re])
                rend <- rstart - 1 + dim.nvar[re]
                cstart <- (corner + (j - 1) * dim.nvar[re])
                cend <- cstart - 1 + dim.nvar[re]
                traces[i, j] <- sum(nvarGi[[re]] * C[rstart:rend, cstart:cend])
              }
            }
          }
          if (nrow(traces) > 1) {
            traces <- forceSymmetric(traces)
          }
        }
        corner <- corner + dim.var.corrected * dim.nvar[re]
      } else {
        traces <- C[corner:(corner - 1 + dim.var.corrected), corner:(corner - 1 + dim.var.corrected)]
        corner <- corner + dim.var.corrected
      }
      
      if (n.nvar[re] > 0) {
        basisScores <- -0.5 * (dim.nvar[re] * sigmai[[re]] - sigmai[[re]] %*% traces %*% sigmai[[re]] - crossprod(umat_sigmai[[re]], nvarGi[[re]] %*% umat_sigmai[[re]])) # Last term maybe on the GPU?
      } else {
        basisScores <- -0.5 * (sigmai[[re]] - sigmai[[re]] %*% traces %*% sigmai[[re]] - crossprod(umat_sigmai[[re]])) # Last term maybe on the GPU?
      }
      for (i in 1:length(derivs[[re]])) {
        if (!isFixed[param]) {
          scores[param, 1] <- sum(derivs[[re]][[i]] * basisScores)
        }
        param <- param + 1
      }
    } # End of loop over random effects.
    #toc(log = T)
    
    #### Residual scores ====
    #tic("Computing residual scores")
    # Can probably be more effecient for diagonal matrices
    evari <- (1 / k[nrow(k), iter])
    if (!isFixed[nrow(k)]) {
      WSW <- t(W) %*% S %*% W
      scores[nrow(scores), 1] <- as.numeric(-0.5 * (sum(S * evari) - sum(evari^2 * C[cbind(WSW@i + 1, rep(seq_len(ncol(WSW)), diff(WSW@p)))] * WSW@x) - (evari^2 * t(e) %*% S %*% e)))
    }
    #toc(log = T)
    
    ### Updating parameters ====
    #tic("Updating parameters and checking constraints")
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
    con <- c() # Vector to keep the total set of constrained parameters across cycles. We need this for recomputing delta for unconstrained params.
    vars <- which(ur.types == "V") # Which of the parameters in k are variances?
    cors <- which(ur.types == "R") # Which of the parameters in k are correlations?
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
      con.bws <- c() # Vector with indices of bandwidths that were constrained in the current cycle due to the deltas being too large.
      con.loadings <- c() # Vector with indices of loadings that were constrained in the current cycle due to the deltas being too large.
      bound.vars <- c() # Vector with indices of variances that were constrained in the current cycle due to the deltas pushing them across boundaries.
      bound.cors <- c() # Vector with indices of correlations that were constrained in the current cycle due to the deltas pushing them across boundaries.
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
            constrained <- check$constrained
            checkConstraints <- check$checkConstraints
            constraintsUsed <- check$constraintsUsed
            hitBoundary <- check$hitBoundary
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
          constrained <- check$constrained
          checkConstraints <- check$checkConstraints
          constraintsUsed <- check$constraintsUsed
          con.loadings <- sort(unique(c(con.loadings, check$con.loadings)))
          bound.loadings <- sort(unique(c(bound.loadings, check$bound.loadings)))
        } # End of loop over fa/exfa structures for non-stack random effects.
      } # End of loop over random effects.
      boundary[bound.loadings] <- boundary[bound.loadings] + 1
      
      # Recomputing update:
      kUpdate <- k[, iter, drop = FALSE] + delta
      
      # Checking whether all unstructured correlation matrices are PD:
      for (re in 1:nU) {
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
      if (constrained | hitBoundary | nonPD) {
        # con consists of indices of params that were constrained in the previous cycles, with the new constrained params added.
        con <- sort(unique(c(con, con.vars, bound.vars, con.cors, bound.cors, npd.cors, con.bws, bound.bws, con.loadings, bound.loadings)))
        noncon <- setdiff(1:nrow(delta), con)
        if (length(setdiff(noncon, which(isFixed))) > 0) {
          nConstrained <- length(con)
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
    #toc(log = T)
    
    # Printing:
    #tic("Printing info")
    stop <- Sys.time()
    duration <- as.numeric(stop - start, units = "secs")
    duration <- sprintf("%02d:%02d:%06.3f", duration %/% 3600, (duration %% 3600) %/% 60, duration %% 60)
    if (iter == 1) {
      ansi <- interactive() && !is.null(Sys.getenv("TERM"))
      if (ansi) {
        cat("\033[1m", sprintf("%-11s %-15s %-20s", "Iteration", "LogLik", "Duration (hh:mm:ss)"), "\033[0m", "\n", sep = "")
      } else {
        cat(sprintf("%-11s %-15s %-20s \n", "Iteration", "LogLik", "Duration (hh:mm:ss)"))
      }
    }
    if (constraintsUsed) {
      msg <- paste0("(", paste(c(
        if (constraintsUsed) sprintf("%d constrained", nConstrained)
      ), collapse = " and "), ")")
      cat(sprintf("%-11d %-15.3f %-20s %-5s %-s\n", iter, LL[iter], duration, "", msg), sep = "")
    } else {
      cat(sprintf("%-11d %-15.3f %-20s\n", iter, LL[iter], duration))
    }
    #toc(log = T)
    
    # Checking for convergence:
    if (iter > 1) {
      deltaLL <- LL[iter] - LL[iter - 1]
      if (deltaLL > 0 & deltaLL <= tolConvLL & !(constrained | hitBoundary)) {
        LL <- LL[1:(iter + 1)]
        k <- k[, 1:(iter + 1)]
        LL[iter + 1] <- LL[iter]
        k[, iter + 1] <- k[, iter, drop = FALSE] + delta
        converged <- TRUE
      }
    }
    
    # Final update to parameters if not converged and not reached maxit:
    if (iter <= maxit & !(converged)) {
      k[, iter + 1] <- k[, iter, drop = FALSE] + delta
      iter <- iter + 1
    }
  }
  return(list(LL = LL, k = k, AI = AI, scores = scores, b = b, u = u, isFixed = isFixed, converged = converged))
}

