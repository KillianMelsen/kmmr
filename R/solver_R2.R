# A simpler wrapper for the Takahashi equations:
takahashi <- function(L) {
  # sparseinv_wrapper <- utils::getFromNamespace(".sparseinv_wrapper", "sparseinv")
  n <- nrow(L)
  d <- diag(L)
  L <- tril(L %*% Diagonal(x = 1 / d), -1)
  d <- d^2
  D <- Diagonal(x = d)
  dp <- diff(L@p)
  jj <- rep(seq_along(dp), dp)
  Zpattern <- sparseMatrix(c(L@i + 1, jj, 1:n), c(jj, L@i + 1, 1:n))
  Z <- forceSymmetric(sparseinv:::.sparseinv_wrapper(L, d, L, Zpattern))
  # Z <- forceSymmetric(sparseinv_wrapper(L, d, L, Zpattern))
  return(Z)
}

getCormat_corg <- function(p, params, doCormat = TRUE, doCormati = TRUE, doLogDet = TRUE) {
  if (p != (1 + sqrt(1 + 8 * length(params))) / 2) {
    stop("Number of factor levels and parameters for corg model do not match.")
  }
  cormat <- logDet <- cormati <- NULL
  if (doCormat) {
    cormat <- diag(p)
    cormat[lower.tri(cormat)] <- params
    cormat <- cormat + t(cormat) - diag(p)
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
    cormat <- diag(p)
    cormat[lower.tri(cormat)] <- params[1]
    cormat <- cormat + t(cormat) - diag(p)
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
    cormat <- diag(p)
  }
  if (doLogDet) {
    cormatChol <- diag(p)
    logDet <- 0
  }
  if (doCormati) {
    cormati <- diag(p)
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
    L <- matrix(params, nrow = p)
    m <- ncol(L)
    LLt <- Matrix::tcrossprod(L)
    cormat <- rbind(cbind(LLt, -L),
                    cbind(-t(L), diag(m)))
    diag(cormat) <- 1
  }
  if (doLogDet) {
    cormatChol <- chol(cormat)
    logDet <- 2 * sum(log(diag(cormatChol)))
  }
  if (doCormati) {
    PSIi <- 1 / (1 - diag(LLt)) # diag(solve(diag(p) - diag(p) * Matrix::tcrossprod(L)))
    PSIi_L <- PSIi * L
    cormati <- rbind(cbind(diag(PSIi), PSIi_L),
                     cbind(t(PSIi_L), diag(m) + Matrix::crossprod(L, PSIi_L)))
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
  deriv <- matrix(0, p, p)
  deriv[lower.tri(deriv)][cormatIndex] <- 1
  deriv <- deriv + t(deriv)
  return(deriv)
}

getDeriv_gk <- function(params, base) {
  deriv <- -base * exp(-params[1] * base)
  return(deriv)
}

getDeriv_fa <- function(cormatIndex, params, p) {
  m <- length(params) / p
  L <- Matrix(params, nrow = p, sparse = FALSE, doDiag = FALSE)
  dL <- Matrix(0, nrow = p, ncol = m, sparse = TRUE, doDiag = FALSE)
  dL[cormatIndex] <- 1
  # This is equal to dLL^T/dl + dPSI/dl, where dLL^T/dl = L %*% (dL/dl)^T + dL/dl %*% L^T and dPSI/dl = -I * dLL^t/dl, easier just to set the diagonal to 0...
  deriv <- tcrossprod(L, dL) + tcrossprod(dL, L) # L %*% t(dL) + dL %*% t(L)
  diag(deriv) <- 0
  return(deriv)
}

getDeriv_exfa <- function(cormatIndex, params, p, inverse) {
  L <- matrix(params, nrow = p)
  m <- ncol(L)
  dL <- matrix(0, nrow = p, ncol = m)
  dL[cormatIndex] <- 1
  if (!inverse) {
    # See getDeriv_fa, ignore dPSI/dl and just set the diagonal elements to 0:
    deriv <- rbind(cbind(Matrix::tcrossprod(L, dL) + Matrix::tcrossprod(dL, L), -dL),
                   cbind(-t(dL), matrix(0, m, m)))
    diag(deriv) <- 0
  } else { # If the derivative of the inverse extended factor analytic correlation matrix is needed:
    PSIi <- diag(1 / (1 - diag(Matrix::tcrossprod(L))))
    dPSI <- diag(-diag(Matrix::tcrossprod(L, dL) + Matrix::tcrossprod(dL, L)))
    dPSIi <- -PSIi %*% dPSI %*% PSIi
    # dPSIi <- diag(PSIi)^2 * dPSI
    deriv <- rbind(cbind(dPSIi, PSIi %*% dL + dPSIi %*% L),
                   cbind(t(PSIi %*% dL + dPSIi %*% L), t(dL) %*% PSIi %*% L + t(L) %*% dPSIi %*% L + t(L) %*% PSIi %*% dL))
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
#'
#' @returns A list containing the log-likelihood, parameter estimates, average information matrix,
#' scores, BLUEs for fixed effects, and BLUPs for random effects.
#' @keywords internal
solveLMM_R <- function(y, X, Z, S, init, weights, useWeights, maxit, stepSize, avInfReg, verbose, tolConvLL,
                       indices, u.applyvars, ur.types, u.nlevels, u.nparams, u.forFactor, u.indexCor,
                       bases, struc, i.var, i.nvar, isFixed, u.isStack, u.nparams.list, sparse) {
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
      C <- forceSymmetric(crossprod(W, Ri %*% W))
    } else {
      C <- as(forceSymmetric(crossprod(W, Ri %*% W)), "unpackedMatrix") # Store as a dense matrix since C will be dense after adding Gi.
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
    C.start <- nX + 1
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
          Gi <- kronecker(sigmai[[re]], nvarGi[[re]]) # This is a bit slow as well.
        } else {
          Gi <- sigmai[[re]]
        }

        # Now add to C:
        if (sparse) {
          # If we want C to be sparse but still have a dense random effect, we need to store the Gi for that
          # dense random effect in sparse format.
          if (!is(Gi, "sparseMatrix")) {
            Gi <- as(Gi, "CsparseMatrix")
          }
          # Padding Gi with zeros to make it conformable with the sparse base C.
          Gi <- sparseMatrix(i = Gi@i + 1 + C.start - 1,
                             j = rep(seq_len(ncol(Gi)), diff(Gi@p)) + C.start - 1,
                             x = Gi@x,
                             dims = dim(C),
                             symmetric = TRUE)
          C <- C + Gi
        } else {
          dimGi <- prod(u.nlevels[[re]])
          C[C.start:(C.start - 1 + dimGi), C.start:(C.start - 1 + dimGi)] <- C[C.start:(C.start - 1 + dimGi), C.start:(C.start - 1 + dimGi)] + Gi
          C.start <- C.start + dimGi
        }
      } else { # If there is no stack in the random effect
        # Always the outer product of standard deviations first:
        vars[[re]] <- tcrossprod(sqrt(k[indices[[re]][[1]], iter, drop = FALSE]))
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
            u.nlevels.applyvar <- u.nlevels[[re]][u.applyvars[[re]]] # Numbers of levels for (variable, obviously) factors that have variances applied to them
            u.nlevels.napplyvar <- u.nlevels.var[[re]][-u.applyvars[[re]]] # Number of levels for variable factors that do not have variances applied to them
            logDetvarGtotal[re] <- sum(log(diag(vars[[re]])))
            for (i in u.applyvars[[re]]) {
              logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * prod(u.nlevels.applyvar[-i]) # prod(numeric(0)) = 1
            }
            logDetvarGtotal[re] <- logDetvarGtotal[re] * prod(u.nlevels.napplyvar)
            # Now we can form the full log determinant:
            for (i in i.var[[re]][-u.applyvars[[re]]]) {
              logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * (dim.var[re] / u.nlevels.var[[re]][i])
            }
          }
        } else { # if there are no variable factors.
          sigmai[[re]] <- 1 / vars[[re]]
          logDetvarGtotal[re] <- -1 * log(sigmai[[re]])@x
        }
        if (n.nvar[re] > 0) {
          Gi <- kronecker(sigmai[[re]], nvarGi[[re]]) # This is a bit slow as well.
        } else {
          Gi <- sigmai[[re]]
        }
        
        # Now add to C:
        if (sparse) {
          # If we want C to be sparse but still have a dense random effect, we need to store the Gi for that
          # dense random effect in sparse format.
          if (!is(Gi, "sparseMatrix")) {
            Gi <- as(Gi, "CsparseMatrix")
          }
          # Padding Gi with zeros to make it conformable with the sparse base C.
          Gi <- sparseMatrix(i = Gi@i + 1 + C.start - 1,
                             j = rep(seq_len(ncol(Gi)), diff(Gi@p)) + C.start - 1,
                             x = Gi@x,
                             dims = dim(C),
                             symmetric = TRUE)
          C <- C + Gi
        } else {
          dimGi <- prod(u.nlevels[[re]])
          C[C.start:(C.start - 1 + dimGi), C.start:(C.start - 1 + dimGi)] <- C[C.start:(C.start - 1 + dimGi), C.start:(C.start - 1 + dimGi)] + Gi # This one is also slow
          C.start <- C.start + dimGi
        }
      }
    } # End loop over random effects.
    #toc(log = T)
    
    ### Cholesky decomposition of C ====
    #tic(msg = "Cholesky decomp of C")
    if (sparse) {
      Cchol <- Cholesky(C, LDL = FALSE, perm = TRUE, super = FALSE)
      # if (iter == 1) {
      #   Cchol <- Cholesky(C, LDL = FALSE, perm = TRUE, super = FALSE)
      # } else {
      #   # Update the cholesky factor
      # }
    } else {
      Cchol <- chol(C) # Good for the GPU?
    }
    #toc(log = T)
    
    ### BLUEs/BLUPs ====
    #tic("Solving for BLUEs and BLUPs")
    if (sparse) {
      # Apply permutation to the right-hand side of the MME: P %*% t(W) %*% Ri %*% y:
      permRHS <- solve(Cchol, t(W) %*% Ri %*% y, system = "P")
      # Sparse triangular forward solve: L %*% z = permRHS where L = Cchol:
      z <- solve(Cchol, permRHS, system = "L")
      # Sparse triangular backward solve: t(L) %*% bu = z where l = Cchol to obtain BLUEs and BLUPs:
      bu <- solve(Cchol, z, system = "Lt")
      # Final reverse permutation t(P) %*% bu to get the BLUEs and BLUPs into model order:
      bu <- solve(Cchol, bu, system = "Pt")
    } else {
      RHS <- t(W) %*% Ri %*% y
      z <- Matrix(backsolve(as.matrix(Cchol), RHS, transpose = TRUE))
      bu <- Matrix(backsolve(as.matrix(Cchol), z))
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
      logDetC <- determinant(Cchol, sqrt = FALSE, logarithm = TRUE)$modulus
    } else {
      logDetC <- 2 * sum(log(diag(Cchol)))
    }
    
    # logDetG:
    logDetG <- 0
    for (re in 1:nU) { # Start loop over random effects
      if (n.nvar[re] > 0) {
        logDetG <- logDetG + dim.nvar[re] * logDetvarGtotal[re] + dim.var[re] * logDetnvarGtotal[re]
      } else {
        logDetG <- logDetG + dim.nvar[re] * logDetvarGtotal[re]
      }
    } # End loop over random effects
    
    # logDetR:
    if (useWeights) {
      # For now just assume a diagonal structure:
      # logDetR <- sum(log(diag(Ri))) * -1 # diag() doesn't work with dgCMatrix...
      logDetR <- sum(log(Ri@x)) * -1
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
      dim.total <- dim.var[re] * dim.nvar[re]
      if (n.nvar[re] > 0) {
        umat <- Matrix(u[u.start:(u.start - 1 + dim.total), 1], nrow = dim.nvar[re], ncol = dim.var[re])
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
          if (isFixed[param] == 1) {
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
          if (isFixed[param] == 1) {
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
                A <- Matrix(0, prod(u.nlevels[[re]][u.applyvars[[re]]]), prod(u.nlevels[[re]][u.applyvars[[re]]]), sparse = FALSE, doDiag = FALSE)
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
    if (isFixed[nrow(k)] == 1) {
      # Wu[, nrow(k)] <- 0
    } else {
      Wu[, nrow(k)] <- (1 / k[nrow(k), iter]) * e@x # I think this is always true for diagonal R (also for non-diagonal fixed R, if that would ever be the case), even with weights.
    }
    #toc(log = T)
    
    ### Computing the AI matrix ====
    #tic("Computing the AI matrix and its inverse")
    WiWj <- forceSymmetric(crossprod(Wu, Ri %*% Wu)) # Similar to t(y) %*% Ri %*% y
    if (sparse) {
      XWjxZWj <- solve(Cchol, crossprod(W, Ri %*% Wu), system = "P") # Similar to permRHS
      zWu <- solve(Cchol, XWjxZWj, system = "L") # Similar to z
      buWu <- solve(Cchol, zWu, system = "Lt") # Similar to obtaining bu using a backsolve
    } else {
      XWjxZWj <- crossprod(W, Ri %*% Wu) # t(W) %*% Ri %*% Wu
      zWu <- backsolve(Cchol, XWjxZWj, transpose = TRUE)
      buWu <- Matrix(backsolve(Cchol, zWu))
    }
    buWuXWjxZWj <- forceSymmetric(crossprod(buWu, XWjxZWj))
    AI <- (WiWj - buWuXWjxZWj) / 2
    if (avInfReg[iter] > 1) { # Potential inflation of diagonal elements of the AI matrix
      diag(AI) <- diag(AI) * avInfReg[iter]
      # diag(AI)[ur.types == "L"] <- diag(AI)[ur.types == "L"] * avInfReg[iter] # Alternative that only inflates loadings.
    }
    if (any(isFixed == 1)) { # If there are any fixed parameters, remove them from the AI matrix before inverting.
      AIi <- Matrix(0, nrow(AI), ncol(AI), sparse = FALSE, doDiag = FALSE)
      notFixed <- which(isFixed == 0)
      AIi[notFixed, notFixed] <- solve(AI[notFixed, notFixed])
    } else {
      AIi <- solve(AI)
    }
    #toc(log = T)
    
    ### Computing the inverse of the coefficient matrix ====
    #tic("Inversion of C from cholesky factor")
    if (sparse) {
      # Sparse inverse subset and changing it back to model order:
      inv_Perm <- invPerm(Cchol@perm, zero.p = TRUE, zero.res = FALSE)
      Ci <- takahashi(as(Cchol, "dtCMatrix"))[inv_Perm, inv_Perm]
    } else {
      Ci <- chol2inv(Cchol)
    }
    #toc(log = TRUE)
    
    ### Computing scores ====
    #tic("Computing genetic scores")
    # See equations (18) and (15) from Jensen, Mantysaari, Madsen, & Thompson 1997:
    # "Residual maximum likelihood estimation of (co) variance components in multivariate
    # mixed linear models using average information"
    #### Genetic scores ====
    scores <- Matrix(0, nrow(k), 1, sparse = FALSE)
    corner <- nX + 1
    param <- 1
    for (re in 1:nU) { # Start of loop over random effects
      if (n.nvar[re] > 0) {
        traces <- Matrix(0, nrow(sigmai[[re]]), ncol(sigmai[[re]]), sparse = FALSE, doDiag = FALSE)
        
        for (i in 1:nrow(sigmai[[re]])) {
          for (j in 1:ncol(sigmai[[re]])) {
            if (j >= i) { # Only do the upper triangle
              rstart <- (corner + (i - 1) * dim.nvar[re])
              rend <- rstart - 1 + dim.nvar[re]
              cstart <- (corner + (j - 1) * dim.nvar[re])
              cend <- cstart - 1 + dim.nvar[re]
              # I suspect this can be done faster for a sparse Ci by selecting the relevant (nonzero) elements from the current block of Ci
              # and the corresonding elements of nvarGi. But for now let's keep the simple subsetting and hadamard product.
              traces[i, j] <- sum(nvarGi[[re]] * Ci[rstart:rend, cstart:cend])
            }
          }
        }
        corner <- corner + dim.var[re] * dim.nvar[re]
        if (nrow(traces) > 1) {
          traces <- forceSymmetric(traces)
        }
      } else {
        traces <- Ci[corner:(corner - 1 + dim.var[re]), corner:(corner - 1 + dim.var[re])]
        corner <- corner + dim.var[re]
      }
      
      if (n.nvar[re] > 0) {
        basisScores <- -0.5 * (dim.nvar[re] * sigmai[[re]] - sigmai[[re]] %*% traces %*% sigmai[[re]] - crossprod(umat_sigmai[[re]], nvarGi[[re]] %*% umat_sigmai[[re]])) # Last term maybe on the GPU?
      } else {
        basisScores <- -0.5 * (sigmai[[re]] - sigmai[[re]] %*% traces %*% sigmai[[re]] - crossprod(umat_sigmai[[re]])) # Last term maybe on the GPU?
      }
      for (i in 1:length(derivs[[re]])) {
        if (isFixed[param] == 0) {
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
    if (isFixed[nrow(k)] == 0) {
      WSW <- t(W) %*% S %*% W
      scores[nrow(scores), 1] <- as.numeric(-0.5 * (sum(S * evari) - sum(evari^2 * Ci[cbind(WSW@i + 1, rep(seq_len(ncol(WSW)), diff(WSW@p)))] * WSW@x) - (evari^2 * t(e) %*% S %*% e)))
    }
    #toc(log = T)
    
    ### Updating parameters ====
    #tic("Updating parameters and checking constraints")
    delta <- as.matrix(AIi %*% scores * stepSize[iter])
    
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
    loadings <- which(ur.types == "L") # Which of the parameters in k are bandwidths?
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
      bound.vars <- c() # Vector with indices of variances that were constrained in the current cycle due to the deltas pushing them across boundaries.
      bound.cors <- c() # Vector with indices of correlations that were constrained in the current cycle due to the deltas pushing them across boundaries.
      bound.bws <- c() # Vector with indices of bandwidths that were constrained in the current cycle due to the deltas pushing them across boundaries.
      npd.cors <- c() # Vector with indices of correlations that were constrained in the current cycle due to the presence of non-PD corg matrices.
      
      # Computing the proposed update so we can check for constraints:
      kUpdate <- as.matrix(k[, iter, drop = FALSE]) + delta
      
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
      kUpdate <- as.matrix(k[, iter, drop = FALSE]) + delta
      
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
      
      # Recomputing the proposed update so we can check whether all corg matrices are PD:
      kUpdate <- as.matrix(k[, iter, drop = FALSE]) + delta
      
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
      
      # Old code, which still assumed only a single random effect:
      # if (any(struc == "corg")) {
      #   for (i in which(struc == "corg")) {
      #     cormats[[i]] <- diag(nrow(cormats[[i]]))
      #     cormats[[i]][lower.tri(cormats[[i]])] <- kUpdate[indices[[1 + i]], 1]
      #     cormats[[i]] <- cormats[[i]] + t(cormats[[i]]) - diag(nrow(cormats[[i]]))
      #     ev <- eigen(cormats[[i]], only.values = TRUE)$values
      #     if (min(ev) <= 0) {
      #       nonPD <- TRUE
      #       checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
      #       constraintsUsed <- TRUE
      #       npd <- indices[[1 + i]]
      #       npd.cors <- c(npd.cors, npd)
      #       cormats[[i]] <- as.matrix(nearPD(cormats[[i]], corr = TRUE, posd.tol = 1e-3)$mat)
      #       delta[npd, 1] <- cormats[[i]][lower.tri(cormats[[i]])] - k[npd, iter]
      #     }
      #   }
      # }
      
      # Now that we have (1) constrained for excessive changes, (2) constrained for boundary
      # issues, and (3) constrained non-PD correlation matrices, we can recompute delta for any
      # unconstrained parameters:
      AI <- as.matrix(AI)
      scores <- as.matrix(scores)
      if (constrained | hitBoundary | nonPD) {
        # con consists of indices of params that were constrained in the previous cycles, with the new constrained params added.
        con <- sort(unique(c(con, con.vars, bound.vars, con.cors, bound.cors, npd.cors, con.bws, bound.bws)))
        noncon <- setdiff(1:nrow(delta), con)
        if (length(setdiff(noncon, which(isFixed == 1))) > 0) {
          nConstrained <- length(con)
          AIuu <- AI[noncon, noncon, drop = FALSE]
          uuisFixed <- isFixed[noncon]
          AIff <- AI[con, con, drop = FALSE]
          AIuf <- AI[noncon, con, drop = FALSE]
          scoresuu <- scores[noncon, 1, drop = FALSE]
          scoresff <- scores[con, 1, drop = FALSE]
          deltaff <- delta[con, 1, drop = FALSE]
          if (any(uuisFixed == 1)) { # If there are any fixed parameters in the unconstrained part, remove them from the AIuu matrix before inverting.
            AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
            uunotFixed <- which(uuisFixed == 0)
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
    if (any(boundary == 5)) {
      isFixed[which(boundary == 5)] <- 1
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

