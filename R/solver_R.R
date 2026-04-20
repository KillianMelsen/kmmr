#' formCormatChol
#'
#' @param struc Character indicating which correlation structure must be formed from the parameters.
#' @param p Numeric indicating the number of levels of this factor.
#' @param params Numeric vector containing the parameters needed to form the correlation matrix.
#' @param base Optional argument containing a numeric matrix to be used in forming the correlation matrix (distance matrix for Gaussian kernel, correlation matrix or kinship for linear kernel).
#' @param doChol Logical indicating whether a cholesky decomposition of the correlation matrix should be computed (for obtaining the log determinant and inverse).
#'
#' @returns A list containing the formed correlation matrix as well as its upper triangular Cholesky decomposition.
#' @keywords internal
formCormatChol <- function(struc, p = NULL, params = NULL, base = NULL, doChol = TRUE) {
  if (struc == "corg") {
    if (p != (1 + sqrt(1 + 8 * length(params))) / 2) {
      stop("Number of factor levels and parameters for corg model do not match.")
    }
    cormat <- diag(p)
    cormat[lower.tri(cormat)] <- params
    cormat <- cormat + t(cormat) - diag(p)
  } else if (struc == "cors") {
    cormat <- diag(p)
    cormat[lower.tri(cormat)] <- params[1]
    cormat <- cormat + t(cormat) - diag(p)
  } else if (struc == "gk") {
    cormat <- exp(-params[1] * base)
  } else if (struc == "lk") {
    cormat <- base
  } else if (struc == "id") {
    cormat <- diag(p)
  }
  if (doChol) {
    cormatChol <- chol(cormat)
  } else {
    cormatChol <- NULL
  }
  return(list(cormat, cormatChol))
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
#' @param H Numeric matrix of weights, i.e., standard errors or numbers of samples.
#' @param useH Logical indicating whether to use weights or not.
#' @param maxit Numeric value indicating the maximum number of iterations for the AI algorithm.
#' @param stepsize Numeric vector indicating the scaling factor for parameter updates for each iteration.
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
#'
#' @returns A list containing the log-likelihood, parameter estimates, average information matrix,
#' scores, BLUEs for fixed effects, and BLUPs for random effects.
#' @keywords internal
#' @importFrom methods as
#' @importFrom Matrix nearPD
solveLMM_R <- function(y, X, Z, S, init, H, useH, maxit, stepsize, verbose, tolConvLL,
                       indices, u.applyvars, ur.types, u.nlevels, u.nparams, u.forFactor, u.indexCor,
                       bases, struc, i.var, i.nvar, isFixed, u.isStack, u.nparams.list) {
  ## Prelims ====
  k <- matrix(NA, length(init), maxit + 1) # Matrix of parameter estimates.
  k[, 1] <- init # Filling with initial values.
  boundary <- numeric(length(init)) # Vector to hold how many times each parameter has crossed the boundary.
  n <- nrow(y)
  
  nU <- length(Z) # Number of random effects
  nR <- length(S) # Number of residual effects
  nUR <- nU + nR # Number of random effects + residual
  nX <- ncol(X) # Number of fixed effects
  intercept <- sum(X[, 1]) == nrow(X) # Is there an intercept?
  
  # Assume for now we have just one random and residual effect:
  # Z <- Z[[1]]
  S <- S[[1]]
  # indices <- indices[[1]]
  # u.applyvars <- u.applyvars[[1]]
  # ur.types <- ur.types[[1]] # Already unlisted
  # u.nlevels <- u.nlevels[[1]]
  # u.forFactor <- u.forFactor[[1]]
  # u.indexCor <- u.indexCor[[1]]
  # bases <- bases[[1]]
  # struc <- struc[[1]]
  # i.var <- i.var[[1]]
  # i.nvar <- i.nvar[[1]]
  # isFixed <- isFixed[[1]] # Already unlisted
  
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
    # n.var[re] <- length(i.var[[re]]) # How many variable factors are there?
    # n.nvar[re] <- length(i.nvar[[re]]) # How many non-variable factors are there?
    if (n.nvar[re] > 0) { # All of this can only be done if there are non-variable factors (such as a kinship/linear kernel for genotypes).
      logDetnvarG <- numeric(n.nvar[re]) # A log determinant for each non-variable base matrix (linear kernel).
      for (i in i.nvar[[re]]) {
        logDetnvarG[i - n.var[re]] <- -1 * sum(log(eigen(bases[[re]][[i]])$values)) # Times -1 because the non-variable base matrices (linear kernels) are in inverse form. # Maybe on the GPU for large matrices?
      }
      
      # If the non-variable part of G consists of more than 1 factor we need to compute the total log determinant:
      if (n.nvar[re] > 1) {
        # u.levels.nvar <- u.nlevels[[re]][i.nvar[[re]]]
        logDetnvarGtotal[re] <- 0
        for (i in 1:n.nvar[re]) {
          # logDetnvarGtotal[re] <- logDetnvarGtotal[re] + logDetnvarG[i] * prod(u.levels.nvar[-i])
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
    if (useH) {
      # Add the weights 
    } else {
      Ri <- 1 / k[nrow(k), iter] * S # The residual variance is always last in kappa for now
    }
    #toc(log = T)
    
    ### Forming W = [X Z y] ====
    #tic("Forming Wy and M = t(Wy) %*% Ri %*% Wy")
    if (iter == 1) {
      W <- X
      for (re in 1:nU) {
        W <- cbind(W, Z[[re]])
      }
      Wy <- cbind(W, y)
    }
    
    ### Computing M = t(Wy) %*% Ri %*% Wy
    if (useH) {
      # We need to compute the entire thing every time? Probably not for diagonal R...
    } else {
      if (iter == 1) {
        M0 <- Matrix::t(Wy) %*% Wy # Slow, but very sparse so faster in C++, maybe use something like syrk()?
      }
      M <- as.matrix(M0 * (1 / k[nrow(k), iter]))
    }
    #toc(log = T)
    
    ### Adding Gi to M ====
    #tic("Computing Gi")
    # cormatis <- vector("list", nU)
    cormats <- vector("list", nU)
    sigmai <- vector("list", nU)
    logDetvarGtotal <- numeric(nU)
    vars <- vector("list", nU)
    M.start <- nX + 1
    for (re in 1:nU) { # Start loop over random effects
      # Computing Gi every iteration:
      if (any(u.isStack[[re]])) { # If there is a stack in the random effect
        # Determining which factor is a stack and the size of the stack:
        i.stack <- which(u.isStack[[re]])
        sos <- length(u.nparams.list[[re]][[i.stack + 1]]) # + 1 due to the variances
        # Always the outer products of standard deviations for each stack first:
        vars[[re]] <- vector("list", sos)
        nvars <- length(indices[[re]][[1]]) / sos
        for (s in 1:sos) {
          vars[[re]][[s]] <- outer(sqrt(k[indices[[re]][[1]][(s * nvars - nvars + 1):(s * nvars)], iter]),
                                   sqrt(k[indices[[re]][[1]][(s * nvars - nvars + 1):(s * nvars)], iter]))
        }
        if (n.var[re] > 0) { # This should always be the case for a random effect with a stack
          cormatis <- vector("list", n.var[re]) # One inverse correlation matrix for every variable factor.
          cormats[[re]] <- vector("list", n.var[re]) # One correlation matrix for every variable factor.
          logDetvarG <- numeric(n.var[re])
          for (i in i.var[[re]]) { # Loop to form the covariance stack and correlation matrices
            # The first variable factor has to be the one with a stack:
            if (i == i.stack) {
              # We need to calculate the complete weighed sum of the structures using the variances:
              covmat.stack <- matrix(0, u.nlevels.var[[re]][i.var[[re]][i]], u.nlevels.var[[re]][i.var[[re]][i]])
              cormats[[re]][[i]] <- vector("list", sos)
              indices.start <- 1
              for (s in 1:sos) {
                indices.stack <- (indices.start):(indices.start - 1 + u.nparams.list[[re]][[1 + i]][s])
                result <- formCormatChol(struc = struc[[re]][[i]][s],
                                         p = u.nlevels[[re]][i],
                                         params = k[indices[[re]][[1 + i]][indices.stack], iter],
                                         base = bases[[re]][[i]][[s]],
                                         doChol = FALSE)
                indices.start <- indices.start + u.nparams.list[[re]][[1 + i]][s]
                cormats[[re]][[i]][[s]] <- result[[1]] # Correlation matrix
                covmat.stack <- covmat.stack + vars[[re]][[s]] * result[[1]] # We first construct sigma by summing all stack layers, then invert
              }
              cormatis[[i]] <- chol(covmat.stack)
              logDetvarG[i] <- 2 * sum(log(diag(cormatis[[i]])))
              cormatis[[i]] <- chol2inv(cormatis[[i]])
            } else {
              # Do the remaining variable structures that are not part of the stack (there can only be a single stack in a random effect):
              result <- formCormatChol(struc = struc[[re]][[i]],
                                       p = u.nlevels[[re]][i],
                                       params = k[indices[[re]][[1 + i]], iter],
                                       base = bases[[re]][[i]],
                                       doChol = TRUE)
              cormats[[re]][[i]] <- result[[1]] # Correlation matrix
              cormatis[[i]] <- result[[2]] # Cholesky factor of correlation matrix
              logDetvarG[i] <- 2 * sum(log(diag(cormatis[[i]]))) # Compute log determinants of all correlation matrices from cholesky factors
              cormatis[[i]] <- chol2inv(cormatis[[i]]) # Invert each correlation matrix using cholesky factor.
            }
          }# End of loop to form the covariance stack and variable correlation matrices
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
        } # End of if (n.var[re] > 0), which should always be TRUE for random effects with a stack
        if (n.nvar[re] > 0) {
          Gi <- kronecker(sigmai[[re]], nvarGi[[re]]) # This is a bit slow as well.
        } else {
          Gi <- sigmai[[re]]
        }

        # Now add to M:
        dimGi <- prod(u.nlevels[[re]])
        M[M.start:(M.start - 1 + dimGi), M.start:(M.start - 1 + dimGi)] <- M[M.start:(M.start - 1 + dimGi), M.start:(M.start - 1 + dimGi)] + Gi # This one is also slow
        M.start <- M.start + dimGi
        
      } else { # If there is no stack in the random effect
        # Always the outer product of standard deviations first:
        vars[[re]] <- outer(sqrt(k[indices[[re]][[1]], iter]), sqrt(k[indices[[re]][[1]], iter]))
        if (n.var[re] > 0) { # Only form the correlation matrices for variable factors if there are variable factors.
          cormatis <- vector("list", n.var[re]) # One inverse correlation matrix for every variable factor.
          cormats[[re]] <- vector("list", n.var[re]) # One correlation matrix for every variable factor.
          logDetvarG <- numeric(n.var[re])
          for (i in i.var[[re]]) { # Loop to form the correlations matrices
            result <- formCormatChol(struc = struc[[re]][i], # Can we directly form the cholesky factor from the parameters and, if applicable, base matrices?
                                     p = u.nlevels[[re]][i],
                                     params = k[indices[[re]][[1 + i]], iter],
                                     base = bases[[re]][[i]],
                                     doChol = TRUE)
            cormats[[re]][[i]] <- result[[1]] # Correlation matrix
            cormatis[[i]] <- result[[2]] # Cholesky factor of correlation matrix
            logDetvarG[i] <- 2 * sum(log(diag(cormatis[[i]]))) # Compute log determinants of all correlation matrices from cholesky factors
            cormatis[[i]] <- chol2inv(cormatis[[i]]) # Invert each correlation matrix using cholesky factor.
          } # End of loop to form the correlations matrices
          if (u.applyvars[[re]][1] == 0) { # If there is just a single variance.
            sigmai[[re]] <- (1 / as.numeric(vars[[re]])) * Reduce(kronecker, cormatis) # Form the inverse covariance matrix.
            # Form the log determinant of sigma:
            # u.nlevels.var <- u.nlevels[[re]][i.var[[re]]]
            if (n.var[re] > 1) {
              logDetvarGtotal[re] <- 0
              for (i in i.var[[re]]) {
                # logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * prod(u.nlevels.var[-i])
                logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * (dim.var[re] / u.nlevels.var[[re]][i])
              }
            } else {
              logDetvarGtotal[re] <- logDetvarG
            }
            # logDetvarGtotal[re] <- prod(u.nlevels.var) * log(as.numeric(vars)) + logDetvarGtotal[re]
            logDetvarGtotal[re] <- dim.var[re] * log(as.numeric(vars[[re]])) + logDetvarGtotal[re]
          } else if (all(u.applyvars[[re]] == i.var[[re]])) { # If variances are applied to all variable factors.
            sigmai[[re]] <- (1 / vars[[re]]) * Reduce(kronecker, cormatis) # Form the inverse covariance matrix.
            # Form the log determinant of sigma:
            # u.nlevels.var <- u.nlevels[[re]][i.var[[re]]]
            logDetvarGtotal[re] <- sum(log(diag(vars[[re]])))
            for (i in i.var[[re]]) {
              # logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * prod(u.nlevels.var[-i])
              logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * (dim.var[re] / u.nlevels.var[[re]][i])
            }
          } else { # If variances are applied to 1 or more factors, but not all factors.
            sigmai[[re]] <- kronecker((1 / vars[[re]]) * Reduce(kronecker, cormatis[u.applyvars[[re]]]), Reduce(kronecker, cormatis[-u.applyvars[[re]]])) # Form the inverse covariance matrix.
            # Form the log determinant of sigma:
            # First form the log determinant of the part with variances:
            # u.nlevels.var <- u.nlevels[[re]][i.var[[re]]]
            u.nlevels.applyvar <- u.nlevels[[re]][u.applyvars[[re]]] # Numbers of levels for (variable, obviously) factors that have variances applied to them
            u.nlevels.napplyvar <- u.nlevels.var[[re]][-u.applyvars[[re]]] # Number of levels for variable factors that do not have variances applied to them
            logDetvarGtotal[re] <- sum(log(diag(vars[[re]])))
            for (i in u.applyvars[[re]]) {
              logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * prod(u.nlevels.applyvar[-i]) # prod(numeric(0)) = 1
            }
            logDetvarGtotal[re] <- logDetvarGtotal[re] * prod(u.nlevels.napplyvar)
            # Now we can form the full log determinant:
            for (i in i.var[[re]][-u.applyvars[[re]]]) {
              # logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * prod(u.nlevels.var[-i])
              logDetvarGtotal[re] <- logDetvarGtotal[re] + logDetvarG[i] * (dim.var[re] / u.nlevels.var[[re]][i])
            }
          }
        } else { # if there are no variable factors.
          sigmai[[re]] <- 1 / vars[[re]]
          logDetvarGtotal[re] <- -1 * as.numeric(log(sigmai[[re]]))
        }
        if (n.nvar[re] > 0) {
          Gi <- kronecker(sigmai[[re]], nvarGi[[re]]) # This is a bit slow as well.
        } else {
          Gi <- sigmai[[re]]
        }
        
        # Now add to M:
        dimGi <- prod(u.nlevels[[re]])
        M[M.start:(M.start - 1 + dimGi), M.start:(M.start - 1 + dimGi)] <- M[M.start:(M.start - 1 + dimGi), M.start:(M.start - 1 + dimGi)] + Gi # This one is also slow
        M.start <- M.start + dimGi
      }
    } # End loop over random effects.
    #toc(log = T)
    
    # Temporary code to verify whether the fit from dev/dev_stacked_models_lkGen_plus_stackEnv-gk-corg_by_corgMan_by_lkGen.RData is processed correctly:
    # G1 <- 3.0 * K
    # G2 <- kronecker(kronecker(outer(sqrt(c(1.1, 1.2, 1.3)), sqrt(c(1.1, 1.2, 1.3))) * exp(-0.1*ED) + outer(sqrt(c(2.1, 2.2, 2.3)), sqrt(c(2.1, 2.2, 2.3))) * matrix(c(1.0, 0.4, 0.5, 0.4, 1.0, 0.6, 0.5, 0.6, 1.0), 3, 3), matrix(c(1.0, 0.8, 0.8, 1.0), 2, 2)), K)
    # all.equal(Gi, solve(G2))
    # all.equal(logDetvarGtotal[1], sum(log(eigen(matrix(3.0))$values)))
    # all.equal(logDetvarGtotal[2], sum(log(eigen(kronecker(outer(sqrt(c(1.1, 1.2, 1.3)), sqrt(c(1.1, 1.2, 1.3))) * exp(-0.1*ED) + outer(sqrt(c(2.1, 2.2, 2.3)), sqrt(c(2.1, 2.2, 2.3))) * matrix(c(1.0, 0.4, 0.5, 0.4, 1.0, 0.6, 0.5, 0.6, 1.0), 3, 3), matrix(c(1.0, 0.8, 0.8, 1.0), 2, 2)))$values)))
    # Zmat1 <- as.matrix(Z[[1]])
    # Zmat2 <- as.matrix(Z[[2]])
    # Zmat <- cbind(Zmat1, Zmat2)
    # Xmat <- as.matrix(X)
    # Cmat <- rbind(cbind(t(Xmat) %*% as.matrix(Ri) %*% Xmat,    t(Xmat) %*% as.matrix(Ri) %*% Zmat1,                t(Xmat) %*% as.matrix(Ri) %*% Zmat2),
    #               cbind(t(Zmat1) %*% as.matrix(Ri) %*% Xmat,   t(Zmat1) %*% as.matrix(Ri) %*% Zmat1 + solve(G1),   t(Zmat1) %*% as.matrix(Ri) %*% Zmat2),
    #               cbind(t(Zmat2) %*% as.matrix(Ri) %*% Xmat,   t(Zmat2) %*% as.matrix(Ri) %*% Zmat1,               t(Zmat2) %*% as.matrix(Ri) %*% Zmat2 + solve(G2)))
    # all.equal(Cmat, M[-nrow(M), -ncol(M)])
    # G <- rbind(cbind(G1, matrix(0, nrow(G1), ncol(G2))),
    #            cbind(matrix(0, nrow(G2), ncol(G1)), G2))
    # R <- solve(as.matrix(Ri))
    # V <- Zmat %*% G %*% t(Zmat) + R
    # Vi <- solve(V)
    # P <- Vi - Vi %*% Xmat %*% solve(t(Xmat) %*% Vi %*% Xmat) %*% t(Xmat) %*% Vi
    # G2.1 <- kronecker(kronecker(outer(sqrt(c(1.1, 1.2, 1.3)), sqrt(c(1.1, 1.2, 1.3))) * exp(-0.1*ED), matrix(c(1.0, 0.8, 0.8, 1.0), 2, 2)), K)
    # G2.2 <- kronecker(kronecker(outer(sqrt(c(2.1, 2.2, 2.3)), sqrt(c(2.1, 2.2, 2.3))) * matrix(c(1.0, 0.4, 0.5, 0.4, 1.0, 0.6, 0.5, 0.6, 1.0), 3, 3), matrix(c(1.0, 0.8, 0.8, 1.0), 2, 2)), K)
    # Znaive <- cbind(as.matrix(Z[[1]]), as.matrix(Z[[2]]), as.matrix(Z[[2]]))
    # Gnaive <- rbind(cbind(G1, matrix(0, nrow(G1), ncol(G2.1)), matrix(0, nrow(G1), ncol(G2.2))),
    #                 cbind(matrix(0, nrow(G2.1), ncol(G1)), G2.1, matrix(0, nrow(G2.1), ncol(G2.2))),
    #                 cbind(matrix(0, nrow(G2.2), ncol(G1)), matrix(0, nrow(G2.2), ncol(G2.1)), G2.2))
    # Vnaive <- Znaive %*% Gnaive %*% t(Znaive) + R
    # Vnaivei <- solve(Vnaive)
    # Pnaive <- Vnaivei - Vnaivei %*% Xmat %*% solve(t(Xmat) %*% Vnaivei %*% Xmat) %*% t(Xmat) %*% Vnaivei
    # all.equal(yPy, as.numeric(t(y) %*% P %*% y))
    # all.equal(yPy, as.numeric(t(y) %*% Pnaive %*% y))
    # all.equal(P, Pnaive)
    # all.equal(Vi, Vnaivei)
    # dim(G); dim(Gnaive)
    
    ### Cholesky decomposition of M ====
    #tic(msg = "Cholesky decomp of M")
    Mchol <- chol(M) # Upper triangular # This is good for the GPU
    #toc(log = T)
    yPy <- Mchol[nrow(Mchol), ncol(Mchol)]^2
    Mchol_XZ <- Mchol[1:(nrow(Mchol) - 1), 1:(ncol(Mchol) - 1)]
    My <- Mchol[1:(nrow(Mchol) - 1), ncol(Mchol), drop = FALSE]
    
    ### Computing the log-likelihood ====
    #tic("Computing log-likelihood")
    # Computing logDetC:
    logDetC <- 2 * sum(log(diag(Mchol_XZ)))
    
    # Computing logDetG:
    logDetG <- 0
    for (re in 1:nU) { # Start loop over random effects
      # dim.var <- prod(u.nlevels[[re]][i.var[[re]]])
      # dim.nvar <- prod(u.nlevels[[re]][i.nvar[[re]]])
      if (n.nvar[re] > 0) {
        logDetG <- logDetG + dim.nvar[re] * logDetvarGtotal[re] + dim.var[re] * logDetnvarGtotal[re]
      } else {
        logDetG <- logDetG + dim.nvar[re] * logDetvarGtotal[re]
      }
    } # End loop over random effects
    
    # Computing logDetR:
    if (useH) {
      # Do something here.
    } else {
      logDetR <- n * log(k[nrow(k), iter])
    }
    
    # Full LL:
    LL[iter] <- -0.5 * (logDetC + logDetR + logDetG + yPy)
    #toc(log = T)
    
    ### Solving for BLUEs and BLUPs ====
    #tic("Solving for BLUEs and BLUPs")
    bu <- backsolve(Mchol_XZ, My)
    b <- bu[1:nX, 1, drop = FALSE]
    u <- bu[(nX + 1):nrow(bu), 1, drop = FALSE]
    #toc(log = T)
    
    ### Working variates ====
    #### Genetic working variates ====
    #tic("Genetic working variates")
    Wu <- matrix(NA, n, nrow(k)) # Matrix for the working variates.
    umat_sigmai <- vector("list", nU) # We will store these matrices for computing the scores later.
    derivs <- vector("list", nU) # Same thing for all derivatives.
    u.start <- 1
    param <- 1 # Keeping track of which column/parameter of Wu we are working on.
    for (re in 1:nU) { # Start of loop over random effects
      dim.total <- dim.var[re] * dim.nvar[re]
      if (n.nvar[re] > 0) {
        umat <- matrix(u[u.start:(u.start - 1 + dim.total), 1], nrow = dim.nvar[re], ncol = dim.var[re])
      } else {
        umat <- matrix(u[u.start:(u.start - 1 + dim.total), 1], nrow = 1, ncol = dim.var[re])
      }
      u.start <- u.start + dim.total
      umat_sigmai[[re]] <- umat %*% sigmai[[re]]
      derivs[[re]] <- vector("list", u.nparams[[re]])
      
      for (i in 1:u.nparams[[re]]) {
        if (isFixed[param] == 1) {
          Wu[, param] <- 0
          param <- param + 1
        } else {
          if (ur.types[param] == "V") {
            if (u.applyvars[[re]][1] == 0) { # If there is only a single variance.
              if (n.var[re] > 0) {
                derivs[[re]][[i]] <- Reduce(kronecker, cormats[[re]])
              } else {
                derivs[[re]][[i]] <- matrix(1)
              }
            } else { # If variances are applied to 1 or more factors.
              A <- matrix(0, prod(u.nlevels[[re]][u.applyvars[[re]]]), prod(u.nlevels[[re]][u.applyvars[[re]]]))
              A[i,] <- 0.5
              if (all(u.applyvars[[re]] == i.var[[re]])) { # If variances are applied to all variable factors.
                derivs[[re]][[i]] <- (A + t(A)) * vars[[re]] * Reduce(kronecker, cormats[[re]]) / k[param, iter]
              } else { # If variances are applied to 1 or more factors, but not all factors.
                derivs[[re]][[i]] <- kronecker((A + t(A)) * vars[[re]] * Reduce(kronecker, cormats[[re]][u.applyvars[[re]]]) / k[param, iter], Reduce(kronecker, cormats[[re]][-u.applyvars[[re]]]))
              }
            }
          } else if (ur.types[param] %in% c("R", "B")) {
            j <- u.forFactor[[re]][i] # Which factor does this correlation or bandwidth belong to?
            cormats1 <- cormats[[re]] # Making a copy of the correlation matrices so we can replace the matrix of factor j with the indicator matrix or derivative of the GK.
            if (ur.types[param] == "R") {
              icor <- u.indexCor[[re]][[i]] # This correlation is the icor-th correlation in factor j.
              B <- matrix(0, u.nlevels[[re]][j], u.nlevels[[re]][j]) # Making the indicator matrix.
              B[lower.tri(B)][icor] <- 1 # Filling the indicator matrix.
              B <- B + t(B)
              cormats1[[j]] <- B # Replacing with the indicator matrix.
            } else if (ur.types[param] == "B") {
              cormats1[[j]] <- -1 * bases[[re]][[j]] * cormats1[[j]] # d/dh exp(-h * D) = -D * exp(-h * D)
            }
            if (u.applyvars[[re]][1] == 0) {
              derivs[[re]][[i]] <- as.numeric(vars[[re]]) * Reduce(kronecker, cormats1)
            } else if (all(u.applyvars[[re]] == i.var[[re]])) {
              derivs[[re]][[i]] <- vars[[re]] * Reduce(kronecker, cormats1)
            } else {
              derivs[[re]][[i]] <- kronecker(vars[[re]] * Reduce(kronecker, cormats1[u.applyvars[[re]]]), Reduce(kronecker, cormats1[-u.applyvars[[re]]]))
            }
          }
          Wu[, param] <- as.matrix(Z[[re]] %*% matrix(umat_sigmai[[re]] %*% derivs[[re]][[i]]))
          param <- param + 1
        }
      }
    } # End of loop over random effects
    #toc(log = T)
    
    #### Residual working variates ====
    #tic("Residual working variates")
    e <- y - W %*% bu
    if (isFixed[nrow(k)] == 1) {
      Wu[, nrow(k)] <- 0
    } else {
      Wu[, nrow(k)] <- (1 / k[nrow(k), iter]) * e@x # I think this is always true for diagonal R, even with weights.
    }
    #toc(log = T)
    
    ### Computing the AI matrix ====
    #tic("Computing the AI matrix and its inverse")
    # XWjxZWj <- Matrix::t(W) %*% Ri %*% Wu
    XWjxZWj <- Matrix::t(W) %*% Wu * (1 / k[nrow(k), iter])
    # WiWj <- t(Wu) %*% Ri %*% Wu
    WiWj <- t(Wu) %*% Wu * (1 / k[nrow(k), iter])
    # We are trying to solve C %*% buWu = XWjxZWj
    # We already have R where R^T %*% R = C (cholesky decomposition of C)
    # So instead of :
    # buWu <- solve(M[1:(nrow(M) - 1), 1:(ncol(M) - 1)], XWjxZWj)
    # we do:
    buWu <- backsolve(Mchol_XZ, forwardsolve(t(Mchol_XZ), XWjxZWj)) # This is a bit slow too
    AI <- as((WiWj - (t(buWu) %*% XWjxZWj)) / 2, "matrix")
    if (any(isFixed == 1)) { # If there are any fixed parameters, remove them from the AI matrix before inverting.
      AIi <- matrix(0, nrow(AI), ncol(AI))
      notFixed <- which(isFixed == 0)
      AIi[notFixed, notFixed] <- solve(AI[notFixed, notFixed])
    } else {
      AIi <- solve(AI)
    }
    #toc(log = T)
    
    ### Computing the inverse of the coefficient matrix ====
    #tic("Inversion of C from cholesky factor")
    Ci <- chol2inv(Mchol_XZ) # Good point to use the GPU!
    #toc(log = TRUE)
    
    ### Computing scores ====
    #tic("Computing genetic scores")
    # See equations (18) and (15) from Jensen, Mantysaari, Madsen, & Thompson 1997:
    # "Residual maximum likelihood estimation of (co) variance components in multivariate
    # mixed linear models using average information"
    #### Genetic scores ====
    scores <- matrix(0, nrow(k), 1)
    corner <- nX + 1
    param <- 1
    for (re in 1:nU) { # Start of loop over random effects
      if (n.nvar[re] > 0) {
        traces <- matrix(0, nrow(sigmai[[re]]), ncol(sigmai[[re]]))
        
        for (i in 1:nrow(sigmai[[re]])) {
          for (j in 1:ncol(sigmai[[re]])) {
            if (j >= i) { # Only do the upper triangle
              rstart <- (corner + (i - 1) * dim.nvar[re])
              rend <- rstart - 1 + dim.nvar[re]
              cstart <- (corner + (j - 1) * dim.nvar[re])
              cend <- cstart - 1 + dim.nvar[re]
              traces[i, j] <- sum(nvarGi[[re]] * Ci[rstart:rend, cstart:cend])
            }
          }
        }
        corner <- corner + dim.var[re] * dim.nvar[re]
        if (nrow(traces) > 1) {
          traces <- traces + t(traces) - diag(diag(traces))
        }
      } else {
        traces <- Ci[corner:(corner - 1 + dim.var[re]), corner:(corner - 1 + dim.var[re])]
        corner <- corner + dim.var[re]
      }
      
      if (n.nvar[re] > 0) {
        basisScores <- -0.5 * (dim.nvar[re] * sigmai[[re]] - sigmai[[re]] %*% traces %*% sigmai[[re]] - t(umat_sigmai[[re]]) %*% nvarGi[[re]] %*% umat_sigmai[[re]]) # Last term maybe on the GPU?
      } else {
        basisScores <- -0.5 * (sigmai[[re]] - sigmai[[re]] %*% traces %*% sigmai[[re]] - t(umat_sigmai[[re]]) %*% umat_sigmai[[re]]) # Last term maybe on the GPU?
      }
      for (i in 1:length(derivs[[re]])) {
        if (isFixed[param] == 1) {
          scores[param, 1] <- 0
        } else {
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
    if (isFixed[nrow(k)] == 1) {
      scores[nrow(scores), 1] <- 0
    } else {
      WSW <- Matrix::t(W) %*% S %*% W
      scores[nrow(scores), 1] <- as.numeric(-0.5 * (sum(S * evari) - sum(evari^2 * Ci[cbind(WSW@i + 1, rep(seq_len(ncol(WSW)), diff(WSW@p)))] * WSW@x) - (evari^2 * Matrix::t(e) %*% S %*% e)))
    }
    #toc(log = T)
    
    ### Updating parameters ====
    #tic("Updating parameters and checking constraints")
    delta <- as.matrix(AIi %*% scores * stepsize[iter])
    
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
      bound.vars <- c() # Vector with indices of variances that were constrained in the current cycle due to the deltas pushing them across boundaries.
      bound.cors <- c() # Vector with indices of correlations that were constrained in the current cycle due to the deltas pushing them across boundaries.
      bound.bws <- c() # Vector with indices of bandwidths that were constrained in the current cycle due to the deltas pushing them across boundaries.
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
      
      # Recomputing the proposed update so we can check whether all corg matrices are PD:
      kUpdate <- k[, iter, drop = FALSE] + delta
      
      # Checking whether all unstructured correlation matrices are PD:
      if (any(struc == "corg")) {
        for (i in which(struc == "corg")) {
          cormats[[i]] <- diag(nrow(cormats[[i]]))
          cormats[[i]][lower.tri(cormats[[i]])] <- kUpdate[indices[[1 + i]], 1]
          cormats[[i]] <- cormats[[i]] + t(cormats[[i]]) - diag(nrow(cormats[[i]]))
          ev <- eigen(cormats[[i]], only.values = TRUE)$values
          if (min(ev) <= 0) {
            nonPD <- TRUE
            checkConstraints <- TRUE # A new parameter was constrained, so we need to update delta and check again.
            constraintsUsed <- TRUE
            npd <- indices[[1 + i]]
            npd.cors <- c(npd.cors, npd)
            cormats[[i]] <- as.matrix(nearPD(cormats[[i]], corr = TRUE, posd.tol = 1e-3)$mat)
            delta[npd, 1] <- cormats[[i]][lower.tri(cormats[[i]])] - k[npd, iter]
          }
        }
      }
      
      # Now that we have (1) constrained for excessive changes, (2) constrained for boundary
      # issues, and (3) constrained non-PD correlation matrices, we can recompute delta for any
      # unconstrained parameters:
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

