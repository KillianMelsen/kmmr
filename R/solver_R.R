#' formCormatChol
#'
#' @param struc Character indicating which correlation structure must be formed from the parameters.
#' @param p Numeric indicating the number of levels of this factor.
#' @param params Numeric vector containing the parameters needed to form the correlation matrix.
#' @param base Optional argument containing a numeric matrix to be used in forming the correlation matrix (distance matrix for Gaussian kernel, correlation matrix or kinship for linear kernel).
#'
#' @returns A list containing the formed correlation matrix as well as its upper triangular Cholesky decomposition.
#' @keywords internal
formCormatChol <- function(struc, p = NULL, params = NULL, base = NULL) {
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
  cormatChol <- chol(cormat)
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
#'
#' @returns A list containing the log-likelihood, parameter estimates, average information matrix,
#' scores, BLUEs for fixed effects, and BLUPs for random effects.
#' @keywords internal
#' @importFrom methods as
#' @importFrom Matrix nearPD
solveLMM_R <- function(y, X, Z, S, init, H, useH, maxit, stepsize, verbose, tolConvLL,
                       indices, u.applyvars, ur.types, u.nlevels, u.forFactor, u.indexCor,
                       bases, struc, i.var, i.nvar, isFixed) {
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
  Z <- Z[[1]]
  S <- S[[1]]
  indices <- indices[[1]]
  u.applyvars <- u.applyvars[[1]]
  # ur.types <- ur.types[[1]] # Already unlisted
  u.nlevels <- u.nlevels[[1]]
  u.forFactor <- u.forFactor[[1]]
  u.indexCor <- u.indexCor[[1]]
  bases <- bases[[1]]
  struc <- struc[[1]]
  i.var <- i.var[[1]]
  i.nvar <- i.nvar[[1]]
  # isFixed <- isFixed[[1]] # Already unlisted
  
  ## Computing determinants of matrices associated with non-variable terms before iterating ====
  n.var <- length(i.var) # How many variable factors are there?
  n.nvar <- length(i.nvar) # How many non-variable factors are there?
  
  if (n.nvar > 0) { # All of this can only be done if there are non-variable factors (such as a kinship/linear kernel for genotypes).
    logDetnvarG <- numeric(n.nvar) # A log determinant for each non-variable base matrix (linear kernel).
    for (i in i.nvar) {
      logDetnvarG[i - n.var] <- -1 * sum(log(eigen(bases[[i]])$values)) # Times -1 because the non-variable base matrices (linear kernels) are in inverse form. # Maybe on the GPU for large matrices?
    }
    
    # If the non-variable part of G consists of more than 1 factor we need to compute the total log determinant:
    if (n.nvar > 1) {
      u.levels.nvar <- u.nlevels[i.nvar]
      logDetnvarGtotal <- 0
      for (i in 1:n.nvar) {
        logDetnvarGtotal <- logDetnvarGtotal + logDetnvarG[i] * prod(u.levels.nvar[-i])
      }
    } else {
      logDetnvarGtotal <- logDetnvarG
    }
    
    # Form the non-variable part of Gi using kronecker products before iterating if it exists:
    nvarGi <- Reduce(kronecker, bases[i.nvar])
  }
  
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
      for (i in 1:nU) {
        W <- cbind(W, Z)
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
    # Computing Gi every iteration:
    #tic("Computing Gi")
    # Always the outer product of standard deviations first:
    vars <- outer(sqrt(k[indices[[1]], iter]), sqrt(k[indices[[1]], iter]))
    if (n.var > 0) { # Only form the correlation matrices for variable factors if there are variable factors.
      cormatis <- vector("list", n.var) # One inverse correlation matrix for every variable factor.
      cormats <- vector("list", n.var) # One correlation matrix for every variable factor.
      logDetvarG <- numeric(n.var)
      for (i in i.var) {
        result <- formCormatChol(struc = struc[i], # Can we directly form the cholesky factor from the parameters and, if applicable, base matrices?
                                 p = u.nlevels[i],
                                 params = k[indices[[1 + i]], iter],
                                 base = bases[[i]])
        cormats[[i]] <- result[[1]] # Correlation matrix
        cormatis[[i]] <- result[[2]] # Cholesky factor of correlation matrix
        logDetvarG[i] <- 2 * sum(log(diag(cormatis[[i]]))) # Compute log determinants of all correlation matrices from cholesky factors
        cormatis[[i]] <- chol2inv(cormatis[[i]]) # Invert each correlation matrix using cholesky factor.
      }
      if (u.applyvars[1] == 0) { # If there is just a single variance.
        sigmai <- (1 / as.numeric(vars)) * Reduce(kronecker, cormatis) # Form the inverse covariance matrix.
        # Form the log determinant of sigma:
        u.nlevels.var <- u.nlevels[i.var]
        if (n.var > 1) {
          logDetvarGtotal <- 0
          for (i in i.var) {
            logDetvarGtotal <- logDetvarGtotal + logDetvarG[i] * prod(u.nlevels.var[-i])
          }
        } else {
          logDetvarGtotal <- logDetvarG
        }
        logDetvarGtotal <- prod(u.nlevels.var) * log(as.numeric(vars)) + logDetvarGtotal
      } else if (all(u.applyvars == i.var)) { # If variances are applied to all variable factors.
        sigmai <- (1 / vars) * Reduce(kronecker, cormatis) # Form the inverse covariance matrix.
        # Form the log determinant of sigma:
        u.nlevels.var <- u.nlevels[i.var]
        logDetvarGtotal <- sum(log(diag(vars)))
        for (i in i.var) {
          logDetvarGtotal <- logDetvarGtotal + logDetvarG[i] * prod(u.nlevels.var[-i])
        }
      } else { # If variances are applied to 1 or more factors, but not all factors.
        sigmai <- kronecker((1 / vars) * Reduce(kronecker, cormatis[u.applyvars]), Reduce(kronecker, cormatis[-u.applyvars])) # Form the inverse covariance matrix.
        # Form the log determinant of sigma:
        # First form the log determinant of the part with variances:
        u.nlevels.var <- u.nlevels[i.var]
        u.nlevels.applyvar <- u.nlevels[u.applyvars]
        u.nlevels.napplyvar <- u.nlevels.var[-u.applyvars]
        logDetvarGtotal <- sum(log(diag(vars)))
        for (i in u.applyvars) {
          logDetvarGtotal <- logDetvarGtotal + logDetvarG[i] * prod(u.nlevels.applyvar[-i]) # prod(numeric(0)) = 1
        }
        logDetvarGtotal <- logDetvarGtotal * prod(u.nlevels.var[-u.applyvars])
        # No we can form the full log determinant:
        for (i in i.var[-u.applyvars]) {
          logDetvarGtotal <- logDetvarGtotal + logDetvarG[i] * prod(u.nlevels.var[-i])
        }
      }
    } else { # if there are no variable factors.
      sigmai <- 1 / vars
      logDetvarGtotal <- -1 * as.numeric(log(sigmai))
    }
    if (n.nvar > 0) {
      Gi <- kronecker(sigmai, nvarGi) # This is a bit slow as well.
    } else {
      Gi <- sigmai
    }
    
    # Now add to M:
    dimGi <- prod(u.nlevels)
    M[(nX + 1):(nX + dimGi), (nX + 1):(nX + dimGi)] <- M[(nX + 1):(nX + dimGi), (nX + 1):(nX + dimGi)] + Gi # This one is also slow
    #toc(log = T)
    
    ### Cholesky decomposition of M ====
    #tic(msg = "Cholesky decomp of M")
    Mchol <- chol(M) # Upper triangular # This is good for the GPU
    #toc(log = T)
    yPy <- Mchol[nrow(Mchol), ncol(Mchol)]^2
    Mchol_XZ <- Mchol[1:(nrow(Mchol) - 1), 1:(ncol(Mchol) - 1)]
    My <- Mchol[1:(nrow(Mchol) - 1), ncol(Mchol):ncol(Mchol), drop = FALSE]
    
    ### Computing the log-likelihood ====
    #tic("Computing log-likelihood")
    # Computing logDetC:
    logDetC <- 2 * sum(log(diag(Mchol_XZ)))
    
    # Computing logDetG:
    dim.var <- prod(u.nlevels[i.var])
    dim.nvar <- prod(u.nlevels[i.nvar])
    if (n.nvar > 0) {
      logDetG <- dim.nvar * logDetvarGtotal + dim.var * logDetnvarGtotal
    } else {
      logDetG <- dim.nvar * logDetvarGtotal
    }
    
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
    Wu <- matrix(NA, n, nrow(k))
    if (n.nvar > 0) {
      umat <- matrix(u, nrow = dim.nvar, ncol = dim.var)
    } else {
      umat <- matrix(u, nrow = 1, ncol = dim.var)
    }
    umat_sigmai <- umat %*% sigmai
    derivs <- vector("list", nrow(k) - 1)
    for (i in 1:(nrow(k) - 1)) { # Assuming a single residual variance for now.
      if (isFixed[i] == 1) {
        Wu[, i] <- 0
      } else {
        if (ur.types[i] == "V") {
          if (u.applyvars[1] == 0) { # If there is only a single variance.
            if (n.var > 0) {
              derivs[[i]] <- Reduce(kronecker, cormats)
            } else {
              derivs[[i]] <- matrix(1)
            }
          } else { # If variances are applied to 1 or more factors.
            A <- matrix(0, prod(u.nlevels[u.applyvars]), prod(u.nlevels[u.applyvars]))
            A[i,] <- 0.5
            if (all(u.applyvars == i.var)) { # If variances are applied to all variable factors.
              derivs[[i]] <- (A + t(A)) * vars * Reduce(kronecker, cormats) / k[i, iter]
            } else { # If variances are applied to 1 or more factors, but not all factors.
              derivs[[i]] <- kronecker((A + t(A)) * vars * Reduce(kronecker, cormats[u.applyvars]) / k[i, iter], Reduce(kronecker, cormats[-u.applyvars]))
            }
          }
        } else if (ur.types[i] %in% c("R", "B")) {
          j <- u.forFactor[i] # Which factor does this correlation or bandwidth belong to?
          cormats1 <- cormats # Making a copy of the correlation matrices so we can replace the matrix of factor j with the indicator matrix or derivative of the GK.
          if (ur.types[i] == "R") {
            icor <- u.indexCor[[i]] # This correlation is the icor-th correlation in factor j.
            B <- matrix(0, u.nlevels[j], u.nlevels[j]) # Making the indicator matrix.
            B[lower.tri(B)][icor] <- 1 # Filling the indicator matrix.
            B <- B + t(B)
            cormats1[[j]] <- B # Replacing with the indicator matrix.
          } else if (ur.types[i] == "B") {
            cormats1[[j]] <- -1 * bases[[j]] * cormats1[[j]] # d/dh exp(-h * D) = -D * exp(-h * D)
          }
          if (u.applyvars[1] == 0) {
            derivs[[i]] <- as.numeric(vars) * Reduce(kronecker, cormats1)
          } else if (all(u.applyvars == i.var)) {
            derivs[[i]] <- vars * Reduce(kronecker, cormats1)
          } else {
            derivs[[i]] <- kronecker(vars * Reduce(kronecker, cormats1[u.applyvars]), Reduce(kronecker, cormats1[-u.applyvars]))
          }
        }
        Wu[, i] <- as.matrix(Z %*% matrix(umat_sigmai %*% derivs[[i]]))
      }
    }
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
    XWjxZWj <- Matrix::t(W) %*% Wu * (1 / k[nrow(k), iter])
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
    if (n.nvar > 0) {
      traces <- matrix(0, nrow(sigmai), ncol(sigmai))
      for (i in 1:nrow(sigmai)) {
        for (j in 1:ncol(sigmai)) {
          if (j >= i) { # Only do the upper triangle
            rstart <- (nX + 1 + (i - 1) * dim.nvar)
            rend <- rstart - 1 + dim.nvar
            cstart <- (nX + 1 + (j - 1) * dim.nvar)
            cend <- cstart - 1 + dim.nvar
            traces[i, j] <- sum(nvarGi * Ci[rstart:rend, cstart:cend])
          }
        }
      }
      # traces[lower.tri(traces)] <- traces[upper.tri(traces)] # Copy upper to lower triangle. This becomes non-symmetric?????
      if (nrow(traces) > 1) {
        traces <- traces + t(traces) - diag(diag(traces))
      }
    } else {
      traces <- Ci[(nX + 1):nrow(Ci), (nX + 1):ncol(Ci)]
    }
    scores <- matrix(0, nrow(k), 1)
    if (n.nvar > 0) {
      basisScores <- -0.5 * (dim.nvar * sigmai - sigmai %*% traces %*% sigmai - t(umat_sigmai) %*% nvarGi %*% umat_sigmai) # Last term maybe on the GPU?
    } else {
      basisScores <- -0.5 * (sigmai - sigmai %*% traces %*% sigmai - t(umat_sigmai) %*% umat_sigmai) # Last term maybe on the GPU?
    }
    for (i in 1:length(derivs)) {
      if (isFixed[i] == 1) {
        scores[i, 1] <- 0
      } else {
        scores[i, 1] <- sum(derivs[[i]] * basisScores)
      }
    }
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
      cat(sprintf("%9s %15s %20s \n", "Iteration", "LogLik", "Duration (hh:mm:ss)"))
    }
    if (constraintsUsed) {
      msg <- paste0("(", paste(c(
        if (constraintsUsed) sprintf("%d constrained", nConstrained)
      ), collapse = " and "), ")")
      cat(sprintf("%9d %15.3f %20s %5s %-s\n", iter, LL[iter], duration, "", msg))
    } else {
      cat(sprintf("%9d %15.3f %20s\n", iter, LL[iter], duration))
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
  return(list(LL = LL, k = k, AI = AI, scores = scores, b = b, u = u))
}

