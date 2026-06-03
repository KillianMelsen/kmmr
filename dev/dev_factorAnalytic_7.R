# The only check/constrained we need to estimate loadings on the correlation scale is that
# for each environment, the sum of squared loadings is smaller than say (1 - 1e-6). This ensures
# that each unique variance is at least 1e-6. Ensuring the previous condition automatically ensures that the
# inner products between loadings of different environments (correlations) are within -1 to 1.

# Custom functions:
rotateInitL <- function(L) {
  L.rotated <- L
  m <- ncol(L.rotated)
  p <- nrow(L.rotated)
  coords <- which(upper.tri(matrix(0, p, m)), arr.ind = TRUE)
  givens <- vector("list", nrow(coords))
  for (i in 1:nrow(coords)) {
    row <- coords[i, "row"]
    col <- coords[i, "col"]
    givens[[i]] <- diag(m)
    h <- sqrt(sum(c(L.rotated[row, row], L.rotated[row, col])^2))
    c <- L.rotated[row, row] / h
    s <- L.rotated[row, col] / h
    diag(givens[[i]])[c(row, col)] <- c
    givens[[i]][row, col] <- -s
    givens[[i]][col, row] <- s
    L.rotated <- L.rotated %*% givens[[i]]
  }
  rotmat <- Reduce("%*%", givens)
  return(list(L.rotated = L.rotated, rotmat = rotmat))
}

getInitFA <- function(d, m, response, factor) {
  covmat.phen <- cov(unstack(d, as.formula(paste0(response, " ~ ", factor))), use = "complete")
  # cormat.phen <- cov2cor(covmat.phen)
  # covmat.init <- matrix(0, nrow(covmat.phen), ncol(covmat.phen))
  # rownames(covmat.init) <- colnames(covmat.init) <- colnames(covmat.phen)
  # covmat.init[lower.tri(covmat.init)] <- cov2cor(covmat.phen)[lower.tri(covmat.phen)]
  # covmat.init <- covmat.init + t(covmat.init)
  # diag(covmat.init) <- diag(covmat.phen) / 10
  vars <- diag(covmat.phen) / 10
  fa.fit <- factanal(factors = m, covmat = covmat.phen)
  L <- fa.fit$loadings[,]
  # L.rotated <- L
  L.rotated <- rotateInitL(L)$L.rotated
  # L.rotated <- diag(sqrt(vars)) %*% L.rotated
  PSI <- vars
  # PSI <- vars * fa.fit$uniquenesses
  init <- c(as.numeric(L.rotated),
            as.numeric(PSI))
  return(init)
}

getLoadingDerivs <- function(rowL, colL, L) {
  deriv <- matrix(0, nrow(L), nrow(L))
  deriv[rowL,] <- deriv[,rowL] <- L[, colL]
  # deriv[rowL, rowL] <- 2 * L[rowL, colL]
  deriv[rowL, rowL] <- 0
  return(deriv)
}

runIteration <- function(k, y, K, Z, X, percPenalty, isFixed, stepSize, m, LK) {
  p <- (nrow(k) - 1) / (m + 2)
  R <- k[nrow(k), 1] * diag(nrow(y))
  logDetR <- sum(log(eigen(R)$values))
  Ri <- solve(R)
  L <- matrix(k[(p + 1):(p + p * m), 1], p, m)
  PSI <- diag(1 - diag(L %*% t(L)))
  varsFA <- outer(sqrt(k[(p + p * m + 1):(p + p * m + p), 1]), sqrt(k[(p + p * m + 1):(p + p * m + p), 1]))
  varsLK <- outer(sqrt(k[1:p, 1]), sqrt(k[1:p, 1]))
  sigmaLK <- varsLK * LK
  sigmaFA <- varsFA * (L %*% t(L) + PSI)
  sigma <- sigmaLK + sigmaFA
  G <- kronecker(sigma, K)
  logDetG <- sum(log(eigen(G)$values))
  Gi <- solve(G)
  V <- Z %*% G %*% t(Z) + R
  Vi <- solve(V)
  P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi
  C <- rbind(cbind(t(X) %*% Ri %*% X, t(X) %*% Ri %*% Z),
             cbind(t(Z) %*% Ri %*% X, t(Z) %*% Ri %*% Z + Gi))
  Ci <- solve(C)
  logDetC <- sum(log(eigen(C)$values))
  LL <- -0.5 * (logDetC + logDetR + logDetG + t(y) %*% P %*% y)
  RHS <- rbind(t(X) %*% Ri %*% y,
               t(Z) %*% Ri %*% y)
  bu <- Ci %*% RHS
  b <- bu[1:p, 1, drop = F]
  u <- bu[(p + 1):nrow(bu), 1, drop = F]
  resids <- y - W %*% bu
  
  # Matrix of working variates:
  Wu <- matrix(0, nrow(d), nrow(k))
  
  # Derivatives wrt to LK variances:
  i <- 1
  dGdVLK <- vector("list", p)
  for (e in 1:p) {
    if (isFixed[i]) {
      dGdVLK[[i]] <- 0
      Wu[, i] <- 0
      i <- i + 1
    } else {
      deriv <- matrix(0, nrow(LK), ncol(LK))
      deriv[e,] <- 0.5
      deriv <- (deriv + t(deriv)) * sigmaLK / k[i, 1]
      deriv <- kronecker(deriv, K)
      dGdVLK[[i]] <- deriv
      Wu[, i] <- deriv %*% Gi %*% u
      i <- i + 1
    }
  }
  
  # Derivatives wrt loadings:
  dGdL <- vector("list", p * m)
  for (f in 1:m) {
    for (e in 1:p) {
      if (isFixed[i]) {
        dGdL[[i - p]] <- 0
        Wu[, i] <- 0
        i <- i + 1
      } else {
        dGdL[[i - p]] <- kronecker(getLoadingDerivs(e, f, L) * varsFA, K)
        Wu[, i] <- dGdL[[i - p]] %*% Gi %*% u
        i <- i + 1
      }
    }
  }
  
  # Derivatives wrt FA variances:
  dGdVFA <- vector("list", p)
  for (e in 1:p) {
    if (isFixed[i]) {
      dGdVFA[[i - p - p * m]] <- 0
      Wu[, i] <- 0
      i <- i + 1
    } else {
      deriv <- matrix(0, nrow(L), nrow(L))
      deriv[e,] <- 0.5
      deriv <- (deriv + t(deriv)) * sigmaFA / k[i, 1]
      deriv <- kronecker(deriv, K)
      dGdVFA[[i - p - p * m]] <- deriv
      Wu[, i] <- deriv %*% Gi %*% u
      i <- i + 1
    }
  }
  
  # Residual working variate:
  Wu[, i] <- Ri %*% resids
  
  AI <- t(Wu) %*% P %*% Wu / 2
  rownames(AI) <- colnames(AI) <- rownames(k)
  diag(AI)[1:(p * m)] <- diag(AI)[1:(p * m)] + percPenalty * diag(AI)[1:(p * m)]
  AIi <- matrix(0, nrow(AI), ncol(AI))
  AIi[!isFixed, !isFixed] <- solve(AI[!isFixed, !isFixed])
  
  # Scores:
  derivs <- c(dGdVLK, dGdL, dGdVFA)
  scores <- matrix(0, nrow(k), 1)
  rownames(scores) <- rownames(k)
  for (i in 1:length(derivs)) {
    if (isFixed[i]) {
      scores[i, 1] <- 0
    } else {
      scores[i, 1] <- -0.5 * (sum(diag(derivs[[i]] %*% Gi)) - sum(diag((Gi %*% derivs[[i]] %*% Gi) %*% Ci[(p + 1):nrow(Ci), (p + 1):nrow(Ci)])) - t(y) %*% P %*% Wu[, i]) 
    }
  }
  
  scores[nrow(scores), 1] <- -0.5 * (sum(diag(Ri)) - sum(diag((t(W) %*% Ri %*% Ri %*% W) %*% Ci)) - t(y) %*% P %*% Wu[, ncol(Wu)])
  
  # Computing delta:
  delta <- stepSize * AIi %*% scores
  rownames(delta) <- rownames(k)
  return(list(LL = LL, delta = delta, AI = AI, scores = scores, Wu = Wu, u = u))
}

library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = TRUE, score = TRUE, rotate.fa = FALSE)

d <- readRDS("~/repos/kmmr/dev/weights_data.rds")
EC <- d$EC
ED <- d$ED
d <- d$ydata
K <- readRDS("~/repos/kmmr/dev/weights_K.rds")
Ki <- solve(K)
attr(Ki, "INVERSE") <- TRUE

d <- droplevels(d[d$Man == "LN",])
m <- 2
fit.asr <- asreml(fixed = GY ~ -1 + Env,
                       random = ~ corgh(Env):vm(Gen, Ki) + sfa(Env, m):vm(Gen, Ki),
                       residual = ~ idv(units),
                       data = d,
                       maxit = 20,
                  start.values = T)
vc <- fit.asr$vparameters.table
vc$Value[1:91] <- EC[upper.tri(EC)]
vc$Constraint[1:91] <- "F"
fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ corgh(Env):vm(Gen, Ki) + sfa(Env, m):vm(Gen, Ki),
                  residual = ~ idv(units),
                  data = d,
                  maxit = 20,
                  G.param = vc)

# Manual:
p <- length(levels(d$Env))
init <- as.numeric(fit.asr$trace[96:152, 1, drop=F])
# init <- c(getInitFA(d = d, m = 2, response = "GY", factor = "Env"), var(d$GY) / 10 / 2)
# init <- c(init[29:42], init)
isFixed <- abs(init) < 1e-16
k <- matrix(NA, length(init), 20)
rownames(k) <- c(sprintf("LK-VAR%d", 1:14),
                 sprintf("L%d1", 1:14),
                 sprintf("L%d2", 1:14),
                 sprintf("FA-VAR%d", 1:14),
                 "rvar")
names(isFixed) <- rownames(k)
LL <- numeric(20)
k[, 1] <- init
k["L12", 1] <- 0
X <- model.matrix(GY ~ -1 + Env, d)
Z <- diag(nrow(d))
W <- cbind(X, Z)
y <- matrix(d$GY)

# Iteration 1 ====
iter <- 1
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.10, isFixed = isFixed, stepSize = 0.10, m = m, LK = EC)
LL[iter] <- result$LL

# Some checks:
all.equal(LL[1], fit.asr$loglik)
all.equal(as.numeric(fit.asr$score)[92:148], as.numeric(result$scores))
all.equal(as.matrix(fit.asr$wrkv)[, c(92:119, 121:148)], result$Wu[, c(1:28, 30:57)]) # Exclude the working variate for L12 because that one is nonzero for asreml even though L12 is fixed...

delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# We need to constrain the loadings for environment 8:
con <- c(22, 36) # L81 and L82
con.total <- sort(unique(c(con.total, con)))
s <- sqrt((1 - 1e-6) / sum(kUpdate[con]^2)) # Scaling factor to ensure sum of squared loadings remains lower than 1-1e-6
delta[con,] <- c(s * kUpdate[con] - k[con, iter]) # Sum of squared loadings must be lower than 1-1e-6
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# We need to constrain the loadings for environment 3, 4, 5, 6, 9, 13, and LK-VAR6, LK-VAR2:
con <- c(17, 31,
         18, 32,
         19, 33,
         20, 34,
         23, 37,
         27, 41,
         2,
         6)
con.total <- sort(unique(c(con.total, con)))
s3 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2))
s4 <- sqrt((1 - 1e-6) / sum(kUpdate[con[3:4]]^2))
s5 <- sqrt((1 - 1e-6) / sum(kUpdate[con[5:6]]^2))
s6 <- sqrt((1 - 1e-6) / sum(kUpdate[con[7:8]]^2))
s9 <- sqrt((1 - 1e-6) / sum(kUpdate[con[9:10]]^2))
s13 <- sqrt((1 - 1e-6) / sum(kUpdate[con[11:12]]^2))
delta[con,] <- c(s3 * kUpdate[con[1:2]] - k[con[1:2], iter],
                 s4 * kUpdate[con[3:4]] - k[con[3:4], iter],
                 s5 * kUpdate[con[5:6]] - k[con[5:6], iter],
                 s6 * kUpdate[con[7:8]] - k[con[7:8], iter],
                 s9 * kUpdate[con[9:10]] - k[con[9:10], iter],
                 s13 * kUpdate[con[11:12]] - k[con[11:12], iter],
                 1 / 15 * k[2, iter] - k[2, iter],
                 1 / 15 * k[6, iter] - k[6, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Now constrain LK-VAR12:
con <- c(12)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1 / 15 * k[12, iter] - k[12, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 2 ====
iter <- 2
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Need to constrain loadings for 5, 8, and 13:
con <- c(19, 33,
         22, 36,
         27, 41)
con.total <- sort(unique(c(con.total, con)))
s5 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2))
s8 <- sqrt((1 - 1e-6) / sum(kUpdate[con[3:4]]^2))
s13 <- sqrt((1 - 1e-6) / sum(kUpdate[con[5:6]]^2))
delta[con,] <- c(s5 * kUpdate[con[1:2]] - k[con[1:2], iter],
                 s8 * kUpdate[con[3:4]] - k[con[3:4], iter],
                 s13 * kUpdate[con[5:6]] - k[con[5:6], iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Now constrain loadings for E6 and LK-VAR5:
con <- c(20, 34,
         5)
con.total <- sort(unique(c(con.total, con)))
s6 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2))
delta[con,] <- c(s6 * kUpdate[con[1:2]] - k[con[1:2], iter],
                 1 / 15 * k[5, iter] - k[5, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 3 ====
iter <- 3
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Need to constrain loading rows 8 and 13:
con <- c(22, 36,
         27, 41)
con.total <- sort(unique(c(con.total, con)))
s8 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2))
s13 <- sqrt((1 - 1e-6) / sum(kUpdate[con[3:4]]^2))
delta[con,] <- c(s8 * kUpdate[con[1:2]] - k[con[1:2], iter],
                 s13 * kUpdate[con[3:4]] - k[con[3:4], iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Now constrain LK-VAR2, LK-VAR3 and LK-VAR6:
con <- c(2, 3, 6)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1 / 15 * k[2, iter] - k[2, iter],
                 1 / 15 * k[3, iter] - k[3, iter],
                 1 / 15 * k[6, iter] - k[6, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 4 ====
iter <- 4
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Need to constrain LK-VAR2. Also constrain and fix L8 and L13:
con <- c(22, 36,
         27, 41,
         2)
con.total <- sort(unique(c(con.total, con)))
s8 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2))
s13 <- sqrt((1 - 1e-6) / sum(kUpdate[con[3:4]]^2))
delta[con,] <- c(s8 * kUpdate[con[1:2]] - k[con[1:2], iter],
                 s13 * kUpdate[con[3:4]] - k[con[3:4], iter],
                 1 / 15 * k[2, iter] - k[2, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Now constrain L122, L5 and L6:
con <- c(19, 33,
         20, 34,
         40)
con.total <- sort(unique(c(con.total, con)))
s5 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2))
s6 <- sqrt((1 - 1e-6) / sum(kUpdate[con[3:4]]^2))
delta[con,] <- c(s5 * kUpdate[con[1:2]] - k[con[1:2], iter],
                 s6 * kUpdate[con[3:4]] - k[con[3:4], iter],
                 -2 * k[40, iter] - k[40, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
isFixed[c("L81", "L82", "L131", "L132")] <- TRUE
k[, iter + 1] <- kUpdate

# Iteration 5 ====
iter <- 5
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Constrain and fix LK-VAR2:
con <- c(2)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1 / 15 * k[2, iter] - k[2, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Now constrain LK-VAR12 and L122:
con <- c(12, 40)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1 / 15 * k[12, iter] - k[12, iter],
                 -2 * k[40, iter] - k[40, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
isFixed["LK-VAR2"] <- TRUE
k[, iter + 1] <- kUpdate

# Iteration 6 ====
iter <- 6
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 7 ====
iter <- 7
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.5, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Constrain and fix L6:
con <- c(20, 34)
con.total <- sort(unique(c(con.total, con)))
s6 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2))
delta[con,] <- c(s6 * kUpdate[con[1:2]] - k[con[1:2], iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
isFixed[c("L61", "L62")] <- TRUE
k[, iter + 1] <- kUpdate

# Iteration 8 ====
iter <- 8
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.5, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Need to constrain LK-VAR12
con <- c(12)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1 / 15 * k[12, iter] - k[12, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Now constrain L2:
con <- c(16, 30)
con.total <- sort(unique(c(con.total, con)))
s2 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2))
delta[con,] <- c(s2 * kUpdate[con[1:2]] - k[con[1:2], iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 9 ====
iter <- 9
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.5, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Constrain and fix LK-VAR12:
con <- c(12)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1e-6 - k[con, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Now constrain LK-VAR14 and L2 as well:
con <- c(16, 30,
         14)
con.total <- sort(unique(c(con.total, con)))
s2 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2))
delta[con,] <- c(s2 * kUpdate[con[1:2]] - k[con[1:2], iter],
                 1 / 15 * k[14, iter] - k[14, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
isFixed["LK-VAR12"] <- TRUE
k[, iter + 1] <- kUpdate

# Iteration 10 ====
iter <- 10
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.5, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Constrain LK-VAR14 again
con <- c(14)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1 / 15 * k[14, iter] - k[14, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 11 ====
iter <- 11
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.5, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Constrain LK-VAR14 for the 3rd time:
con <- c(14)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1 / 15 * k[14, iter] - k[14, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 12 ====
iter <- 12
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.5, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Constrain and fix LK-VAR14:
con <- c(14)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1e-6 - k[14, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
isFixed["LK-VAR14"] <- TRUE
k[, iter + 1] <- kUpdate

# Iteration 13 ====
iter <- 13
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 14 ====
iter <- 14
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Constrain L112:
con <- c(39)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(-2 * k[39, iter] - k[39, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 15 ====
iter <- 15
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 16 ====
iter <- 16
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.02, isFixed = isFixed, stepSize = 1, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 17 ====
iter <- 17
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.00, isFixed = isFixed, stepSize = 1, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Constrain L4:
con <- c(18, 32)
con.total <- sort(unique(c(con.total, con)))
s4 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2))
delta[con,] <- c(s4 * kUpdate[con[1:2]] - k[con[1:2], iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 18 ====
iter <- 18
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.00, isFixed = isFixed, stepSize = 1, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Constrain L4 for the 3rd time:
con <- c(18, 32)
con.total <- sort(unique(c(con.total, con)))
s4 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2))
delta[con,] <- c(s4 * kUpdate[con[1:2]] - k[con[1:2], iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 19 ====
iter <- 19
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.00, isFixed = isFixed, stepSize = 1, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Constrain LK-VAR11 and fix L4:
con <- c(18, 32,
         11)
con.total <- sort(unique(c(con.total, con)))
s4 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2))
delta[con,] <- c(s4 * kUpdate[con[1:2]] - k[con[1:2], iter],
                 1 / 15 * k[11, iter] - k[11, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
isFixed[c("L41", "L42")] <- TRUE
k[, iter + 1] <- kUpdate

# Iteration 20 ====
iter <- 20
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.00, isFixed = isFixed, stepSize = 1, m = m, LK = EC)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Constrain and fix LK-VAR11:
con <- c(11)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1e-6 - k[11, iter])
# Recompute delta for the unconstrained parameters:
noncon <- setdiff(1:nrow(k), con.total)
uuisFixed <- isFixed[noncon]
AIuu <- AI[noncon, noncon, drop = FALSE]
AIff <- AI[con.total, con.total, drop = FALSE]
AIuf <- AI[noncon, con.total, drop = FALSE]
scoresuu <- scores[noncon, 1, drop = FALSE]
scoresff <- scores[con.total, 1, drop = FALSE]
deltaff <- delta[con.total, 1, drop = FALSE]
AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
uunotFixed <- which(uuisFixed == 0)
AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
delta[noncon, 1] <- deltauu
# Recompute kUpdate:
kUpdate <- k[, iter, drop = F] + delta
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[15:42, iter]) * sign(kUpdate[15:42,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[15:42, , drop = F]
(meanMags <- mean(mags[!isFixed[15:42]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[c(1:14, 43:57), , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[15:42, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Converged! sort of...
vars.LK.asr <- summary(fit.asr)$varcomp[92:105, 1]
vars.LK.kmmr <- kUpdate[1:14]

L.asr <- matrix(summary(fit.asr)$varcomp[106:133, 1], 14, 2)
L.kmmr <- matrix(kUpdate[15:42], 14, 2)

vars.FA.asr <- as.numeric(summary(fit.asr)$varcomp[134:147, 1])
vars.FA.kmmr <- kUpdate[43:56]

cor.asr <- L.asr %*% t(L.asr); diag(cor.asr) <- 1
cor.kmmr <- L.kmmr %*% t(L.kmmr); diag(cor.kmmr) <- 1

rownames(L.kmmr) <- levels(d$Env)
colnames(L.kmmr) <- paste0("F", 1:m)
L <- varimax(L.kmmr)$loadings

library(ggplot2)
# 2. Convert to a data frame for ggplot
plot_data <- data.frame(Feature = rownames(L),
                        F1 = L[, "F1"],
                        F2 = L[, "F2"],
                        F3 = L[, "F3"])

# Create points for a perfect unit circle
angle <- seq(0, 2 * pi, length.out = 100)
circle_data <- data.frame(x = cos(angle),
                          y = sin(angle))

p2 <- ggplot(data = plot_data) +
  # 1. Add a unit circle for reference
  geom_path(data = circle_data, aes(x = x, y = y), color = "red", linetype = "dashed") +
  # 2. Add baseline axes passing through the origin
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
  # 3. Plot the feature loadings as vectors (arrows)
  geom_segment(aes(x = 0, y = 0, xend = F1, yend = F2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black", lwd = 0.8) +
  # 4. Add text labels for each feature
  geom_text(aes(x = F1, y = F2, label = Feature), 
            hjust = -0.2, vjust = 0.5, size = 3.5, check_overlap = TRUE) +
  # 5. Force a 1:1 aspect ratio so the circle doesn't look like an ellipse
  coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
  # Labels and Styling
  labs(x = "F1",
       y = "F2") +
  theme_minimal()

# p3 <- ggplot(data = plot_data) +
#   # 1. Add a unit circle for reference
#   geom_path(data = circle_data, aes(x = x, y = y), color = "red", linetype = "dashed") +
#   # 2. Add baseline axes passing through the origin
#   geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
#   geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
#   # 3. Plot the feature loadings as vectors (arrows)
#   geom_segment(aes(x = 0, y = 0, xend = F1, yend = F3),
#                arrow = arrow(length = unit(0.2, "cm")), 
#                color = "black", lwd = 0.8) +
#   # 4. Add text labels for each feature
#   geom_text(aes(x = F1, y = F3, label = Feature), 
#             hjust = -0.2, vjust = 0.5, size = 3.5, check_overlap = TRUE) +
#   # 5. Force a 1:1 aspect ratio so the circle doesn't look like an ellipse
#   coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
#   # Labels and Styling
#   labs(x = "F1",
#        y = "F3") +
#   theme_minimal()

# p4 <- ggplot(data = plot_data) +
#   # 1. Add a unit circle for reference
#   geom_path(data = circle_data, aes(x = x, y = y), color = "red", linetype = "dashed") +
#   # 2. Add baseline axes passing through the origin
#   geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
#   geom_vline(xintercept = 0, linetype = "dotted", color = "red") +
#   # 3. Plot the feature loadings as vectors (arrows)
#   geom_segment(aes(x = 0, y = 0, xend = F2, yend = F3),
#                arrow = arrow(length = unit(0.2, "cm")), 
#                color = "black", lwd = 0.8) +
#   # 4. Add text labels for each feature
#   geom_text(aes(x = F2, y = F3, label = Feature), 
#             hjust = -0.2, vjust = 0.5, size = 3.5, check_overlap = TRUE) +
#   # 5. Force a 1:1 aspect ratio so the circle doesn't look like an ellipse
#   coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
#   # Labels and Styling
#   labs(x = "F2",
#        y = "F3") +
#   theme_minimal()

p1 <- ggplotify::as.ggplot(expression(
  heatmap(cor.kmmr)))

p1 | p2# | p3 | p4

library(plotly)