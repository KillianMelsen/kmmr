# The only check/constrained we need to estimate loadings on the correlation scale is that
# for each environment, the sum of squared loadings is smaller than say (1 - 1e-6). This ensures
# that each unique variance is at least 1e-6. Ensuring the previous condition automatically ensures that the
# inner products between loadings of different environments (correlations) are within -1 to 1.

# CONTINUE WITH THIS FILE, CURRENT POSITION IS AT RUNNING ITERATION 6.

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

runIteration <- function(k, y, K, Z, X, percPenalty, isFixed, stepSize, m) {
  p <- (nrow(k) - 1) / (m + 1)
  R <- k[nrow(k), 1] * diag(nrow(y))
  logDetR <- sum(log(eigen(R)$values))
  Ri <- solve(R)
  L <- matrix(k[1:(p * m), 1], p, m)
  PSI <- diag(1 - diag(L %*% t(L)))
  vars <- outer(sqrt(k[(p * m + 1):(p * m + p), 1]), sqrt(k[(p * m + 1):(p * m + p), 1]))
  sigma <- vars * (L %*% t(L) + PSI)
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
  
  # Derivatives wrt loadings:
  dGdL <- vector("list", p * m)
  i <- 1
  for (f in 1:m) {
    for (e in 1:p) {
      if (isFixed[i]) {
        dGdL[[i]] <- 0
        Wu[, i] <- 0
        i <- i + 1
      } else {
        dGdL[[i]] <- kronecker(getLoadingDerivs(e, f, L) * vars, K)
        Wu[, i] <- dGdL[[i]] %*% Gi %*% u
        i <- i + 1
      }
    }
  }
  
  # Derivatives wrt variances:
  dGdV <- vector("list", p)
  for (e in 1:p) {
    if (isFixed[i]) {
      dGdV[[i - p * m]] <- 0
      Wu[, i - p * m] <- 0
      i <- i + 1
    } else {
      deriv <- matrix(0, nrow(L), nrow(L))
      deriv[e,] <- 0.5
      deriv <- (deriv + t(deriv)) * sigma / k[i, 1]
      deriv <- kronecker(deriv, K)
      dGdV[[i - p * m]] <- deriv
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
  derivs <- c(dGdL, dGdV)
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
                       random = ~ sfa(Env, 3):vm(Gen, Ki),
                       residual = ~ idv(units),
                       data = d,
                       maxit = 20)

# Manual:
p <- length(levels(d$Env))
init <- as.numeric(fit.asr$trace[5:47, 1, drop=F])
# init <- c(getInitFA(d = d, m = 2, response = "GY", factor = "Env"), var(d$GY) / 10 / 2)
isFixed <- abs(init) < 1e-16
k <- matrix(NA, length(init), 20)
rownames(k) <- c(sprintf("L%d1", 1:14),
                 sprintf("L%d2", 1:14),
                 sprintf("VAR%d", 1:14),
                 "rvar")
names(isFixed) <- rownames(k)
LL <- numeric(20)
k[, 1] <- init
X <- model.matrix(GY ~ -1 + Env, d)
Z <- diag(nrow(d))
W <- cbind(X, Z)
y <- matrix(d$GY)

# Iteration 1 ====
iter <- 1
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25, m = m)
LL[iter] <- result$LL
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 2 ====
iter <- 2
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25, m = m)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 3 ====
iter <- 3
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25, m = m)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 4 ====
iter <- 4
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25, m = m)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 5 ====
iter <- 5
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.5, m = m)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# # We need to constrain the loadings in row 3 to ensure PSI3 remains positive:
con <- c(3, 17) # L31 and L32
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# We now also need to shrink row 13, and constrain L42 because it flips sign and increases in magnitude by more than 2x:
con <- c(13, 27, 18) # L131, L132, and L42
con.total <- sort(unique(c(con.total, con)))
s <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2)) # Scaling factor to ensure sum of squared loadings remains lower than 1-1e-6
delta[con,] <- c(s * kUpdate[con[1:2]] - k[con[1:2], iter], # Sum of squared loadings must be lower than 1-1e-6
                 -2 * k[con[3], iter] - k[con[3], iter]) 
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 6 ====
iter <- 6
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.5, m = m)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# We have to constrain the 3rd and 13th row of L again:
con <- c(3, 17, 13, 27)
con.total <- sort(unique(c(con.total, con)))
sL3 <- sqrt((1 - 1e-6) / sum(kUpdate[con[1:2]]^2)) # Scaling factor to ensure sum of squared loadings remains lower than 1-1e-6
sL13 <- sqrt((1 - 1e-6) / sum(kUpdate[con[3:4]]^2)) # Scaling factor to ensure sum of squared loadings remains lower than 1-1e-6
delta[con,] <- c(sL3 * kUpdate[con[1:2]] - k[con[1:2], iter],
                 sL13 * kUpdate[con[3:4]] - k[con[3:4], iter]) # Sum of squared loadings must be lower than 1-1e-6
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 7 ====
iter <- 7
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.5, m = m)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
# Same issue again for the third time with row 13, shrink one more time and fix at those values:
con <- c(13, 27)
con.total <- sort(unique(c(con.total, con)))
sL13 <- sqrt((1 - 1e-6) / sum(kUpdate[con]^2)) # Scaling factor to ensure sum of squared loadings remains lower than 1-1e-6
delta[con,] <- c(sL13 * kUpdate[con] - k[con, iter]) # Sum of squared loadings must be lower than 1-1e-6
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
isFixed[con] <- TRUE
k[, iter + 1] <- kUpdate

# Iteration 8 ====
iter <- 8
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1, m = m)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 9 ====
iter <- 9
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1, m = m)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate

# Iteration 10 ====
iter <- 10
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1, m = m)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:43, , drop = F]) # All should be between 1/15 (0.067) and 15
L <- matrix(kUpdate[1:28, 1], 14, 2)
matrix(rowSums(L^2)) # Should all be lower than 1 - 1e-6 (0.999999)
k[, iter + 1] <- kUpdate
# Converged!

L.asr <- matrix(summary(fit.asr)$varcomp[1:28, 1], 14, 2)
L.kmmr <- matrix(k[1:28, iter + 1], 14, 2)

vars.asr <- as.numeric(summary(fit.asr)$varcomp[29:43, 1])
vars.kmmr <- as.numeric(k[29:43, iter + 1])

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