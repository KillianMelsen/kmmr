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
  covmat.init <- matrix(0, nrow(covmat.phen), ncol(covmat.phen))
  rownames(covmat.init) <- colnames(covmat.init) <- colnames(covmat.phen)
  covmat.init[lower.tri(covmat.init)] <- cov2cor(covmat.phen)[lower.tri(covmat.phen)]
  covmat.init <- covmat.init + t(covmat.init)
  diag(covmat.init) <- diag(covmat.phen) / 10
  vars <- diag(covmat.init)
  fa.fit <- factanal(factors = m, covmat = covmat.init)
  L <- fa.fit$loadings[,]
  L.rotated <- rotateInitL(L)$L.rotated
  L.rotated <- diag(sqrt(vars)) %*% L.rotated
  PSI <- vars * fa.fit$uniquenesses
  init <- c(as.numeric(L.rotated),
            as.numeric(PSI))
  return(init)
}

getLoadingDerivs <- function(rowL, colL, L) {
  deriv <- matrix(0, nrow(L), nrow(L))
  deriv[rowL,] <- deriv[,rowL] <- L[, colL]
  deriv[rowL, rowL] <- 2 * L[rowL, colL]
  return(deriv)
}

runIteration <- function(k, y, K, Z, X, percPenalty, isFixed, stepSize) {
  R <- k[nrow(k), 1] * diag(nrow(y))
  logDetR <- sum(log(eigen(R)$values))
  Ri <- solve(R)
  L <- matrix(k[1:10, 1], 5, 2)
  PSI <- diag(k[11:15, 1])
  G <- kronecker(L %*% t(L) + PSI, K)
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
  b <- bu[1:5, 1, drop = F]
  u <- bu[6:1015, 1, drop = F]
  e <- y - W %*% bu
  
  {
    dGdL11 <- kronecker(getLoadingDerivs(1, 1, L), K)
    dGdL21 <- kronecker(getLoadingDerivs(2, 1, L), K)
    dGdL31 <- kronecker(getLoadingDerivs(3, 1, L), K)
    dGdL41 <- kronecker(getLoadingDerivs(4, 1, L), K)
    dGdL51 <- kronecker(getLoadingDerivs(5, 1, L), K)
    
    dGdL12 <- kronecker(getLoadingDerivs(1, 2, L), K)
    dGdL22 <- kronecker(getLoadingDerivs(2, 2, L), K)
    dGdL32 <- kronecker(getLoadingDerivs(3, 2, L), K)
    dGdL42 <- kronecker(getLoadingDerivs(4, 2, L), K)
    dGdL52 <- kronecker(getLoadingDerivs(5, 2, L), K)
    
    dGdPSI1 <- matrix(0, nrow(L), nrow(L)); dGdPSI1[1, 1] <- 1; dGdPSI1 <- kronecker(dGdPSI1, K)
    dGdPSI2 <- matrix(0, nrow(L), nrow(L)); dGdPSI2[2, 2] <- 1; dGdPSI2 <- kronecker(dGdPSI2, K)
    dGdPSI3 <- matrix(0, nrow(L), nrow(L)); dGdPSI3[3, 3] <- 1; dGdPSI3 <- kronecker(dGdPSI3, K)
    dGdPSI4 <- matrix(0, nrow(L), nrow(L)); dGdPSI4[4, 4] <- 1; dGdPSI4 <- kronecker(dGdPSI4, K)
    dGdPSI5 <- matrix(0, nrow(L), nrow(L)); dGdPSI5[5, 5] <- 1; dGdPSI5 <- kronecker(dGdPSI5, K)
    
    Wu <- matrix(0, nrow(d), nrow(k))
    Wu[, 1] <- dGdL11 %*% Gi %*% u
    Wu[, 2] <- dGdL21 %*% Gi %*% u
    Wu[, 3] <- dGdL31 %*% Gi %*% u
    Wu[, 4] <- dGdL41 %*% Gi %*% u
    Wu[, 5] <- dGdL51 %*% Gi %*% u
    
    Wu[, 6] <-  0
    Wu[, 7] <-  dGdL22 %*% Gi %*% u
    Wu[, 8] <-  dGdL32 %*% Gi %*% u
    Wu[, 9] <-  dGdL42 %*% Gi %*% u
    Wu[, 10] <- dGdL52 %*% Gi %*% u
    
    Wu[, 11] <- dGdPSI1 %*% Gi %*% u
    Wu[, 12] <- dGdPSI2 %*% Gi %*% u
    Wu[, 13] <- dGdPSI3 %*% Gi %*% u
    Wu[, 14] <- dGdPSI4 %*% Gi %*% u
    Wu[, 15] <- dGdPSI5 %*% Gi %*% u
    
    Wu[, 16] <- Ri %*% e
  }
  
  AI <- t(Wu) %*% P %*% Wu / 2
  rownames(AI) <- colnames(AI) <- rownames(k)
  diag(AI)[1:10] <- diag(AI)[1:10] + percPenalty * diag(AI)[1:10]
  AIi <- matrix(0, nrow(AI), ncol(AI))
  AIi[!isFixed, !isFixed] <- solve(AI[!isFixed, !isFixed])
  
  scores <- matrix(0, nrow(k), 1)
  rownames(scores) <- rownames(k)
  {
    scores[1, 1] <-  -0.5 * (sum(diag(dGdL11 %*% Gi)) - sum(diag((Gi %*% dGdL11 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 1])
    scores[2, 1] <-  -0.5 * (sum(diag(dGdL21 %*% Gi)) - sum(diag((Gi %*% dGdL21 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 2])
    scores[3, 1] <-  -0.5 * (sum(diag(dGdL31 %*% Gi)) - sum(diag((Gi %*% dGdL31 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 3])
    scores[4, 1] <-  -0.5 * (sum(diag(dGdL41 %*% Gi)) - sum(diag((Gi %*% dGdL41 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 4])
    scores[5, 1] <-  -0.5 * (sum(diag(dGdL51 %*% Gi)) - sum(diag((Gi %*% dGdL51 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 5])
    
    scores[6, 1] <-  0
    scores[7, 1] <-  -0.5 * (sum(diag(dGdL22 %*% Gi)) - sum(diag((Gi %*% dGdL22 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 7])
    scores[8, 1] <-  -0.5 * (sum(diag(dGdL32 %*% Gi)) - sum(diag((Gi %*% dGdL32 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 8])
    scores[9, 1] <-  -0.5 * (sum(diag(dGdL42 %*% Gi)) - sum(diag((Gi %*% dGdL42 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 9])
    scores[10, 1] <- -0.5 * (sum(diag(dGdL52 %*% Gi)) - sum(diag((Gi %*% dGdL52 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 10])
    
    scores[11, 1] <- -0.5 * (sum(diag(dGdPSI1 %*% Gi)) - sum(diag((Gi %*% dGdPSI1 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 11])
    scores[12, 1] <- -0.5 * (sum(diag(dGdPSI2 %*% Gi)) - sum(diag((Gi %*% dGdPSI2 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 12])
    scores[13, 1] <- -0.5 * (sum(diag(dGdPSI3 %*% Gi)) - sum(diag((Gi %*% dGdPSI3 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 13])
    scores[14, 1] <- -0.5 * (sum(diag(dGdPSI4 %*% Gi)) - sum(diag((Gi %*% dGdPSI4 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 14])
    scores[15, 1] <- -0.5 * (sum(diag(dGdPSI5 %*% Gi)) - sum(diag((Gi %*% dGdPSI5 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 15])
    
    scores[16, 1] <- -0.5 * (sum(diag(Ri)) - sum(diag((t(W) %*% Ri %*% Ri %*% W) %*% Ci)) - t(y) %*% P %*% Wu[, 16])
  }
  
  # Investigating differences between asreml AI matrix and my AI matrix:
  # AIi.asr <- as.matrix(fit.asr$ai)[-17, -17]
  # rownames(AIi.asr) <- colnames(AIi.asr) <- colnames(AI)
  # AI.asr <- matrix(0, 16, 16)
  # rownames(AI.asr) <- colnames(AI.asr) <- colnames(AI)
  # AI.asr[-6, -6] <- solve(AIi.asr[-6, -6])
  # Wu.asr <- as.matrix(fit.asr$wrkv[, -17])
  # all.equal((t(Wu.asr) %*% P %*% Wu.asr / 2)[-6,-6], AI[-6,-6], check.attributes = FALSE)
  # all.equal(as.numeric(fit.asr$score)[-17], as.numeric(scores))
  
  # Computing delta:
  # delta <- fit.asr$trace[4,1]^2 * AIi %*% scores
  delta <- stepSize * AIi %*% scores
  rownames(delta) <- rownames(k)
  return(list(LL = LL, delta = delta, AI = AI, scores = scores))
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

d <- droplevels(d[d$Man == "LN" & d$Env %in% c("HAN_2015", "KAL_2015", "KIE_2015", "QLB_2015", "RHH_2015"),])

fit.asr <- asreml(fixed = GY ~ -1 + Env,
                       random = ~ facv(Env, 2):vm(Gen, Ki),
                       residual = ~ idv(units),
                       data = d,
                       maxit = 1)

# Manual:
p <- length(levels(d$Env))
# init <- as.numeric(fit.asr$trace[5:20,1,drop=F])
init <- c(getInitFA(d = d, m = 2, response = "GY", factor = "Env"), var(d$GY) / 10 / 2)
isFixed <- abs(init) < 1e-16
k <- matrix(NA, 16, 20)
rownames(k) <- c(sprintf("L%d1", 1:5),
                 sprintf("L%d2", 1:5),
                 sprintf("PSI%d", 1:5),
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
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.00, isFixed = isFixed, stepSize = 0.25)
LL[iter] <- result$LL
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
# We need to constrain L22:
con <- 7 # L22 is the 7th parameter
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(-2 * k[7, iter] - k[7, iter]) # Sign-flipped loading can only go x2 in absolute value!
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
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
# We now need to shrink all updates to other loadings, and also constrain PSI2 and PSI5:
con <- c(1:6, 8:10, 12, 15)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(0.25 * delta[c(1:6, 8:10),], # Remaining loadings
                 1 / 15 * k[12, iter] - k[12, iter],
                 1 / 15 * k[15, iter] - k[15, iter])
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
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate

# Iteration 2 ====
iter <- 2
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
# We need to constrain PSI5 again:
con <- 15
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1 / 15 * k[15, iter] - k[15, iter]) # Variances can't go negative
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
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
# We now also need to constrain L22 and PSI1:
con <- c(7, 11)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(-2 * k[7, iter] - k[7, iter], # Maximum x2 for sign-flipped loadings
                 1 / 15 * k[11, iter] - k[11, iter]) # Variances can't go negative
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
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate

# Iteration 3 ====
iter <- 3
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
# PSI5 has gone negative for the third time now, fix to 1e-6:
con <- 15
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1e-6 - k[15, iter]) # Compute delta so PSI5 goes to 1e-6 and can be fixed.
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
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15, ignore PSI5 because it will be fixed anyway
isFixed[15] <- TRUE
k[, iter + 1] <- kUpdate

# Iteration 4 ====
iter <- 4
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate

# Iteration 5 ====
iter <- 5
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.5)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
# # We need to constrain PSI5 to stop it from going negative:
# con <- 15 # PSI5 is the 15th parameter
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(1 / 15 * k[15, iter] - k[15, iter]) # Variances can't go negative
# # Recompute delta for the unconstrained parameters:
# noncon <- setdiff(1:nrow(k), con.total)
# uuisFixed <- isFixed[noncon]
# AIuu <- AI[noncon, noncon, drop = FALSE]
# AIff <- AI[con.total, con.total, drop = FALSE]
# AIuf <- AI[noncon, con.total, drop = FALSE]
# scoresuu <- scores[noncon, 1, drop = FALSE]
# scoresff <- scores[con.total, 1, drop = FALSE]
# deltaff <- delta[con.total, 1, drop = FALSE]
# AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
# uunotFixed <- which(uuisFixed == 0)
# AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
# deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
# delta[noncon, 1] <- deltauu
# # Recompute kUpdate:
# kUpdate <- k[, iter, drop = F] + delta
# cbind(k[, iter, drop = F], kUpdate)
# signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate

# Iteration 6 ====
iter <- 6
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
# # We need to constrain PSI5 again to stop it from going negative:
# con <- 15 # PSI5 is the 15th parameter
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(1 / 15 * k[15, iter] - k[15, iter]) # Variances can't go negative
# # Recompute delta for the unconstrained parameters:
# noncon <- setdiff(1:nrow(k), con.total)
# uuisFixed <- isFixed[noncon]
# AIuu <- AI[noncon, noncon, drop = FALSE]
# AIff <- AI[con.total, con.total, drop = FALSE]
# AIuf <- AI[noncon, con.total, drop = FALSE]
# scoresuu <- scores[noncon, 1, drop = FALSE]
# scoresff <- scores[con.total, 1, drop = FALSE]
# deltaff <- delta[con.total, 1, drop = FALSE]
# AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
# uunotFixed <- which(uuisFixed == 0)
# AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
# deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
# delta[noncon, 1] <- deltauu
# # Recompute kUpdate:
# kUpdate <- k[, iter, drop = F] + delta
# cbind(k[, iter, drop = F], kUpdate)
# signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate

# Iteration 7 ====
iter <- 7
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
# # PSI5 went negative for the third time, we fix it at 1e-6:
# con <- 15 # PSI5 is the 15th parameter
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(1e-6 - k[15, iter]) # Variances can't go negative
# # Recompute delta for the unconstrained parameters:
# noncon <- setdiff(1:nrow(k), con.total)
# uuisFixed <- isFixed[noncon]
# AIuu <- AI[noncon, noncon, drop = FALSE]
# AIff <- AI[con.total, con.total, drop = FALSE]
# AIuf <- AI[noncon, con.total, drop = FALSE]
# scoresuu <- scores[noncon, 1, drop = FALSE]
# scoresff <- scores[con.total, 1, drop = FALSE]
# deltaff <- delta[con.total, 1, drop = FALSE]
# AIuui <- matrix(0, nrow(AIuu), ncol(AIuu))
# uunotFixed <- which(uuisFixed == 0)
# AIuui[uunotFixed, uunotFixed] <- solve(AIuu[uunotFixed, uunotFixed])
# deltauu <- AIuui %*% (scoresuu - AIuf %*% deltaff)
# delta[noncon, 1] <- deltauu
# # Recompute kUpdate:
# kUpdate <- k[, iter, drop = F] + delta
# cbind(k[, iter, drop = F], kUpdate)
# signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15, PSI5 is lower but will be fixed at 1e-6 from now on.
# isFixed[15] <- TRUE
k[, iter + 1] <- kUpdate

# Iteration 8 ====
iter <- 8
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate

# Iteration 9 ====
iter <- 9
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate

# Iteration 10 ====
iter <- 10
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate

# Iteration 11 ====
iter <- 11
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate

# Iteration 12 ====
iter <- 12
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 1)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate

# Iteration 13 ====
iter <- 13
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0, isFixed = isFixed, stepSize = 1)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate

# Iteration 14 ====
iter <- 14
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0, isFixed = isFixed, stepSize = 1)
LL[iter] <- result$LL
deltaLL <- LL[iter] - LL[iter - 1]
delta <- result$delta
AI <- result$AI
scores <- result$scores
kUpdate <- k[, iter, drop = F] + delta
con.total <- c()
cbind(k[, iter, drop = F], kUpdate)
signflip <- matrix(sign(k[1:10, iter]) * sign(kUpdate[1:10,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:10, , drop = F]
(meanMags <- mean(mags[!isFixed[1:10]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[11:16, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate

# Converged!










