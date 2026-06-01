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
  L <- matrix(k[1:28, 1], 14, 2)
  PSI <- diag(k[29:42, 1])
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
  b <- bu[1:14, 1, drop = F]
  u <- bu[15:2842, 1, drop = F]
  e <- y - W %*% bu
  
  {
    dGdL11 <- kronecker(getLoadingDerivs(1, 1, L), K)
    dGdL21 <- kronecker(getLoadingDerivs(2, 1, L), K)
    dGdL31 <- kronecker(getLoadingDerivs(3, 1, L), K)
    dGdL41 <- kronecker(getLoadingDerivs(4, 1, L), K)
    dGdL51 <- kronecker(getLoadingDerivs(5, 1, L), K)
    dGdL61 <- kronecker(getLoadingDerivs(6, 1, L), K)
    dGdL71 <- kronecker(getLoadingDerivs(7, 1, L), K)
    dGdL81 <- kronecker(getLoadingDerivs(8, 1, L), K)
    dGdL91 <- kronecker(getLoadingDerivs(9, 1, L), K)
    dGdL101 <- kronecker(getLoadingDerivs(10, 1, L), K)
    dGdL111 <- kronecker(getLoadingDerivs(11, 1, L), K)
    dGdL121 <- kronecker(getLoadingDerivs(12, 1, L), K)
    dGdL131 <- kronecker(getLoadingDerivs(13, 1, L), K)
    dGdL141 <- kronecker(getLoadingDerivs(14, 1, L), K)
    
    dGdL22 <- kronecker(getLoadingDerivs(2, 2, L), K)
    dGdL32 <- kronecker(getLoadingDerivs(3, 2, L), K)
    dGdL42 <- kronecker(getLoadingDerivs(4, 2, L), K)
    dGdL52 <- kronecker(getLoadingDerivs(5, 2, L), K)
    dGdL62 <- kronecker(getLoadingDerivs(6, 2, L), K)
    dGdL72 <- kronecker(getLoadingDerivs(7, 2, L), K)
    dGdL82 <- kronecker(getLoadingDerivs(8, 2, L), K)
    dGdL92 <- kronecker(getLoadingDerivs(9, 2, L), K)
    dGdL102 <- kronecker(getLoadingDerivs(10, 2, L), K)
    dGdL112 <- kronecker(getLoadingDerivs(11, 2, L), K)
    dGdL122 <- kronecker(getLoadingDerivs(12, 2, L), K)
    dGdL132 <- kronecker(getLoadingDerivs(13, 2, L), K)
    dGdL142 <- kronecker(getLoadingDerivs(14, 2, L), K)
    
    # dGdL33 <- kronecker(getLoadingDerivs(3, 3, L), K)
    # dGdL43 <- kronecker(getLoadingDerivs(4, 3, L), K)
    # dGdL53 <- kronecker(getLoadingDerivs(5, 3, L), K)
    # dGdL63 <- kronecker(getLoadingDerivs(6, 3, L), K)
    # dGdL73 <- kronecker(getLoadingDerivs(7, 3, L), K)
    # dGdL83 <- kronecker(getLoadingDerivs(8, 3, L), K)
    # dGdL93 <- kronecker(getLoadingDerivs(9, 3, L), K)
    # dGdL103 <- kronecker(getLoadingDerivs(10, 3, L), K)
    # dGdL113 <- kronecker(getLoadingDerivs(11, 3, L), K)
    # dGdL123 <- kronecker(getLoadingDerivs(12, 3, L), K)
    # dGdL133 <- kronecker(getLoadingDerivs(13, 3, L), K)
    # dGdL143 <- kronecker(getLoadingDerivs(14, 3, L), K)
    
    dGdPSI1 <- matrix(0, nrow(L), nrow(L)); dGdPSI1[1, 1] <- 1; dGdPSI1 <- kronecker(dGdPSI1, K)
    dGdPSI2 <- matrix(0, nrow(L), nrow(L)); dGdPSI2[2, 2] <- 1; dGdPSI2 <- kronecker(dGdPSI2, K)
    dGdPSI3 <- matrix(0, nrow(L), nrow(L)); dGdPSI3[3, 3] <- 1; dGdPSI3 <- kronecker(dGdPSI3, K)
    dGdPSI4 <- matrix(0, nrow(L), nrow(L)); dGdPSI4[4, 4] <- 1; dGdPSI4 <- kronecker(dGdPSI4, K)
    dGdPSI5 <- matrix(0, nrow(L), nrow(L)); dGdPSI5[5, 5] <- 1; dGdPSI5 <- kronecker(dGdPSI5, K)
    dGdPSI6 <- matrix(0, nrow(L), nrow(L)); dGdPSI6[6, 6] <- 1; dGdPSI6 <- kronecker(dGdPSI6, K)
    dGdPSI7 <- matrix(0, nrow(L), nrow(L)); dGdPSI7[7, 7] <- 1; dGdPSI7 <- kronecker(dGdPSI7, K)
    dGdPSI8 <- matrix(0, nrow(L), nrow(L)); dGdPSI8[8, 8] <- 1; dGdPSI8 <- kronecker(dGdPSI8, K)
    dGdPSI9 <- matrix(0, nrow(L), nrow(L)); dGdPSI9[9, 9] <- 1; dGdPSI9 <- kronecker(dGdPSI9, K)
    dGdPSI10 <- matrix(0, nrow(L), nrow(L)); dGdPSI10[10, 10] <- 1; dGdPSI10 <- kronecker(dGdPSI10, K)
    dGdPSI11 <- matrix(0, nrow(L), nrow(L)); dGdPSI11[11, 11] <- 1; dGdPSI11 <- kronecker(dGdPSI11, K)
    dGdPSI12 <- matrix(0, nrow(L), nrow(L)); dGdPSI12[12, 12] <- 1; dGdPSI12 <- kronecker(dGdPSI12, K)
    dGdPSI13 <- matrix(0, nrow(L), nrow(L)); dGdPSI13[13, 13] <- 1; dGdPSI13 <- kronecker(dGdPSI13, K)
    dGdPSI14 <- matrix(0, nrow(L), nrow(L)); dGdPSI14[14, 14] <- 1; dGdPSI14 <- kronecker(dGdPSI14, K)
    
    Wu <- matrix(0, nrow(d), nrow(k))
    Wu[, 1] <- dGdL11 %*% Gi %*% u
    Wu[, 2] <- dGdL21 %*% Gi %*% u
    Wu[, 3] <- dGdL31 %*% Gi %*% u
    Wu[, 4] <- dGdL41 %*% Gi %*% u
    Wu[, 5] <- dGdL51 %*% Gi %*% u
    Wu[, 6] <- dGdL61 %*% Gi %*% u
    Wu[, 7] <- dGdL71 %*% Gi %*% u
    Wu[, 8] <- dGdL81 %*% Gi %*% u
    Wu[, 9] <- dGdL91 %*% Gi %*% u
    Wu[, 10] <- dGdL101 %*% Gi %*% u
    Wu[, 11] <- dGdL111 %*% Gi %*% u
    Wu[, 12] <- dGdL121 %*% Gi %*% u
    Wu[, 13] <- dGdL131 %*% Gi %*% u
    Wu[, 14] <- dGdL141 %*% Gi %*% u
    
    Wu[, 15] <- 0
    Wu[, 16] <- dGdL22 %*% Gi %*% u
    Wu[, 17] <- dGdL32 %*% Gi %*% u
    Wu[, 18] <- dGdL42 %*% Gi %*% u
    Wu[, 19] <- dGdL52 %*% Gi %*% u
    Wu[, 20] <- dGdL62 %*% Gi %*% u
    Wu[, 21] <- dGdL72 %*% Gi %*% u
    Wu[, 22] <- dGdL82 %*% Gi %*% u
    Wu[, 23] <- dGdL92 %*% Gi %*% u
    Wu[, 24] <- dGdL102 %*% Gi %*% u
    Wu[, 25] <- dGdL112 %*% Gi %*% u
    Wu[, 26] <- dGdL122 %*% Gi %*% u
    Wu[, 27] <- dGdL132 %*% Gi %*% u
    Wu[, 28] <- dGdL142 %*% Gi %*% u
    
    # Wu[, 29] <- 0
    # Wu[, 30] <- 0
    # Wu[, 31] <- dGdL33 %*% Gi %*% u
    # Wu[, 32] <- dGdL43 %*% Gi %*% u
    # Wu[, 33] <- dGdL53 %*% Gi %*% u
    # Wu[, 34] <- dGdL63 %*% Gi %*% u
    # Wu[, 35] <- dGdL73 %*% Gi %*% u
    # Wu[, 36] <- dGdL83 %*% Gi %*% u
    # Wu[, 37] <- dGdL93 %*% Gi %*% u
    # Wu[, 38] <- dGdL103 %*% Gi %*% u
    # Wu[, 39] <- dGdL113 %*% Gi %*% u
    # Wu[, 40] <- dGdL123 %*% Gi %*% u
    # Wu[, 41] <- dGdL133 %*% Gi %*% u
    # Wu[, 42] <- dGdL143 %*% Gi %*% u
    
    Wu[, 29] <- dGdPSI1 %*% Gi %*% u
    Wu[, 30] <- dGdPSI2 %*% Gi %*% u
    Wu[, 31] <- dGdPSI3 %*% Gi %*% u
    Wu[, 32] <- dGdPSI4 %*% Gi %*% u
    Wu[, 33] <- dGdPSI5 %*% Gi %*% u
    Wu[, 34] <- dGdPSI6 %*% Gi %*% u
    Wu[, 35] <- dGdPSI7 %*% Gi %*% u
    Wu[, 36] <- dGdPSI8 %*% Gi %*% u
    Wu[, 37] <- dGdPSI9 %*% Gi %*% u
    Wu[, 38] <- dGdPSI10 %*% Gi %*% u
    Wu[, 39] <- dGdPSI11 %*% Gi %*% u
    Wu[, 40] <- dGdPSI12 %*% Gi %*% u
    Wu[, 41] <- dGdPSI13 %*% Gi %*% u
    Wu[, 42] <- dGdPSI14 %*% Gi %*% u
    
    Wu[, 43] <- Ri %*% e
  }
  
  AI <- t(Wu) %*% P %*% Wu / 2
  rownames(AI) <- colnames(AI) <- rownames(k)
  diag(AI)[1:28] <- diag(AI)[1:28] + percPenalty * diag(AI)[1:28]
  AIi <- matrix(0, nrow(AI), ncol(AI))
  AIi[!isFixed, !isFixed] <- solve(AI[!isFixed, !isFixed])
  
  scores <- matrix(0, nrow(k), 1)
  rownames(scores) <- rownames(k)
  {
    scores[1, 1] <-  -0.5 * (sum(diag(dGdL11 %*% Gi)) - sum(diag((Gi %*% dGdL11 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 1])
    scores[2, 1] <-  -0.5 * (sum(diag(dGdL21 %*% Gi)) - sum(diag((Gi %*% dGdL21 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 2])
    scores[3, 1] <-  -0.5 * (sum(diag(dGdL31 %*% Gi)) - sum(diag((Gi %*% dGdL31 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 3])
    scores[4, 1] <-  -0.5 * (sum(diag(dGdL41 %*% Gi)) - sum(diag((Gi %*% dGdL41 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 4])
    scores[5, 1] <-  -0.5 * (sum(diag(dGdL51 %*% Gi)) - sum(diag((Gi %*% dGdL51 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 5])
    scores[6, 1] <-  -0.5 * (sum(diag(dGdL61 %*% Gi)) - sum(diag((Gi %*% dGdL61 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 6])
    scores[7, 1] <-  -0.5 * (sum(diag(dGdL71 %*% Gi)) - sum(diag((Gi %*% dGdL71 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 7])
    scores[8, 1] <-  -0.5 * (sum(diag(dGdL81 %*% Gi)) - sum(diag((Gi %*% dGdL81 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 8])
    scores[9, 1] <-  -0.5 * (sum(diag(dGdL91 %*% Gi)) - sum(diag((Gi %*% dGdL91 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 9])
    scores[10, 1] <-  -0.5 * (sum(diag(dGdL101 %*% Gi)) - sum(diag((Gi %*% dGdL101 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 10])
    scores[11, 1] <-  -0.5 * (sum(diag(dGdL111 %*% Gi)) - sum(diag((Gi %*% dGdL111 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 11])
    scores[12, 1] <-  -0.5 * (sum(diag(dGdL121 %*% Gi)) - sum(diag((Gi %*% dGdL121 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 12])
    scores[13, 1] <-  -0.5 * (sum(diag(dGdL131 %*% Gi)) - sum(diag((Gi %*% dGdL131 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 13])
    scores[14, 1] <-  -0.5 * (sum(diag(dGdL141 %*% Gi)) - sum(diag((Gi %*% dGdL141 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 14])
    
    scores[15, 1] <-  0
    scores[16, 1] <-  -0.5 * (sum(diag(dGdL22 %*% Gi)) - sum(diag((Gi %*% dGdL22 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 16])
    scores[17, 1] <-  -0.5 * (sum(diag(dGdL32 %*% Gi)) - sum(diag((Gi %*% dGdL32 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 17])
    scores[18, 1] <-  -0.5 * (sum(diag(dGdL42 %*% Gi)) - sum(diag((Gi %*% dGdL42 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 18])
    scores[19, 1] <-  -0.5 * (sum(diag(dGdL52 %*% Gi)) - sum(diag((Gi %*% dGdL52 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 19])
    scores[20, 1] <-  -0.5 * (sum(diag(dGdL62 %*% Gi)) - sum(diag((Gi %*% dGdL62 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 20])
    scores[21, 1] <-  -0.5 * (sum(diag(dGdL72 %*% Gi)) - sum(diag((Gi %*% dGdL72 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 21])
    scores[22, 1] <-  -0.5 * (sum(diag(dGdL82 %*% Gi)) - sum(diag((Gi %*% dGdL82 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 22])
    scores[23, 1] <-  -0.5 * (sum(diag(dGdL92 %*% Gi)) - sum(diag((Gi %*% dGdL92 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 23])
    scores[24, 1] <-  -0.5 * (sum(diag(dGdL102 %*% Gi)) - sum(diag((Gi %*% dGdL102 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 24])
    scores[25, 1] <-  -0.5 * (sum(diag(dGdL112 %*% Gi)) - sum(diag((Gi %*% dGdL112 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 25])
    scores[26, 1] <-  -0.5 * (sum(diag(dGdL122 %*% Gi)) - sum(diag((Gi %*% dGdL122 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 26])
    scores[27, 1] <-  -0.5 * (sum(diag(dGdL132 %*% Gi)) - sum(diag((Gi %*% dGdL132 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 27])
    scores[28, 1] <-  -0.5 * (sum(diag(dGdL142 %*% Gi)) - sum(diag((Gi %*% dGdL142 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 28])
    
    # scores[29, 1] <-  0
    # scores[30, 1] <-  0
    # scores[31, 1] <-  -0.5 * (sum(diag(dGdL33 %*% Gi)) - sum(diag((Gi %*% dGdL33 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 31])
    # scores[32, 1] <-  -0.5 * (sum(diag(dGdL43 %*% Gi)) - sum(diag((Gi %*% dGdL43 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 32])
    # scores[33, 1] <-  -0.5 * (sum(diag(dGdL53 %*% Gi)) - sum(diag((Gi %*% dGdL53 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 33])
    # scores[34, 1] <-  -0.5 * (sum(diag(dGdL63 %*% Gi)) - sum(diag((Gi %*% dGdL63 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 34])
    # scores[35, 1] <-  -0.5 * (sum(diag(dGdL73 %*% Gi)) - sum(diag((Gi %*% dGdL73 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 35])
    # scores[36, 1] <-  -0.5 * (sum(diag(dGdL83 %*% Gi)) - sum(diag((Gi %*% dGdL83 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 36])
    # scores[37, 1] <-  -0.5 * (sum(diag(dGdL93 %*% Gi)) - sum(diag((Gi %*% dGdL93 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 37])
    # scores[38, 1] <-  -0.5 * (sum(diag(dGdL103 %*% Gi)) - sum(diag((Gi %*% dGdL103 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 38])
    # scores[39, 1] <-  -0.5 * (sum(diag(dGdL113 %*% Gi)) - sum(diag((Gi %*% dGdL113 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 39])
    # scores[40, 1] <-  -0.5 * (sum(diag(dGdL123 %*% Gi)) - sum(diag((Gi %*% dGdL123 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 40])
    # scores[41, 1] <-  -0.5 * (sum(diag(dGdL133 %*% Gi)) - sum(diag((Gi %*% dGdL133 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 41])
    # scores[42, 1] <-  -0.5 * (sum(diag(dGdL143 %*% Gi)) - sum(diag((Gi %*% dGdL143 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 42])
  
    scores[29, 1] <- -0.5 * (sum(diag(dGdPSI1 %*% Gi)) - sum(diag((Gi %*% dGdPSI1 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 29])
    scores[30, 1] <- -0.5 * (sum(diag(dGdPSI2 %*% Gi)) - sum(diag((Gi %*% dGdPSI2 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 30])
    scores[31, 1] <- -0.5 * (sum(diag(dGdPSI3 %*% Gi)) - sum(diag((Gi %*% dGdPSI3 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 31])
    scores[32, 1] <- -0.5 * (sum(diag(dGdPSI4 %*% Gi)) - sum(diag((Gi %*% dGdPSI4 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 32])
    scores[33, 1] <- -0.5 * (sum(diag(dGdPSI5 %*% Gi)) - sum(diag((Gi %*% dGdPSI5 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 33])
    scores[34, 1] <- -0.5 * (sum(diag(dGdPSI6 %*% Gi)) - sum(diag((Gi %*% dGdPSI6 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 34])
    scores[35, 1] <- -0.5 * (sum(diag(dGdPSI7 %*% Gi)) - sum(diag((Gi %*% dGdPSI7 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 35])
    scores[36, 1] <- -0.5 * (sum(diag(dGdPSI8 %*% Gi)) - sum(diag((Gi %*% dGdPSI8 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 36])
    scores[37, 1] <- -0.5 * (sum(diag(dGdPSI9 %*% Gi)) - sum(diag((Gi %*% dGdPSI9 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 37])
    scores[38, 1] <- -0.5 * (sum(diag(dGdPSI10 %*% Gi)) - sum(diag((Gi %*% dGdPSI10 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 38])
    scores[39, 1] <- -0.5 * (sum(diag(dGdPSI11 %*% Gi)) - sum(diag((Gi %*% dGdPSI11 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 39])
    scores[40, 1] <- -0.5 * (sum(diag(dGdPSI12 %*% Gi)) - sum(diag((Gi %*% dGdPSI12 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 40])
    scores[41, 1] <- -0.5 * (sum(diag(dGdPSI13 %*% Gi)) - sum(diag((Gi %*% dGdPSI13 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 41])
    scores[42, 1] <- -0.5 * (sum(diag(dGdPSI14 %*% Gi)) - sum(diag((Gi %*% dGdPSI14 %*% Gi) %*% Ci[15:2842, 15:2842])) - t(y) %*% P %*% Wu[, 42])
    
    scores[43, 1] <- -0.5 * (sum(diag(Ri)) - sum(diag((t(W) %*% Ri %*% Ri %*% W) %*% Ci)) - t(y) %*% P %*% Wu[, 43])
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

d <- droplevels(d[d$Man == "LN",])

fit.asr <- asreml(fixed = GY ~ -1 + Env,
                       random = ~ facv(Env, 2):vm(Gen, Ki),
                       residual = ~ idv(units),
                       data = d,
                       maxit = 20)

# Manual:
p <- length(levels(d$Env))
# init <- as.numeric(fit.asr$trace[5:61,1,drop=F])
init <- c(getInitFA(d = d, m = 2, response = "GY", factor = "Env"), var(d$GY) / 10 / 2)
isFixed <- abs(init) < 1e-16
k <- matrix(NA, length(init), 20)
rownames(k) <- c(sprintf("L%d-1", 1:14),
                 sprintf("L%d-2", 1:14),
                 sprintf("PSI%d", 1:14),
                 "rvar")
names(isFixed) <- rownames(k)
LL <- numeric(20)
k[, 1] <- init
k[isFixed, 1] <- 0
X <- model.matrix(GY ~ -1 + Env, d)
Z <- diag(nrow(d))
W <- cbind(X, Z)
y <- matrix(d$GY)

# Iteration 1 ====
iter <- 1
result <- runIteration(k = k[, iter, drop = F], y = y, K = K, Z = Z, X = X, percPenalty = 0.05, isFixed = isFixed, stepSize = 0.25)
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
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15
# We need to constrain L42:
# con <- 18 # L42 is the 18th parameter
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(-2 * k[18, iter] - k[18, iter]) # Sign-flipped loading can only go x2 in absolute value!
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
# signflip <- matrix(sign(k[1:42, iter]) * sign(kUpdate[1:42,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:42, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:42]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[43:57, , drop = F]) # All should be between 1/15 (0.067) and 15
# # We now need to shrink all updates to other loadings, and also constrain PSI3 and PSI8:
# con <- c(1:17, 19:42, 45, 50)
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(0.25 * delta[c(1:17, 19:42),], # Remaining loadings
#                  1 / 15 * k[45, iter] - k[45, iter],
#                  1 / 15 * k[50, iter] - k[50, iter])
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
# signflip <- matrix(sign(k[1:42, iter]) * sign(kUpdate[1:42,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:42, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:42]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[43:56, , drop = F]) # All should be between 1/15 (0.067) and 15
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15
# We need to constrain PSI8 again:
# con <- 50
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(1 / 15 * k[50, iter] - k[50, iter]) # Variances can't go negative
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
# signflip <- matrix(sign(k[1:42, iter]) * sign(kUpdate[1:42,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:42, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:42]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[43:57, , drop = F]) # All should be between 1/15 (0.067) and 15
# # We now also need to constrain PSI3, PSI4, PSI5, PSI6, PSI9, PSI13:
# con <- c(45, 46, 47, 48, 51, 55)
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(1 / 15 * k[45, iter] - k[45, iter], # Variances can't go negative
#                  1 / 15 * k[46, iter] - k[46, iter],
#                  1 / 15 * k[47, iter] - k[47, iter],
#                  1 / 15 * k[48, iter] - k[48, iter],
#                  1 / 15 * k[51, iter] - k[51, iter],
#                  1 / 15 * k[55, iter] - k[55, iter]) # Variances can't go negative
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
# signflip <- matrix(sign(k[1:42, iter]) * sign(kUpdate[1:42,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:42, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:42]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[43:57, , drop = F]) # All should be between 1/15 (0.067) and 15
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15
# PSI3 and PSI8 have both gone negative for the third time now, fix to 1e-6: NOT YET ACTUALLY
# con <- c(45, 50)
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(1 / 15 * k[45, iter] - k[45, iter],
#                  1 / 15 * k[50, iter] - k[50, iter]) # Compute delta so PSI5 goes to 1e-6 and can be fixed.
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
# signflip <- matrix(sign(k[1:42, iter]) * sign(kUpdate[1:42,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:42, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:42]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[43:57, , drop = F]) # All should be between 1/15 (0.067) and 15, ignore PSI5 because it will be fixed anyway
# # We now also need to constrain L33 and PSI13:
# con <- c(31, 55)
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(-2 * k[31, iter] - k[31, iter],
#                  1 / 15 * k[55, iter] - k[55, iter])
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
# signflip <- matrix(sign(k[1:42, iter]) * sign(kUpdate[1:42,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:42, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:42]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[43:57, , drop = F]) # All should be between 1/15 (0.067) and 15, ignore PSI5 because it will be fixed anyway
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15
# We need to constrain PSI3 because it goes up way too much, and PSI8 because it goes negative again. PSI8 will be fixed at 1e-6:
# con <- c(45, 50)
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(15 * k[45, iter] - k[45, iter],
#                  1e-6 - k[50, iter])
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
# signflip <- matrix(sign(k[1:42, iter]) * sign(kUpdate[1:42,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:42, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:42]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[43:57, , drop = F]) # All should be between 1/15 (0.067) and 15, ignore PSI5 because it will be fixed anyway
# # Now PSI13 also needs to be constrained:
# con <- c(55)
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(1 / 15 * k[55, iter] - k[55, iter])
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
# signflip <- matrix(sign(k[1:42, iter]) * sign(kUpdate[1:42,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:42, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:42]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[43:57, , drop = F]) # All should be between 1/15 (0.067) and 15, ignore PSI5 because it will be fixed anyway
# isFixed[50] <- TRUE
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15
# Fourth time that PSI13 has gone negative. We need to fix it to 1e-6:
# con <- c(55)
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(1e-6 - k[55, iter])
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
# signflip <- matrix(sign(k[1:42, iter]) * sign(kUpdate[1:42,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:42, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:42]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[43:57, , drop = F]) # All should be between 1/15 (0.067) and 15, ignore PSI5 because it will be fixed anyway
# # We now also need to constrain PSI3 because it is increasing too fast:
# con <- c(45)
# con.total <- sort(unique(c(con.total, con)))
# delta[con,] <- c(15 * k[45, iter] - k[45, iter])
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
# signflip <- matrix(sign(k[1:42, iter]) * sign(kUpdate[1:42,]) == -1)
# mags <- abs(kUpdate / k[, iter, drop = F])[1:42, , drop = F]
# (meanMags <- mean(mags[!isFixed[1:42]])) # Smaller than 10
# mags[signflip] # All smaller than 2
# (varFactors <- (kUpdate / k[, iter, drop = F])[43:57, , drop = F]) # All should be between 1/15 (0.067) and 15, ignore PSI5 because it will be fixed anyway
# isFixed[55] <- TRUE
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15
# We need to constrain PSI3 and PSI13:
con <- c(31, 41)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1 / 15 * k[31, iter] - k[31, iter],
                 1 / 15 * k[41, iter] - k[41, iter]) # Compute delta so PSI5 goes to 1e-6 and can be fixed.
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
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15, ignore PSI5 because it will be fixed anyway
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15
# Constrain PSI13 again:
con <- c(41)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1 / 15 * k[41, iter] - k[41, iter]) # Compute delta so PSI5 goes to 1e-6 and can be fixed.
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
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15, ignore PSI5 because it will be fixed anyway
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15
# Third time that PSI13 goes negative, fix at 1e-6:
con <- c(41)
con.total <- sort(unique(c(con.total, con)))
delta[con,] <- c(1e-6 - k[41, iter]) # Compute delta so PSI5 goes to 1e-6 and can be fixed.
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
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15, ignore PSI13 because it will be fixed anyway
isFixed["PSI13"] <- TRUE
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15
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
signflip <- matrix(sign(k[1:28, iter]) * sign(kUpdate[1:28,]) == -1)
mags <- abs(kUpdate / k[, iter, drop = F])[1:28, , drop = F]
(meanMags <- mean(mags[!isFixed[1:28]])) # Smaller than 10
mags[signflip] # All smaller than 2
(varFactors <- (kUpdate / k[, iter, drop = F])[29:42, , drop = F]) # All should be between 1/15 (0.067) and 15
k[, iter + 1] <- kUpdate
# Converged!

round(cbind(k[,13], summary(fit.asr)$varcomp[1:43, 1]), 4)
plot(k[,13], summary(fit.asr)$varcomp[1:43, 1])








