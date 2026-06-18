# Final check of the implementation of the dense and extended (sparse) FA models, combined with a dense or FA kinship:
# We do not extend C with y.
# 1: dense x dense
# 2: extended x dense
# 3: dense x FALK
# 4: extended x FALK
library(Matrix)

# Function to rotate loadings matrices to a structure with upper triangle of 0:
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

# Loading in the data:
d <- readRDS("dev/exampledata.rds")
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY")])
d$Units <- factor(1:nrow(d))

# First we check dense by dense ====
# Function for constructing the fa correlation matrix:
formCormatChol <- function(struc = "fa", p = NULL, params = NULL, base = NULL, doChol = TRUE) {
  m <- length(params) / p # Number of factors
  L <- Matrix(params, nrow = p, ncol = m, sparse = FALSE)
  cormat <- L %*% t(L)
  diag(cormat) <- 1
  if (doChol) {
    ret <- list(cormat, chol(cormat))
  } else {
    ret <- list(cormat)
  }
  return(ret)
}

# Function for computing derivatives with respect to loadings:
getLoadingDerivs <- function(indexL, L) {
  rowL <- ((indexL - 1) %% nrow(L)) + 1
  colL <- ((indexL - 1) %/% nrow(L)) + 1
  deriv <- Matrix(0, nrow(L), nrow(L))
  deriv[rowL,] <- deriv[,rowL] <- L[, colL]
  deriv[rowL, rowL] <- 0
  return(deriv)
}

# Dense kinship matrix:
K <- Matrix(readRDS("dev/exampledata_K.rds"))

# Initial values for the fa(Env, ext = TRUE) parameters
vars <- diag(cov(unstack(d, GY ~ Env), use = "complete"))
fit.factanal <- factanal(covmat = cov(unstack(d, GY ~ Env), use = "complete"), factors = 2)
L <- fit.factanal$loadings[,]
L <- Matrix(rotateInitL(L)$L.rotated, sparse = F)
PSI <- Matrix(diag(1 - diag(L %*% t(L))))
resvar <- mean(vars) / 20
genvar <- vars / 20

# Matrix of parameter estimates:
k <- Matrix(0, 14 + 28 + 1, 10, sparse = F)
k[, 1] <- c(genvar, matrix(L), resvar) # k[, 1] <- init # For comparison with asreml
isFixed <- numeric(nrow(k))
isFixed[14 + 15] <- 1

# Design matrices and data vector:
y <- Matrix(d$GY)
X <- sparse.model.matrix(GY ~ -1 + Env, d)
# Where do we need to insert zeros in the design matrix? 1 for no, 2 for yes.
# insertZero <- kronecker(c(rep(1, nrow(L)), rep(2, ncol(L))), c(rep(1, ncol(K)))) > 1
# Z <- matrix(0, nrow(d), (nrow(L) + ncol(L)) * ncol(K))
# Z[, !insertZero] <- model.matrix(GY ~ -1 + Units, d)
# Z <- Matrix(Z)
Z <- sparse.model.matrix(GY ~ -1 + Units, d)

# Preliminaries:
logDetnvarGtotal <- sum(log(eigen(K)$values))
nvarGi <- solve(K)
LL <- numeric(ncol(k))
iter <- 1
sparse <- FALSE # Will C be reasonable sparse?

## Iteration 1 ====
system.time({
Ri <- Matrix(diag(1 / k[nrow(k), iter], nrow(d)))
W <- cbind(X, Z)
if (sparse) {
  C <- t(W) %*% Ri %*% W
} else {
  C <- as(t(W) %*% Ri %*% W, "unpackedMatrix")
}

# Computing Gi:
vars <- sqrt(k[1:14, iter, drop = F]) %*% t(sqrt(k[1:14, iter, drop = F]))
result <- formCormatChol(struc = "fa",
                         p = 14,
                         params = k[15:42, iter],
                         base = NULL,
                         doChol = TRUE)
cormats <- list(result[[1]])
logDetvarG <- 2 * sum(log(diag(result[[2]])))
cormatis <- list(chol2inv(result[[2]]))
sigmai <- (1 / vars) * Reduce(kronecker, cormatis)
logDetvarGtotal <- sum(log(diag(vars)))
logDetvarGtotal <- logDetvarGtotal + logDetvarG * (14 / 14)
Gi <- kronecker(sigmai, nvarGi)
dimGi <- ncol(Gi)
C.start <- ncol(X) + 1
C[C.start:(C.start - 1 + dimGi), C.start:(C.start - 1 + dimGi)] <- C[C.start:(C.start - 1 + dimGi), C.start:(C.start - 1 + dimGi)] + Gi

# (sparse) cholesky decomposition of C:
if (sparse) {
  # Sparse cholesky with permutation
} else {
  Cchol <- chol(C)
}

# Compute the right-hand side of the MME:
rhs <- t(W) %*% Ri %*% y
if (sparse) {
  # Sparse triangular solve
} else {
  z <- Matrix(backsolve(as.matrix(Cchol), rhs, transpose = TRUE))
}

# Compute the LL:
yPy <- as.numeric(t(y) %*% Ri %*% y - t(z) %*% z)
logDetC <- 2 * sum(log(diag(Cchol)))
logDetG <- 246 * logDetvarGtotal + 14 * logDetnvarGtotal
logDetR <- nrow(y) * log(k[nrow(k), iter])
LL[iter] <- -0.5 * (logDetC + logDetR + logDetG + yPy)

# Solving for BLUEs and BLUPs:
bu <- Matrix(backsolve(as.matrix(Cchol), z))
b <- bu[1:ncol(X), 1, drop = FALSE]
u <- bu[(ncol(X) + 1):nrow(bu), 1, drop = FALSE]

# Working variates:
Wu <- Matrix(0, nrow(y), nrow(k), sparse = FALSE)
umat_sigmai <- Matrix(u[1:nrow(u), 1], nrow = 246, ncol = 14) %*% sigmai
derivs <- vector("list", nrow(k) - 1)

# Variances:
for (i in 1:14) {
  if (isFixed[i] == 0) {
    A <- Matrix(0, 14, 14)
    A[i, ] <- 0.5
    derivs[[i]] <- (A + t(A)) * vars * Reduce(kronecker, cormats) / k[i, iter]
  }
}
# Loadings:
L <- Matrix(k[15:42, iter], 14, 2, sparse = FALSE)
for (i in 15:42) {
  if (isFixed[i] == 0) {
    indexL <- i - 14
    derivs[[i]] <- vars * getLoadingDerivs(indexL, L)
  }
}

for (i in 1:42) {
  if (isFixed[i] == 0) {
    Wu[, i] <- Z %*% Matrix((umat_sigmai %*% derivs[[i]])@x)
  }
}

# Residuals:
e <- y - W %*% bu
if (isFixed[nrow(k)] == 0) {
  Wu[, nrow(k)] <- (1 / k[nrow(k), iter]) * e
}

# AI matrix:
XWjxZWj <- t(W) %*% Ri %*% Wu
WiWj <- t(Wu) %*% Ri %*% Wu
buWu <- Matrix(backsolve(Cchol, backsolve(Cchol, XWjxZWj, transpose = TRUE)))
AI <- (WiWj - (t(buWu) %*% XWjxZWj)) / 2
if (any(isFixed == 1)) {
  AIi <- Matrix(0, nrow(AI), ncol(AI), sparse = FALSE, doDiag = FALSE)
  notFixed <- which(isFixed == 0)
  AIi[notFixed, notFixed] <- solve(AI[notFixed, notFixed])
} else {
  AIi <- solve(AI)
}

# Inverse of Ci:
if (sparse) {
  # Sparse inverse subset and permutation
} else {
  Ci <- chol2inv(Cchol)
}

# Random effect scores:
scores <- Matrix(0, nrow(k), 1, sparse = FALSE)
corner <- ncol(X) + 1
traces <- Matrix(0, nrow(sigmai), ncol(sigmai), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigmai)) {
  for (j in 1:ncol(sigmai)) {
    if (j >= i) { # Only do the upper triangle
      rstart <- (corner + (i - 1) * 246)
      rend <- rstart - 1 + 246
      cstart <- (corner + (j - 1) * 246)
      cend <- cstart - 1 + 246
      traces[i, j] <- traces[j, i] <- sum(nvarGi * Ci[rstart:rend, cstart:cend])
    }
  }
}
basisScores <- -0.5 * (246 * sigmai - sigmai %*% traces %*% sigmai - t(umat_sigmai) %*% nvarGi %*% umat_sigmai)
for (i in 1:42) {
  if (isFixed[i] == 0) {
    scores[i, 1] <- sum(derivs[[i]] * basisScores)
  }
}

# Residual score:
if (isFixed[nrow(k)] == 0) {
  evari <- (1 / k[nrow(k), iter])
  S <- Diagonal(nrow(y), 1)
  WSW <- t(W) %*% S %*% W
  scores[nrow(k), 1] <- as.numeric(-0.5 * (sum(S * evari) - sum(evari^2 * Ci[cbind(WSW@i + 1, rep(seq_len(ncol(WSW)), diff(WSW@p)))] * WSW@x) - (evari^2 * t(e) %*% S %*% e)))
}
})

## Checking kmmr against asreml ====
# fit.asr <- asreml(fixed = GY ~ -1 + Env,
#                   random = ~ sfa(Env, 2):vm(Gen, K),
#                   residual = ~ idv(units),
#                   data = d,
#                   maxit = 1)
# init <- c(fit.asr$trace[33:46, 1], fit.asr$trace[5:32, 1], fit.asr$trace[47, 1])
# # Now replace these initial values with our own in k on line 76. Other way around would be easier, but asreml won't do it...
# all.equal(LL[1], fit.asr$trace[1,1])
# all.equal(matrix(fit.asr$score)[c(29:42, 1:28, 43), 1],
#           as.numeric(scores))
# Wu.asr <- as(fit.asr$wrkv[, -ncol(fit.asr$wrkv)], "unpackedMatrix")[, c(29:42, 1:28, 43)]
# all.equal(Wu.asr[, -29], Wu[, -29])

## Saving some stuff to compare agains the sparse x dense model ====
b.dense_dense <- b
u.dense_dense <- u
AI.dense_dense <- AI
Wu.dense_dense <- Wu
scores.dense_dense <- scores
LL.dense_dense <- LL

# Next is the (extended) sparse x dense model ====
# Function for constructing the exfa correlation matrix as well as its inverse:
formCormatChol <- function(struc = "exfa", p = NULL, params = NULL, base = NULL, doChol = TRUE) {
  m <- length(params) / p # Number of factors
  L <- matrix(params, nrow = p, ncol = m)
  LLt <- tcrossprod(L)
  PSIi <- 1 / (1 - diag(LLt))
  diag(LLt) <- 1
  cormat <- forceSymmetric(Matrix(rbind(cbind(LLt, -L),
                                        cbind(-t(L), diag(m))),
                                  sparse = FALSE))
  PSIi_L <- PSIi * L
  cormati <- forceSymmetric(Matrix(rbind(cbind(diag(PSIi), PSIi_L),
                                         cbind(t(PSIi_L), diag(m) + crossprod(L, PSIi_L)))))
  # L <- Matrix(params, nrow = p, ncol = m, sparse = FALSE)
  # cormat <- rbind(cbind(L %*% t(L), -1 * L),
  #                 cbind(-1 * t(L), Diagonal(m, 1)))
  # PSIi <- Diagonal(p, 1 / (1 - diag(cormat)[1:p]))
  # diag(cormat) <- 1
  # cormati <- Matrix(rbind(cbind(PSIi, PSIi %*% L),
  #                         cbind(t(L) %*% PSIi, Diagonal(m, 1) + t(L) %*% PSIi %*% L)))
  if (doChol) {
    ret <- list(cormat, chol(cormat), cormati)
  } else {
    ret <- list(cormat)
  }
  return(ret)
}

# Function for computing derivatives with respect to loadings:
getLoadingDerivs <- function(indexL, p, loadings) {
  m <- length(loadings) / p
  rowL <- ((indexL - 1) %% p) + 1
  colL <- ((indexL - 1) %/% p) + 1
  deriv11 <- matrix(0, p, p)
  # deriv11[rowL,] <- deriv11[, rowL] <- L[, colL]
  deriv11[rowL,] <- deriv11[, rowL] <- loadings[(p * (colL - 1) + 1):(p * colL)]
  deriv11[rowL, rowL] <- 0
  deriv12 <- matrix(0, p, m)
  deriv12[rowL, colL] <- -1
  deriv <- forceSymmetric(Matrix(rbind(cbind(deriv11, deriv12),
                                       cbind(t(deriv12), matrix(0, m, m)))))
  return(deriv)
}

# A simpler wrapper for the Takahashi equations:
takahashi <- function(L) {
  n <- nrow(L)
  d <- diag(L)
  L <- tril(L %*% Diagonal(x = 1 / d), -1)
  d <- d^2
  D <- Diagonal(x = d)
  dp <- diff(L@p)
  jj <- rep(seq_along(dp), dp)
  Zpattern <- sparseMatrix(c(L@i + 1, jj, 1:n), c(jj, L@i + 1, 1:n))
  Z <- forceSymmetric(sparseinv:::.sparseinv_wrapper(L, d, L, Zpattern))
  return(Z)
}

# Dense kinship matrix:
K <- Matrix(readRDS("dev/exampledata_K.rds"))

# Initial values for the fa(Env, ext = TRUE) parameters
vars <- diag(cov(unstack(d, GY ~ Env), use = "complete"))
fit.factanal <- factanal(covmat = cov(unstack(d, GY ~ Env), use = "complete"), factors = 2)
L <- fit.factanal$loadings[,]
L <- Matrix(rotateInitL(L)$L.rotated, sparse = F)
PSI <- Matrix(diag(1 - diag(L %*% t(L))))
resvar <- mean(vars) / 20
genvar <- vars / 20

# Matrix of parameter estimates:
k <- Matrix(0, 14 + 28 + 1, 10, sparse = F)
k[, 1] <- c(genvar, matrix(L), resvar) # k[, 1] <- init # For comparison with asreml
isFixed <- numeric(nrow(k))
isFixed[14 + 15] <- 1

# Design matrices and data vector:
y <- Matrix(d$GY)
X <- sparse.model.matrix(GY ~ -1 + Env, d)
# Where do we need to insert zeros in the design matrix? 1 for no, 2 for yes.
insertZero <- kronecker(c(rep(1, nrow(L)), rep(2, ncol(L))), c(rep(1, ncol(K)))) > 1
Z <- matrix(0, nrow(d), (nrow(L) + ncol(L)) * ncol(K))
Z[, !insertZero] <- model.matrix(GY ~ -1 + Units, d)
Z <- Matrix(Z)

# Preliminaries:
logDetnvarGtotal <- sum(log(eigen(K)$values))
nvarGi <- solve(K)
LL <- numeric(ncol(k))
iter <- 1
sparse <- TRUE # Will C be reasonable sparse?

## Iteration 1 ====
system.time({
Ri <- Matrix(diag(1 / k[nrow(k), iter], nrow(d)))
W <- cbind(X, Z)
if (sparse) {
  C <- forceSymmetric(t(W) %*% Ri %*% W)
} else {
  C <- forceSymmetric(as(t(W) %*% Ri %*% W, "unpackedMatrix")) # Because we will need to add a dense Gi anyway
}

# Computing Gi:
# We need to add unit variances for the factors
vars <- tcrossprod(Matrix(c(sqrt(k[1:14, iter]), rep(1, ncol(L)))))
result <- formCormatChol(struc = "exfa",
                         p = 14,
                         params = k[15:42, iter],
                         base = NULL,
                         doChol = TRUE) # Returns (1) the correlation matrix, (2) the cholesky factor of the correlation matrix, (3) the inverse correlation matrix
cormats <- list(result[[1]])
logDetvarG <- 2 * sum(log(diag(result[[2]])))
cormatis <- list(result[[3]]) # For the exfa structures we explicitly form the inverse in formCormatChol()
sigmai <- (1 / vars) * Reduce(kronecker, cormatis)
# Some code for verification of results:
# Lmat <- matrix(k[15:42, iter], 14, 2)
# PSImat <- diag(1 - diag(Lmat %*% t(Lmat)))
# ses <- diag(sqrt(k[1:14,1]))
# varsmat <- outer(sqrt(c(k[1:14,1], 1, 1)), sqrt(c(k[1:14,1], 1, 1)))
# all.equal(as.matrix(solve(sigmai)), rbind(cbind(ses %*% Lmat %*% t(ses %*% Lmat) + ses %*% PSImat %*% ses, -1 * ses %*% Lmat), cbind(-1 * t(ses %*% Lmat), diag(2))))
# all.equal(as.matrix(sigmai), rbind(cbind(solve(ses %*% PSImat %*% ses), solve(ses %*% PSImat %*% ses) %*% ses %*% Lmat), cbind(t(ses %*% Lmat) %*% solve(ses %*% PSImat %*% ses), diag(2) + t(ses %*% Lmat) %*% solve(ses %*% PSImat %*% ses) %*% ses %*% Lmat)))
# all.equal(as.matrix(solve(sigmai)), varsmat * rbind(cbind(Lmat %*% t(Lmat) + PSImat, -1 * Lmat), cbind(-1 * t(Lmat), diag(2))))
# all.equal(as.matrix(sigmai), (1 / varsmat) * rbind(cbind(solve(PSImat), solve(PSImat) %*% Lmat), cbind(t(Lmat) %*% solve(PSImat), diag(2) + t(Lmat) %*% solve(PSImat) %*% Lmat)))
logDetvarGtotal <- sum(log(diag(vars)))
logDetvarGtotal <- logDetvarGtotal + logDetvarG * (14 / 14)
Gi <- kronecker(sigmai, nvarGi)
dimGi <- ncol(Gi)
C.start <- ncol(X) + 1
if (sparse) {
  # There may be a faster way to do this.
  # Gi <- as(Gi, "TsparseMatrix")
  # Padding Gi with zeros to make it conformable with the sparse base C.
  Gi <- sparseMatrix(i = Gi@i + 1 + C.start - 1,
                     j = rep(seq_len(ncol(Gi)), diff(Gi@p)) + C.start - 1,
                     x = Gi@x,
                     dims = dim(C),
                     symmetric = TRUE)
  C <- C + Gi
} else {
  C[C.start:(C.start - 1 + dimGi), C.start:(C.start - 1 + dimGi)] <- C[C.start:(C.start - 1 + dimGi), C.start:(C.start - 1 + dimGi)] + Gi
}

# (sparse) cholesky decomposition of C:
if (sparse) {
  Cchol <- Cholesky(C, LDL = FALSE, perm = TRUE, super = FALSE)
} else {
  Cchol <- chol(C)
}

# Compute the (permuted) right-hand side of the MME, intermediate result z, and BLUEs/BLUPs:
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

# Compute the LL:
yPy <- as.numeric(crossprod(y, Ri %*% y) - crossprod(z))
if (sparse) {
  logDetC <- determinant(Cchol, sqrt = FALSE, logarithm = TRUE)$modulus
} else {
  logDetC <- 2 * sum(log(diag(Cchol)))
}
logDetG <- nrow(nvarGi) * logDetvarGtotal + nrow(sigmai) * logDetnvarGtotal
logDetR <- nrow(y) * log(k[nrow(k), iter])
LL[iter] <- -0.5 * (logDetC + logDetR + logDetG + yPy)
# Some code for verification of results:
# Vi <- as.matrix(solve(Z %*% solve(kronecker(sigmai, nvarGi)) %*% t(Z) + solve(Ri)))
# Xmat <- as.matrix(X)
# P <- Vi - Vi %*% Xmat %*% solve(t(Xmat) %*% Vi %*% Xmat) %*% t(Xmat) %*% Vi
# t(y) %*% P %*% y

# Splitting BLUEs and BLUPs:
b <- bu[1:ncol(X), 1, drop = FALSE]
u <- bu[(ncol(X) + 1):nrow(bu), 1, drop = FALSE]

# Working variates:
Wu <- Matrix(0, nrow(y), nrow(k), sparse = FALSE)
umat_sigmai <- Matrix(u[1:nrow(u), 1], nrow = nrow(nvarGi), ncol = nrow(sigmai)) %*% sigmai
derivs <- vector("list", nrow(k) - 1)

# Variances:
for (i in 1:14) {
  if (isFixed[i] == 0) {
    A <- Matrix(0, nrow(sigmai), ncol(sigmai))
    A[i, ] <- 0.5
    derivs[[i]] <- forceSymmetric(A + t(A)) * vars * Reduce(kronecker, cormats) / k[i, iter]
  }
}
# Loadings:
for (i in 15:42) {
  if (isFixed[i] == 0) {
    indexL <- i - 14
    derivs[[i]] <- vars * getLoadingDerivs(indexL, p = 14, loadings = k[15:42, iter])
  }
}

for (i in 1:42) {
  if (isFixed[i] == 0) {
    Wu[, i] <- Z %*% Matrix((umat_sigmai %*% derivs[[i]])@x)
  }
}

# Residuals:
e <- y - W %*% bu
if (isFixed[nrow(k)] == 0) {
  Wu[, nrow(k)] <- (1 / k[nrow(k), iter]) * e
}

# AI matrix:
WiWj <- forceSymmetric(crossprod(Wu, Ri %*% Wu)) # Similar to t(y) %*% Ri %*% y
if (sparse) {
  XWjxZWj <- solve(Cchol, t(W) %*% Ri %*% Wu, system = "P") # Similar to permRHS
  zWu <- solve(Cchol, XWjxZWj, system = "L") # Similar to z
  buWu <- solve(Cchol, zWu, system = "Lt") # Similar to obtaining bu using a backsolve
} else {
  XWjxZWj <- t(W) %*% Ri %*% Wu
  zWu <- backsolve(Cchol, XWjxZWj, transpose = TRUE)
  buWu <- Matrix(backsolve(Cchol, zWu))
}
buWuXWjxZWj <- forceSymmetric(t(buWu) %*% XWjxZWj)
AI <- (WiWj - buWuXWjxZWj) / 2
if (any(isFixed == 1)) {
  AIi <- Matrix(0, nrow(AI), ncol(AI), sparse = FALSE, doDiag = FALSE)
  notFixed <- which(isFixed == 0)
  AIi[notFixed, notFixed] <- solve(AI[notFixed, notFixed])
} else {
  AIi <- solve(AI)
}

# Inverse of Ci:
if (sparse) {
  # Sparse inverse subset and changing it back to model order:
  inv_Perm <- invPerm(Cchol@perm, zero.p = TRUE, zero.res = FALSE)
  Ci <- takahashi(as(Cchol, "dtCMatrix"))[inv_Perm, inv_Perm]
} else {
  Ci <- chol2inv(Cchol)
}

# Random effect scores:
scores <- Matrix(0, nrow(k), 1, sparse = FALSE)
corner <- ncol(X) + 1
traces <- Matrix(0, nrow(sigmai), ncol(sigmai), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigmai)) {
  for (j in 1:ncol(sigmai)) {
    if (j >= i) { # Only do the upper triangle
      if (!(sigmai[i, j] == 0)) {
        rstart <- (corner + (i - 1) * nrow(nvarGi))
        rend <- rstart - 1 + nrow(nvarGi)
        cstart <- (corner + (j - 1) * nrow(nvarGi))
        cend <- cstart - 1 + nrow(nvarGi)
        traces[i, j] <- sum(nvarGi * Ci[rstart:rend, cstart:cend]) # The indexing of Ci is a little bit slow, could be faster...
      }
    }
  }
}
traces <- forceSymmetric(traces)
basisScores <- -0.5 * (nrow(nvarGi) * sigmai - sigmai %*% traces %*% sigmai - t(umat_sigmai) %*% nvarGi %*% umat_sigmai)
for (i in 1:42) {
  if (isFixed[i] == 0) {
    scores[i, 1] <- sum(derivs[[i]] * basisScores)
  }
}

# Residual score:
if (isFixed[nrow(k)] == 0) {
  evari <- (1 / k[nrow(k), iter])
  S <- Diagonal(nrow(y), 1)
  WSW <- t(W) %*% S %*% W
  scores[nrow(k), 1] <- as.numeric(-0.5 * (sum(S * evari) - sum(evari^2 * Ci[cbind(WSW@i + 1, rep(seq_len(ncol(WSW)), diff(WSW@p)))] * WSW@x) - (evari^2 * t(e) %*% S %*% e)))
}
})

## Comparing to the saved results from the dense x dense model ====
# All TRUE!
all.equal(as.matrix(b.dense_dense), as.matrix(b), check.attributes = FALSE)
all.equal(as.matrix(u.dense_dense), as.matrix(u[1:3444, 1, drop = FALSE]), check.attributes = FALSE)
all.equal(as.matrix(AI.dense_dense), as.matrix(AI), check.attributes = FALSE)
all.equal(as.matrix(Wu.dense_dense), as.matrix(Wu), check.attributes = FALSE)
all.equal(LL.dense_dense, LL, check.attributes = FALSE)

# Next is the dense x sparse model ====













