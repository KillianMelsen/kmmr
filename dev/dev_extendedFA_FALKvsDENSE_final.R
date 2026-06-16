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
k[, 1] <- c(genvar, matrix(L), resvar)
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
fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ sfa(Env, 2):vm(Gen, K),
                  residual = ~ idv(units),
                  data = d,
                  maxit = 1)
init <- c(fit.asr$trace[33:46, 1], fit.asr$trace[5:32, 1], fit.asr$trace[47, 1])
# Now replace these initial values with our own in k on line 76. Other way around would be easier, but asreml won't do it...
all.equal(LL[1], fit.asr$trace[1,1])
all.equal(matrix(fit.asr$score)[c(29:42, 1:28, 43), 1],
          as.numeric(scores))
Wu.asr <- as(fit.asr$wrkv[, -ncol(fit.asr$wrkv)], "unpackedMatrix")[, c(29:42, 1:28, 43)]
all.equal(Wu.asr[, -29], Wu[, -29])

# Next is the (extended) sparse x dense model ====
# Function for constructing the fa correlation matrix as well as its inverse:
formCormatChol <- function(struc = "exfa", p = NULL, params = NULL, base = NULL, doChol = TRUE) {
  m <- length(params) / p # Number of factors
  L <- Matrix(params, nrow = p, ncol = m, sparse = FALSE)
  cormat <- rbind(cbind(L %*% t(L), -1 * L),
                  cbind(-1 * t(L), Diagonal(m, 1)))
  PSIi <- Diagonal(p, 1 / (1 - diag(cormat)[1:p]))
  diag(cormat) <- 1
  cormati <- Matrix(rbind(cbind(PSIi, PSIi %*% L),
                          cbind(t(L) %*% PSIi, Diagonal(m, 1) + t(L) %*% PSIi %*% L)))
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


















# Making naive and extended versions of the genetic FA structure:
fa.K <- factanal(covmat = K, factors = 20)
LK <- diag(sqrt(diag(K))) %*% rotateInitL(fa.K$loadings[,])$L.rotated
PSIK <- outer(sqrt(diag(K)), sqrt(diag(K))) * diag(fa.K$uniquenesses)
# FALK <- rbind(cbind(LK %*% t(LK) + PSIK, -1 * LK),
#               cbind(-1 * t(LK), diag(ncol(LK))))
# FALKi2 <- solve(FALK); FALKi2[abs(FALKi2) < 1e-13] <- 0
FALK.naive <- LK %*% t(LK) + PSIK
FALKi.naive <- solve(FALK.naive)
FALKi.extended <- rbind(cbind(solve(PSIK), solve(PSIK) %*% LK),
                        cbind(t(LK) %*% solve(PSIK), diag(ncol(LK)) + t(LK) %*% solve(PSIK) %*% LK))

# Making naive and extended versions of the environmental FA structure:
vars <- diag(cov(unstack(d, GY ~ Env), use = "complete"))
fit.factanal <- factanal(covmat = cov(unstack(d, GY ~ Env), use = "complete"), factors = 2)
L <- fit.factanal$loadings[,]
L <- rotateInitL(L)$L.rotated
PSI <- diag(1 - diag(L %*% t(L)))
resvar <- mean(vars) / 20
genvar <- vars / 20
covL <- diag(sqrt(genvar)) %*% L
covPSI <- diag(genvar) %*% PSI
sigma.naive <- covL %*% t(covL) + covPSI
# sigma.extended <- rbind(cbind(covL %*% t(covL) + covPSI,    -1 * covL),
#                         cbind(-1 * t(covL),                 diag(ncol(covL))))
sigmai.naive <- solve(sigma.naive)
sigmai.extended <- rbind(cbind(solve(covPSI),                solve(covPSI) %*% covL),
                         cbind(t(covL) %*% solve(covPSI),    diag(ncol(covL)) + t(covL) %*% solve(covPSI) %*% covL))

# Combinine the environmental and genetic structures to form G:
Gi.naive <- kronecker(sigmai.naive, FALKi.naive)
Gi.extended <- kronecker(sigmai.extended, FALKi.extended)

# Design matrices:
y <- matrix(d$GY)
X <- model.matrix(GY ~ -1 + Env, d)
Z.naive <- diag(nrow(d))
insertZero <- kronecker(c(rep(1, nrow(L)), rep(2, ncol(L))),
                        c(rep(1, nrow(LK)), rep(2, ncol(LK)))) > 1
# Z.extended <- cbind(Z.naive, matrix(0, nrow(Z.naive), 2 * nrow(K)))
Z.extended <- matrix(0, nrow(Z.naive), length(insertZero))
Z.extended[, !insertZero] <- Z.naive
W.naive <- Matrix(cbind(X, Z.naive))
W.extended <- Matrix(cbind(X, Z.extended))

# Residual covariance matrix:
Ri <- Matrix(diag(1 / resvar, nrow(d)))

# Coefficient matrices
C.naive <- as.matrix(t(W.naive) %*% Ri %*% W.naive)
C.naive[15:nrow(C.naive), 15:ncol(C.naive)] <- C.naive[15:nrow(C.naive), 15:ncol(C.naive)] + Gi.naive
C.extended <- as.matrix(t(W.extended) %*% Ri %*% W.extended)
C.extended[15:nrow(C.extended), 15:ncol(C.extended)] <- C.extended[15:nrow(C.extended), 15:ncol(C.extended)] + Gi.extended

# Naive version:
# C.naive.spam <- spam::as.spam(C.naive)
# chC.naive <- spam::chol(C.naive.spam)
# spam::display(C.naive.spam)
# spam::display(chC.naive)
# spam::summary(C.naive.spam)
# spam::summary(chC.naive)

# Extended version:
# C.extended.spam <- spam::as.spam(C.extended)
# chC.extended <- spam::chol(C.extended.spam)
# spam::display(C.extended.spam)
# spam::display(chC.extended)
# spam::summary(C.extended.spam)
# spam::summary(chC.extended)

# Right hand sides of the MME:
RHS.naive <- t(W.naive) %*% Ri %*% y
RHS.extended <- t(W.extended) %*% Ri %*% y

# Full, simple inverses:
Ci.naive <- solve(C.naive)
Ci.extended <- solve(C.extended)

# Solutions (BLUEs/BLUPs)
bu.naive <- Ci.naive %*% RHS.naive
bu.extended <- Ci.extended %*% RHS.extended

# BLUPs only:
u.naive <- as.matrix(bu.naive[15:nrow(bu.naive),, drop = F])
u.extended <- as.matrix(bu.extended[15:nrow(bu.extended),, drop = F])
u.extended.nofactorscores <- u.extended[!insertZero,, drop = F]

all.equal(u.naive, u.extended.nofactorscores)

U.naive <- chol(C.naive)
U.extended <- chol(C.extended)

system.time({
C.extended.sparse <- Matrix(C.extended)
U.extended.sparse <- Matrix::Cholesky(C.extended.sparse, LDL = FALSE, perm = TRUE, super = FALSE)
tmp <- as(U.extended.sparse, "dtCMatrix") # Lower triangular
permMat <- Matrix::Diagonal(ncol(C.extended.sparse))[U.extended.sparse@perm + 1,]
Ci.extended.sparse <- sparseinv::Takahashi_Davis(Q = C.extended.sparse,
                                                 cholQp = tmp,
                                                 P = Matrix::Diagonal(ncol(C.extended.sparse)),
                                                 return_perm_chol = 0)
})

left <- permMat %*% C.extended.sparse %*% t(permMat)
right <- tmp %*% t(tmp)
all.equal(as(left, "generalMatrix"), right)

# Putting the sparse inverse coefficient matrix back into normal order:
Ci.extended.sparse <- t(permMat) %*% Ci.extended.sparse %*% permMat

# Naive alternative:
system.time(Ci.extended <- solve(C.extended))

# Some checks.
# 1. All elements that are non-zero in the coefficient matrix MUST be non-zero in the sparse inverse.
#    As all derivatives of C should have the same sparsity pattern, this ensures all traces can be computed
#    using Hadamard products.
all(Ci.extended.sparse[C.extended.sparse != 0] != 0)
# 2. Percentage of non-zero elements:
nnzero(Ci.extended.sparse) / prod(dim(Ci.extended.sparse))
# 3. All non-zero elements in the sparse inverse should be the same as the corresponding elements
#    in the dense inverse:
all.equal(Ci.extended.sparse[which(Ci.extended.sparse != 0)], Ci.extended[which(Ci.extended.sparse != 0)])

# Now let's see if we can achieve the same using LMMsolver's internal functions.
lP <- list(spam::as.spam(C.extended))
test <- LMMsolver:::ADchol(lP)
test2 <- sparseMatrix()









