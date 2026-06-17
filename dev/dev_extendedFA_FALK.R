library(Matrix)

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

d <- readRDS("dev/exampledata.rds")
K <- readRDS("dev/exampledata_K.rds")
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY")])

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
# lP <- list(spam::as.spam(C.extended))
# test <- LMMsolver:::ADchol(lP)
# test2 <- sparseMatrix()

# Checking score calculations:
corner <- ncol(X) + 1
traces <- matrix(0, nrow(sigmai.extended), ncol(sigmai.extended))
dim.nvar <- ncol(FALKi.extended)
i <- 1
j <- 1
for (i in 1:nrow(sigmai.extended)) {
  for (j in 1:ncol(sigmai.extended)) {
    if (j >= i) { # Only do the upper triangle
      rstart <- (corner + (i - 1) * dim.nvar)
      rend <- rstart - 1 + dim.nvar
      cstart <- (corner + (j - 1) * dim.nvar)
      cend <- cstart - 1 + dim.nvar
      traces[i, j] <- sum(FALKi.extended * Ci.extended.sparse[rstart:rend, cstart:cend])
    }
  }
}







