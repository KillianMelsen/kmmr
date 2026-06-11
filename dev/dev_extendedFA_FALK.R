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
fa.K <- factanal(covmat = K, factors = 5)
LK <- diag(sqrt(diag(K))) %*% rotateInitL(fa.K$loadings[,])$L.rotated
PSIK <- outer(sqrt(diag(K)), sqrt(diag(K))) * diag(fa.K$uniquenesses)
FALK <- rbind(cbind(LK %*% t(LK) + PSIK, -1 * LK),
              cbind(-1 * t(LK), diag(ncol(LK))))
FALKi2 <- solve(FALK); FALKi2[abs(FALKi2) < 1e-13] <- 0
FALKi <- rbind(cbind(solve(PSIK), solve(PSIK) %*% LK),
               cbind(t(LK) %*% solve(PSIK), diag(ncol(LK)) + t(LK) %*% solve(PSIK) %*% LK))

vars <- diag(cov(unstack(d, GY ~ Env), use = "complete"))
fit.factanal <- factanal(covmat = cov(unstack(d, GY ~ Env), use = "complete"), factors = 2)
L <- fit.factanal$loadings[,]
L <- rotateInitL(L)$L.rotated
PSI <- diag(1 - diag(L %*% t(L)))
resvar <- mean(vars) / 20
genvar <- vars / 20

y <- matrix(d$GY)
X <- model.matrix(GY ~ -1 + Env, d)
Z.naive <- diag(nrow(d))
# Z.extended <- cbind(Z.naive, matrix(0, nrow(Z.naive), 2 * nrow(K)))
Ri <- Matrix(diag(1 / resvar, nrow(d)))
covL <- diag(sqrt(genvar)) %*% L
covPSI <- diag(genvar) %*% PSI
sigma.naive <- covL %*% t(covL) + covPSI
sigma.extended <- rbind(cbind(covL %*% t(covL) + covPSI,    -1 * covL),
                        cbind(-1 * t(covL),                 diag(2)))
sigmai.naive <- solve(sigma.naive)
sigmai.extended <- solve(sigma.extended); sigmai.extended[abs(sigmai.extended) < 1e-13] <- 0
Gi.naive <- kronecker(sigmai.naive, FALKi)
Gi.extended <- kronecker(sigmai.extended, FALKi)

W.naive <- Matrix(cbind(X, Z.naive))
W.extended <- Matrix(cbind(X, Z.extended))

C.naive <- as.matrix(t(W.naive) %*% Ri %*% W.naive)
C.naive[15:nrow(C.naive), 15:ncol(C.naive)] <- C.naive[15:nrow(C.naive), 15:ncol(C.naive)] + Gi.naive
C.extended <- as.matrix(t(W.extended) %*% Ri %*% W.extended)
C.extended[15:nrow(C.extended), 15:ncol(C.extended)] <- C.extended[15:nrow(C.extended), 15:ncol(C.extended)] + Gi.extended

# Naive version:
C.naive.spam <- as.spam(C.naive)
chC.naive <- chol(C.naive.spam)
display(C.naive.spam)
display(chC.naive)
summary(C.naive.spam)
summary(chC.naive)

# Extended version:
C.extended.spam <- as.spam(C.extended)
chC.extended <- chol(C.extended.spam)
display(C.extended.spam)
display(chC.extended)
summary(C.extended.spam)
summary(chC.extended)





RHS.naive <- t(W.naive) %*% Ri %*% y
RHS.extended <- t(W.extended) %*% Ri %*% y

Ci.naive <- solve(C.naive)
Ci.extended <- solve(C.extended)

bu.naive <- Ci.naive %*% RHS.naive
bu.extended <- Ci.extended %*% RHS.extended

u.naive <- bu.naive[15:nrow(bu.naive),, drop = F]
u.extended <- bu.extended[15:nrow(bu.extended),, drop = F]

all.equal(u.naive[, 1], u.extended[1:nrow(u.naive), 1])

U.naive <- chol(C.naive)
U.extended <- chol(C.extended)

system.time({
Gi.extended.sparse <- Matrix(Gi.extended)
U.extended.sparse.pivoted <- Matrix::Cholesky(Gi.extended.sparse, LDL = FALSE, perm = FALSE, super = FALSE)
tmp <- as(U.extended.sparse.pivoted, "dtCMatrix")
permMat <- as(diag(ncol(Gi.extended.sparse)), "dgCMatrix")
Gii.extended.sparse <- sparseinv::Takahashi_Davis(Q = Gi.extended.sparse,
                                                  cholQp = tmp,
                                                  P = permMat)
})

system.time(solve(Gi.extended))













C.extended.sparse <- Matrix::Matrix(C.extended)
U.extended.sparse <- as.matrix(Matrix::chol(C.extended.sparse))
system.time({
U.extended.sparse.pivoted <- Matrix::Cholesky(C.extended.sparse, LDL = FALSE, perm = FALSE, super = FALSE)
tmp <- as(U.extended.sparse.pivoted, "dtCMatrix")
permMat <- as(diag(ncol(C.extended.sparse)), "dgCMatrix")
Ci.extended.sparse <- sparseinv::Takahashi_Davis(Q = C.extended.sparse,
                                                 cholQp = tmp,
                                                 P = permMat)
})
tmp2 <- as.matrix(Ci.extended.sparse)
tmp3 <- Ci.extended

