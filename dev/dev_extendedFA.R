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
# load("~/repos/glfBLUP_2025/genotypes/K_hyper.RData")
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY")])
# Ki <- solve(as.matrix(Matrix::nearPD(K)$mat))
Ki <- solve(K)
vars <- diag(cov(unstack(d, GY ~ Env), use = "complete"))
fit.factanal <- factanal(covmat = cov(unstack(d, GY ~ Env), use = "complete"), factors = 2)
L <- fit.factanal$loadings[,]
L <- rotateInitL(L)$L.rotated
PSI <- diag(1 - diag(L %*% t(L)))
resvar <- mean(vars) / 20
genvar <- vars / 20

y <- matrix(d$GY)
# d <- expand.grid(Gen = rownames(K),
#                  Env = paste0("E", 1:14),
#                  GY = 1) # Dummy dataframe
X <- model.matrix(GY ~ -1 + Env, d)
Z.naive <- diag(nrow(d))
Z.extended <- cbind(Z.naive, matrix(0, nrow(Z.naive), 2 * nrow(K)))
Ri <- Matrix(diag(1 / resvar, nrow(d)))
covL <- diag(sqrt(genvar)) %*% L
covPSI <- diag(genvar) %*% PSI
sigma.naive <- covL %*% t(covL) + covPSI
sigma.extended <- rbind(cbind(covL %*% t(covL) + covPSI,    -1 * covL),
                        cbind(-1 * t(covL),                 diag(2)))
sigmai.naive <- solve(sigma.naive)
sigmai.extended <- solve(sigma.extended); sigmai.extended[abs(sigmai.extended) < 1e-13] <- 0
Gi.naive <- kronecker(sigmai.naive, Ki)
Gi.extended <- kronecker(sigmai.extended, Ki)

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

