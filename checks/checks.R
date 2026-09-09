d <- kmmr::dataDROPS
K <- kmmr::kinshipDROPS
Ki <- solve(K)
attr(Ki, "logDet") <- sum(log(eigen(K, only.values = TRUE)$values))

# corg for 2 environments ====
set.seed(1997)
dss <- droplevels(d[d$Man == "R" & d$Env %in% c("Cam12", "Kar13"),])
dss <- dss[sample(1:nrow(dss), nrow(dss)),] # Shuffling data.
dss$GY[sample(1:nrow(dss), floor(92 / 4))] <- NA # Introducing missing values.

fit <- kmmr(fixed = GY ~ -1 + Env,
            random = ~ re(~ corg(Env):lk(Gen, Ki), vars = ~ Env),
            residual = ~ units,
            data = dss, sparse = FALSE, maxit = 1, debugging = TRUE)

R <- Diagonal(nrow(fit$y), fit$k[4, 1])
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), fit$logDetR) # TRUE
sigma <- Matrix(c(fit$k[1, 1], sqrt(fit$k[1, 1] * fit$k[2, 1]) * fit$k[3, 1],
                  sqrt(fit$k[1, 1] * fit$k[2, 1]) * fit$k[3, 1], fit$k[2, 1]),
                2, 2)
all.equal(sigma, solve(fit$sigmai[[1]])) # TRUE
all.equal(sum(log(eigen(sigma, only.values = TRUE)$values)), fit$logDetvarGtotal[1]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit$logDetnvarGtotal[1]) # TRUE
G <- kronecker(sigma, K)
all.equal(sum(log(eigen(G, only.values = TRUE)$values)), fit$logDetG) # TRUE
V <- fit$Z[[1]] %*% G %*% t(fit$Z[[1]]) + R
Vi <- solve(V)
X <- fit$X
P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi
all.equal(as.numeric(t(fit$y) %*% P %*% fit$y), fit$yPy)
Ri <- solve(R)
Gi <- solve(G)
W <- cbind(fit$X, fit$Z[[1]])
C <- as(t(W) %*% Ri %*% W, "unpackedMatrix")
C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] <- C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] + Gi
all.equal(sum(log(eigen(C, only.values = TRUE)$values)), fit$logDetC) # TRUE
Ci <- solve(C)
Cire1 <- Ci[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)]
tr11 <- sum(Cire1[1:50, 1:50] * Ki)
tr12 <- sum(Cire1[1:50, 51:100] * Ki)
tr22 <- sum(Cire1[51:100, 51:100] * Ki)
all.equal(Matrix(c(tr11, tr12, tr12, tr22), 2, 2), fit$traces[[1]]) # TRUE
Wu.check <- Matrix(0, nrow(fit$y), ncol(fit$Wu), sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit$derivs[[1]])) {
  Wu.check[, param] <- fit$Z[[1]] %*% kronecker(fit$derivs[[1]][[i]], K) %*% Gi %*% fit$u
  param <- param + 1
}
for (i in 1:length(fit$derivs[[2]])) {
  Wu.check[, param] <- fit$derivs[[2]][[i]] %*% Ri %*% fit$e
  param <- param + 1
}
all.equal(Wu.check, fit$Wu) # TRUE
scores.check <- Matrix(0, nrow(fit$scores), 1, sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit$derivs[[1]])) {
  s1 <- sum(diag(kronecker(fit$derivs[[1]][[i]], K) %*% Gi))
  s2 <- sum(diag(Gi %*% (kronecker(fit$derivs[[1]][[i]], K) %*% Gi) %*% Cire1))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
for (i in 1:length(fit$derivs[[2]])) {
  s1 <- sum(diag(fit$derivs[[2]][[i]] %*% Ri))
  s2 <- sum(diag((t(W) %*% Ri %*% fit$derivs[[2]][[i]] %*% Ri %*% W) %*% Ci))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
all.equal(scores.check, fit$scores) # TRUE

# fa for 10 environments ====

# exfa for 10 environments ====

# gk for 10 environments ====

# exfa + gk for 10 environments ====

# Single-stage with corg for 2 environments ====

# Single-stage with exfa for 2 environments ====

# Multi-trait with 2 x corg for 2 traits ====







