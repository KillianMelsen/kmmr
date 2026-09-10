d <- kmmr::dataDROPS
K <- kmmr::kinshipDROPS
Ki <- solve(K)
attr(Ki, "logDet") <- sum(log(eigen(K, only.values = TRUE)$values))

# corg for 2 environments ====
set.seed(1997)
dss <- droplevels(d[d$Man == "R" & d$Env %in% c("Cam12", "Kar13"),])
dss <- dss[sample(1:nrow(dss), nrow(dss)),] # Shuffling data.
dss$GY[sample(1:nrow(dss), floor(nrow(dss) / 4))] <- NA # Introducing missing values.

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
all.equal(as.numeric(t(fit$y) %*% P %*% fit$y), fit$yPy) # TRUE
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
rm(list = setdiff(ls(), c("d", "K", "Ki")))
set.seed(1997)
dss <- droplevels(d[d$Man == "R" & d$Env %in% levels(d$Env)[1:10],])
dss <- dss[sample(1:nrow(dss), nrow(dss)),] # Shuffling data.
dss$GY[sample(1:nrow(dss), floor(nrow(dss) / 4))] <- NA # Introducing missing values.

fit <- kmmr(fixed = GY ~ -1 + Env,
            random = ~ re(~ fa(Env, 1):lk(Gen, Ki), vars = ~ Env),
            residual = ~ units,
            data = dss, sparse = FALSE, maxit = 1, debugging = TRUE)

R <- Diagonal(nrow(fit$y), fit$k[21, 1])
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), fit$logDetR) # TRUE
vars <- outer(sqrt(fit$k[1:10, 1]), sqrt(fit$k[1:10, 1]))
all.equal(vars, as.matrix(fit$vars[[1]])) # TRUE
L <- Matrix(fit$k[11:20, 1], length(levels(dss$Env)), 1)
sigma <- forceSymmetric(vars * (L %*% t(L) + diag(1 - diag(L %*% t(L)))))
all.equal(sigma, solve(fit$sigmai[[1]])) # TRUE
all.equal(sum(log(eigen(sigma, only.values = TRUE)$values)), fit$logDetvarGtotal[1]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit$logDetnvarGtotal[1]) # TRUE
G <- kronecker(sigma, K)
all.equal(sum(log(eigen(G, only.values = TRUE)$values)), fit$logDetG) # TRUE
V <- fit$Z[[1]] %*% G %*% t(fit$Z[[1]]) + R
Vi <- solve(V)
X <- fit$X
P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi
all.equal(as.numeric(t(fit$y) %*% P %*% fit$y), fit$yPy) # TRUE
Ri <- solve(R)
Gi <- solve(G)
W <- cbind(fit$X, fit$Z[[1]])
C <- as(t(W) %*% Ri %*% W, "unpackedMatrix")
C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] <- C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] + Gi
all.equal(sum(log(eigen(C, only.values = TRUE)$values)), fit$logDetC) # TRUE
Ci <- solve(C)
Cire1 <- Ci[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)]
traces <- Matrix(0, nrow(sigma), ncol(sigma), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigma)) {
  for (j in i:ncol(sigma)) {
    traces[i, j] <- sum(Cire1[((i - 1) * nrow(Ki) + 1):(i * nrow(Ki)), ((j - 1) * nrow(Ki) + 1):(j * nrow(Ki))] * Ki)
  }
}
traces <- forceSymmetric(traces)
all.equal(traces, fit$traces[[1]]) # TRUE
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

# exfa for 10 environments ====
rm(list = setdiff(ls(), c("d", "K", "Ki", "fit")))
set.seed(1997)
dss <- droplevels(d[d$Man == "R" & d$Env %in% levels(d$Env)[1:10],])
dss <- dss[sample(1:nrow(dss), nrow(dss)),] # Shuffling data.
dss$GY[sample(1:nrow(dss), floor(nrow(dss) / 4))] <- NA # Introducing missing values.

fit2 <- kmmr(fixed = GY ~ -1 + Env,
             random = ~ re(~ exfa(Env, 1):lk(Gen, Ki), vars = ~ Env),
             residual = ~ units,
             data = dss, sparse = TRUE, maxit = 1, debugging = TRUE, denseThreshold = 0.2)

R <- Diagonal(nrow(fit2$y), fit2$k[21, 1])
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), fit2$logDetR) # TRUE
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), fit$logDetR) # TRUE
vars <- outer(sqrt(c(fit2$k[1:10, 1], 1)), sqrt(c(fit2$k[1:10, 1], 1)))
all.equal(vars, as.matrix(fit2$vars[[1]])) # TRUE
L <- Matrix(fit2$k[11:20, 1], length(levels(dss$Env)), 1)
sigma <- vars * rbind(cbind(L %*% t(L) + diag(1 - diag(L %*% t(L))), -L), cbind(-t(L), Matrix(1, 1, 1)))
sigma <- forceSymmetric(sigma)
all.equal(sigma, as(solve(fit2$sigmai[[1]]), "unpackedMatrix")) # TRUE
all.equal(sum(log(eigen(sigma, only.values = TRUE)$values)), fit2$logDetvarGtotal[1]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit2$logDetnvarGtotal[1]) # TRUE
G <- kronecker(sigma, K)
all.equal(sum(log(eigen(G, only.values = TRUE)$values)), fit2$logDetG) # TRUE
V <- fit2$Z[[1]] %*% G %*% t(fit2$Z[[1]]) + R
Vi <- solve(V)
X <- fit2$X
P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi
all.equal(as.numeric(t(fit2$y) %*% P %*% fit2$y), fit2$yPy) # TRUE
Ri <- solve(R)
Gi <- solve(G)
W <- cbind(fit2$X, fit2$Z[[1]])
C <- as(t(W) %*% Ri %*% W, "unpackedMatrix")
C[(ncol(fit2$X) + 1):nrow(C), (ncol(fit2$X) + 1):ncol(C)] <- C[(ncol(fit2$X) + 1):nrow(C), (ncol(fit2$X) + 1):ncol(C)] + Gi
all.equal(sum(log(eigen(C, only.values = TRUE)$values)), as.numeric(fit2$logDetC)) # TRUE
Ci <- solve(C)
Cire1 <- Ci[(ncol(fit2$X) + 1):nrow(C), (ncol(fit2$X) + 1):ncol(C)]
traces <- Matrix(0, nrow(sigma), ncol(sigma), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigma)) {
  for (j in i:ncol(sigma)) {
    traces[i, j] <- sum(Cire1[((i - 1) * nrow(Ki) + 1):(i * nrow(Ki)), ((j - 1) * nrow(Ki) + 1):(j * nrow(Ki))] * Ki)
  }
}
traces <- forceSymmetric(traces)
all.equal(traces[fit2$traces[[1]] != 0], fit2$traces[[1]][fit2$traces[[1]] != 0]) # TRUE
Wu.check <- Matrix(0, nrow(fit2$y), ncol(fit2$Wu), sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit2$derivs[[1]])) {
  Wu.check[, param] <- fit2$Z[[1]] %*% kronecker(fit2$derivs[[1]][[i]], K) %*% Gi %*% fit2$u
  param <- param + 1
}
for (i in 1:length(fit2$derivs[[2]])) {
  Wu.check[, param] <- fit2$derivs[[2]][[i]] %*% Ri %*% fit2$e
  param <- param + 1
}
all.equal(Wu.check, fit2$Wu) # TRUE
scores.check <- Matrix(0, nrow(fit2$scores), 1, sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit2$derivs[[1]])) {
  s1 <- sum(diag(kronecker(fit2$derivs[[1]][[i]], K) %*% Gi))
  s2 <- sum(diag(Gi %*% (kronecker(fit2$derivs[[1]][[i]], K) %*% Gi) %*% Cire1))
  s3 <- as.numeric(t(fit2$y) %*% P %*% fit2$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
for (i in 1:length(fit2$derivs[[2]])) {
  s1 <- sum(diag(fit2$derivs[[2]][[i]] %*% Ri))
  s2 <- sum(diag((t(W) %*% Ri %*% fit2$derivs[[2]][[i]] %*% Ri %*% W) %*% Ci))
  s3 <- as.numeric(t(fit2$y) %*% P %*% fit2$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
all.equal(scores.check, fit$scores) # TRUE

# gk for 10 environments ====
rm(list = setdiff(ls(), c("d", "K", "Ki")))
ED <- kmmr::distmatDROPS
set.seed(1997)
dss <- droplevels(d[d$Man == "R" & d$Env %in% levels(d$Env)[1:10],])
dss <- dss[sample(1:nrow(dss), nrow(dss)),] # Shuffling data.
dss$GY[sample(1:nrow(dss), floor(nrow(dss) / 4))] <- NA # Introducing missing values.
ED <- ED[levels(dss$Env), levels(dss$Env)]

fit <- kmmr(fixed = GY ~ -1 + Env,
            random = ~ re(~ gk(Env, ED):lk(Gen, Ki), vars = ~ Env),
            residual = ~ units,
            data = dss, sparse = FALSE, maxit = 1, debugging = TRUE)

R <- Diagonal(nrow(fit$y), fit$k[12, 1])
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), fit$logDetR) # TRUE
vars <- outer(sqrt(fit$k[1:10, 1]), sqrt(fit$k[1:10, 1]))
all.equal(vars, as.matrix(fit$vars[[1]])) # TRUE
sigma <- forceSymmetric(vars * exp(-fit$k[11, 1] * ED))
all.equal(sigma, solve(fit$sigmai[[1]])) # TRUE
all.equal(sum(log(eigen(sigma, only.values = TRUE)$values)), fit$logDetvarGtotal[1]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit$logDetnvarGtotal[1]) # TRUE
G <- kronecker(sigma, K)
all.equal(sum(log(eigen(G, only.values = TRUE)$values)), fit$logDetG) # TRUE
V <- fit$Z[[1]] %*% G %*% t(fit$Z[[1]]) + R
Vi <- solve(V)
X <- fit$X
P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi
all.equal(as.numeric(t(fit$y) %*% P %*% fit$y), fit$yPy) # TRUE
Ri <- solve(R)
Gi <- solve(G)
W <- cbind(fit$X, fit$Z[[1]])
C <- as(t(W) %*% Ri %*% W, "unpackedMatrix")
C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] <- C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] + Gi
all.equal(sum(log(eigen(C, only.values = TRUE)$values)), fit$logDetC) # TRUE
Ci <- solve(C)
Cire1 <- Ci[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)]
traces <- Matrix(0, nrow(sigma), ncol(sigma), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigma)) {
  for (j in i:ncol(sigma)) {
    traces[i, j] <- sum(Cire1[((i - 1) * nrow(Ki) + 1):(i * nrow(Ki)), ((j - 1) * nrow(Ki) + 1):(j * nrow(Ki))] * Ki)
  }
}
traces <- forceSymmetric(traces)
all.equal(traces, fit$traces[[1]]) # TRUE
all.equal(vars * -ED * exp(-fit$k[11, 1] * ED), fit$derivs[[1]][[11]], check.attributes = FALSE) # TRUE
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

# exfa + gk for 10 environments, fitted as dense ====
rm(list = setdiff(ls(), c("d", "K", "Ki", "ED")))
set.seed(1997)
dss <- droplevels(d[d$Man == "R" & d$Env %in% levels(d$Env)[1:10],])
dss <- dss[sample(1:nrow(dss), nrow(dss)),] # Shuffling data.
dss$GY[sample(1:nrow(dss), floor(nrow(dss) / 4))] <- NA # Introducing missing values.

fit <- kmmr(fixed = GY ~ -1 + Env,
            random = ~ re(~ exfa(Env, 1):lk(Gen, Ki), vars = ~ Env) +
                       re(~ gk(Env, ED):lk(Gen, Ki), vars = ~ Env),
            residual = ~ units,
            data = dss, sparse = FALSE, maxit = 1, debugging = TRUE)

R <- Diagonal(nrow(fit$y), fit$k[32, 1])
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), fit$logDetR) # TRUE
## Random effect 1 ====
vars1 <- outer(sqrt(c(fit$k[1:10, 1], 1)), sqrt(c(fit$k[1:10, 1], 1)))
all.equal(vars1, as.matrix(fit$vars[[1]])) # TRUE
L <- Matrix(fit$k[11:20, 1], length(levels(dss$Env)), 1)
sigma1 <- vars1 * rbind(cbind(L %*% t(L) + diag(1 - diag(L %*% t(L))), -L), cbind(-t(L), Matrix(1, 1, 1)))
sigma1 <- forceSymmetric(sigma1)
all.equal(sigma1, as(solve(fit$sigmai[[1]]), "unpackedMatrix")) # TRUE
all.equal(sum(log(eigen(sigma1, only.values = TRUE)$values)), fit$logDetvarGtotal[1]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit$logDetnvarGtotal[1]) # TRUE
G1 <- kronecker(sigma1, K)
## Random effect 2 ====
vars2 <- outer(sqrt(c(fit$k[21:30, 1])), sqrt(c(fit$k[21:30, 1])))
all.equal(vars2, as.matrix(fit$vars[[2]])) # TRUE
sigma2 <- vars2 * exp(-fit$k[31, 1] * ED)
sigma2 <- forceSymmetric(sigma2)
all.equal(sigma2, solve(fit$sigmai[[2]])) # TRUE
all.equal(sum(log(eigen(sigma2, only.values = TRUE)$values)), fit$logDetvarGtotal[2]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit$logDetnvarGtotal[2]) # TRUE
G2 <- kronecker(sigma2, K)

G <- as(bdiag(G1, G2), "unpackedMatrix")
V <- cbind(fit$Z[[1]], fit$Z[[2]]) %*% G %*% t(cbind(fit$Z[[1]], fit$Z[[2]])) + R
Vi <- solve(V)
X <- fit$X
P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi
all.equal(as.numeric(t(fit$y) %*% P %*% fit$y), fit$yPy) # TRUE
Ri <- solve(R)
Gi <- solve(G)
W <- cbind(fit$X, fit$Z[[1]], fit$Z[[2]])
C <- as(t(W) %*% Ri %*% W, "unpackedMatrix")
C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] <- C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] + Gi
all.equal(sum(log(eigen(C, only.values = TRUE)$values)), fit$logDetC) # TRUE
Ci <- solve(C)
## Random effect 1 ====
Cire1 <- Ci[(ncol(fit$X) + 1):560, (ncol(fit$X) + 1):560]
traces1 <- Matrix(0, nrow(sigma1), ncol(sigma1), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigma1)) {
  for (j in i:ncol(sigma1)) {
    traces1[i, j] <- sum(Cire1[((i - 1) * nrow(Ki) + 1):(i * nrow(Ki)), ((j - 1) * nrow(Ki) + 1):(j * nrow(Ki))] * Ki)
  }
}
traces1 <- forceSymmetric(traces1)
all.equal(traces1, fit$traces[[1]]) # TRUE
## Random effect 2 ====
Cire2 <- Ci[561:nrow(Ci), 561:ncol(Ci)]
traces2 <- Matrix(0, nrow(sigma2), ncol(sigma2), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigma2)) {
  for (j in i:ncol(sigma2)) {
    traces2[i, j] <- sum(Cire2[((i - 1) * nrow(Ki) + 1):(i * nrow(Ki)), ((j - 1) * nrow(Ki) + 1):(j * nrow(Ki))] * Ki)
  }
}
traces2 <- forceSymmetric(traces2)
all.equal(traces2, fit$traces[[2]]) # TRUE

Wu.check <- Matrix(0, nrow(fit$y), ncol(fit$Wu), sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit$derivs[[1]])) {
  Wu.check[, param] <- fit$Z[[1]] %*% kronecker(fit$derivs[[1]][[i]], K) %*% solve(G1) %*% fit$u[1:550, 1, drop = FALSE]
  param <- param + 1
}
for (i in 1:length(fit$derivs[[2]])) {
  Wu.check[, param] <- fit$Z[[2]] %*% kronecker(fit$derivs[[2]][[i]], K) %*% solve(G2) %*% fit$u[551:1050, 1, drop = FALSE]
  param <- param + 1
}
for (i in 1:length(fit$derivs[[3]])) {
  Wu.check[, param] <- fit$derivs[[3]][[i]] %*% Ri %*% fit$e
  param <- param + 1
}
all.equal(Wu.check, fit$Wu) # TRUE
scores.check <- Matrix(0, nrow(fit$scores), 1, sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit$derivs[[1]])) {
  s1 <- sum(diag(kronecker(fit$derivs[[1]][[i]], K) %*% solve(G1)))
  s2 <- sum(diag(solve(G1) %*% (kronecker(fit$derivs[[1]][[i]], K) %*% solve(G1)) %*% Cire1))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
for (i in 1:length(fit$derivs[[2]])) {
  s1 <- sum(diag(kronecker(fit$derivs[[2]][[i]], K) %*% solve(G2)))
  s2 <- sum(diag(solve(G2) %*% (kronecker(fit$derivs[[2]][[i]], K) %*% solve(G2)) %*% Cire2))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
for (i in 1:length(fit$derivs[[3]])) {
  s1 <- sum(diag(fit$derivs[[3]][[i]] %*% Ri))
  s2 <- sum(diag((t(W) %*% Ri %*% fit$derivs[[3]][[i]] %*% Ri %*% W) %*% Ci))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
all.equal(scores.check, fit$scores) # TRUE

# exfa + gk for 10 environments, fitted as sparse ====
rm(list = setdiff(ls(), c("d", "K", "Ki", "ED")))
set.seed(1997)
dss <- droplevels(d[d$Man == "R" & d$Env %in% levels(d$Env)[1:10],])
dss <- dss[sample(1:nrow(dss), nrow(dss)),] # Shuffling data.
dss$GY[sample(1:nrow(dss), floor(nrow(dss) / 4))] <- NA # Introducing missing values.

fit <- kmmr(fixed = GY ~ -1 + Env,
            random = ~ re(~ exfa(Env, 1):lk(Gen, Ki), vars = ~ Env) +
              re(~ gk(Env, ED):lk(Gen, Ki), vars = ~ Env),
            residual = ~ units,
            data = dss, sparse = TRUE, maxit = 1, debugging = TRUE)

R <- Diagonal(nrow(fit$y), fit$k[32, 1])
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), fit$logDetR) # TRUE
## Random effect 1 ====
vars1 <- outer(sqrt(c(fit$k[1:10, 1], 1)), sqrt(c(fit$k[1:10, 1], 1)))
all.equal(vars1, as.matrix(fit$vars[[1]])) # TRUE
L <- Matrix(fit$k[11:20, 1], length(levels(dss$Env)), 1)
sigma1 <- vars1 * rbind(cbind(L %*% t(L) + diag(1 - diag(L %*% t(L))), -L), cbind(-t(L), Matrix(1, 1, 1)))
sigma1 <- forceSymmetric(sigma1)
all.equal(sigma1, as(solve(fit$sigmai[[1]]), "unpackedMatrix")) # TRUE
all.equal(sum(log(eigen(sigma1, only.values = TRUE)$values)), fit$logDetvarGtotal[1]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit$logDetnvarGtotal[1]) # TRUE
G1 <- kronecker(sigma1, K)
## Random effect 2 ====
vars2 <- outer(sqrt(c(fit$k[21:30, 1])), sqrt(c(fit$k[21:30, 1])))
all.equal(vars2, as.matrix(fit$vars[[2]])) # TRUE
sigma2 <- vars2 * exp(-fit$k[31, 1] * ED)
sigma2 <- forceSymmetric(sigma2)
all.equal(sigma2, solve(fit$sigmai[[2]])) # TRUE
all.equal(sum(log(eigen(sigma2, only.values = TRUE)$values)), fit$logDetvarGtotal[2]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit$logDetnvarGtotal[2]) # TRUE
G2 <- kronecker(sigma2, K)

G <- as(bdiag(G1, G2), "unpackedMatrix")
V <- cbind(fit$Z[[1]], fit$Z[[2]]) %*% G %*% t(cbind(fit$Z[[1]], fit$Z[[2]])) + R
Vi <- solve(V)
X <- fit$X
P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi
all.equal(as.numeric(t(fit$y) %*% P %*% fit$y), fit$yPy) # TRUE
Ri <- solve(R)
Gi <- solve(G)
W <- cbind(fit$X, fit$Z[[1]], fit$Z[[2]])
C <- as(t(W) %*% Ri %*% W, "unpackedMatrix")
C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] <- C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] + Gi
all.equal(sum(log(eigen(C, only.values = TRUE)$values)), fit$logDetC) # TRUE
Ci <- solve(C)
## Random effect 1 ====
Cire1 <- Ci[(ncol(fit$X) + 1):560, (ncol(fit$X) + 1):560]
traces1 <- Matrix(0, nrow(sigma1), ncol(sigma1), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigma1)) {
  for (j in i:ncol(sigma1)) {
    traces1[i, j] <- sum(Cire1[((i - 1) * nrow(Ki) + 1):(i * nrow(Ki)), ((j - 1) * nrow(Ki) + 1):(j * nrow(Ki))] * Ki)
  }
}
traces1 <- forceSymmetric(traces1)
all.equal(traces1[fit$traces[[1]] != 0], fit$traces[[1]][fit$traces[[1]] != 0]) # TRUE
## Random effect 2 ====
Cire2 <- Ci[561:nrow(Ci), 561:ncol(Ci)]
traces2 <- Matrix(0, nrow(sigma2), ncol(sigma2), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigma2)) {
  for (j in i:ncol(sigma2)) {
    traces2[i, j] <- sum(Cire2[((i - 1) * nrow(Ki) + 1):(i * nrow(Ki)), ((j - 1) * nrow(Ki) + 1):(j * nrow(Ki))] * Ki)
  }
}
traces2 <- forceSymmetric(traces2)
all.equal(traces2, fit$traces[[2]]) # TRUE

Wu.check <- Matrix(0, nrow(fit$y), ncol(fit$Wu), sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit$derivs[[1]])) {
  Wu.check[, param] <- fit$Z[[1]] %*% kronecker(fit$derivs[[1]][[i]], K) %*% solve(G1) %*% fit$u[1:550, 1, drop = FALSE]
  param <- param + 1
}
for (i in 1:length(fit$derivs[[2]])) {
  Wu.check[, param] <- fit$Z[[2]] %*% kronecker(fit$derivs[[2]][[i]], K) %*% solve(G2) %*% fit$u[551:1050, 1, drop = FALSE]
  param <- param + 1
}
for (i in 1:length(fit$derivs[[3]])) {
  Wu.check[, param] <- fit$derivs[[3]][[i]] %*% Ri %*% fit$e
  param <- param + 1
}
all.equal(Wu.check, fit$Wu) # TRUE
scores.check <- Matrix(0, nrow(fit$scores), 1, sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit$derivs[[1]])) {
  s1 <- sum(diag(kronecker(fit$derivs[[1]][[i]], K) %*% solve(G1)))
  s2 <- sum(diag(solve(G1) %*% (kronecker(fit$derivs[[1]][[i]], K) %*% solve(G1)) %*% Cire1))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
for (i in 1:length(fit$derivs[[2]])) {
  s1 <- sum(diag(kronecker(fit$derivs[[2]][[i]], K) %*% solve(G2)))
  s2 <- sum(diag(solve(G2) %*% (kronecker(fit$derivs[[2]][[i]], K) %*% solve(G2)) %*% Cire2))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
for (i in 1:length(fit$derivs[[3]])) {
  s1 <- sum(diag(fit$derivs[[3]][[i]] %*% Ri))
  s2 <- sum(diag((t(W) %*% Ri %*% fit$derivs[[3]][[i]] %*% Ri %*% W) %*% Ci))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
all.equal(scores.check, fit$scores) # TRUE

# Single-stage with corg for 2 environments, fitted as sparse ====
rm(list = ls())
d <- readRDS("~/repos/kmmr/checks/dataPlotsDROPS.rds")
K <- readRDS("~/repos/kmmr/checks/kinshipPlotsDROPS.rds")
Ki <- solve(K)
attr(Ki, "logDet") <- sum(log(eigen(K, only.values = TRUE)$values))
dss <- droplevels(d[d$Env %in% c("Gai12", "Cra12"),])

fit <- kmmr(fixed = GY ~ 1 + Env,
            random = ~ atEach("Env", ~ re(~ ar1(Row):ar1(Col), vars = ~ 1)) +
                       re(~ corg(Env):lk(Gen, Ki), vars = ~ Env),
                 residual = ~ re(~ id(Env):id(Gen):id(Rep), vars = ~ Env),
                 data = dss, sparse = TRUE, maxit = 1, debugging = TRUE)

all.equal(c(rep("Gai12", 720), rep("Cra12", 721)), as.character(dss$Env))
levels(dss$Env) # Swapped around! Cra12 first, then Gai12!
R <- Diagonal(nrow(fit$y), c(rep(fit$k[11, 1], 720), rep(fit$k[10, 1], 721))) # First 11, then 10!
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), fit$logDetR) # TRUE
## Random effect 1 ====
vars1 <- outer(sqrt(c(fit$k[1:1, 1])), sqrt(c(fit$k[1:1, 1])))
all.equal(vars1, as.matrix(fit$vars[[1]])) # TRUE
sigma1 <- as.numeric(vars1) * kronecker(fit$k[2, 1]^fit$bases[[1]][[1]], fit$k[3, 1]^fit$bases[[1]][[2]])
sigma1 <- forceSymmetric(sigma1)
all.equal(sigma1, as(solve(fit$sigmai[[1]]), "unpackedMatrix")) # TRUE
all.equal(sum(log(eigen(sigma1, only.values = TRUE)$values)), fit$logDetvarGtotal[1]) # TRUE
G1 <- sigma1
## Random effect 2 ====
vars2 <- outer(sqrt(c(fit$k[4:4, 1])), sqrt(c(fit$k[4:4, 1])))
all.equal(vars2, as.matrix(fit$vars[[2]])) # TRUE
sigma2 <- as.numeric(vars2) * kronecker(fit$k[5, 1]^fit$bases[[2]][[1]], fit$k[6, 1]^fit$bases[[2]][[2]])
sigma2 <- forceSymmetric(sigma2)
all.equal(sigma2, as(solve(fit$sigmai[[2]]), "unpackedMatrix")) # TRUE
all.equal(sum(log(eigen(sigma2, only.values = TRUE)$values)), fit$logDetvarGtotal[2]) # TRUE
G2 <- sigma2
## Random effect 3 ====
vars3 <- outer(sqrt(c(fit$k[7:8, 1])), sqrt(c(fit$k[7:8, 1])))
all.equal(vars3, as.matrix(fit$vars[[3]])) # TRUE
sigma3 <- vars3 * Matrix(c(1, fit$k[9, 1], fit$k[9, 1], 1), 2, 2)
sigma3 <- forceSymmetric(sigma3)
all.equal(sigma3, as(solve(fit$sigmai[[3]]), "unpackedMatrix")) # TRUE
all.equal(sum(log(eigen(sigma3, only.values = TRUE)$values)), fit$logDetvarGtotal[3]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit$logDetnvarGtotal[3]) # TRUE
G3 <- kronecker(sigma3, K)

G <- as(bdiag(G1, G2, G3), "unpackedMatrix")
V <- cbind(fit$Z[[1]], fit$Z[[2]], fit$Z[[3]]) %*% G %*% t(cbind(fit$Z[[1]], fit$Z[[2]], fit$Z[[3]])) + R
Vi <- solve(V)
X <- fit$X
P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi
all.equal(as.numeric(t(fit$y) %*% P %*% fit$y), fit$yPy) # TRUE
Ri <- solve(R)
Gi <- solve(G)
W <- cbind(fit$X, fit$Z[[1]], fit$Z[[2]], fit$Z[[3]])
C <- as(t(W) %*% Ri %*% W, "unpackedMatrix")
C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] <- C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] + Gi
all.equal(sum(log(eigen(C, only.values = TRUE)$values)), fit$logDetC) # TRUE
Ci <- solve(C)
## Random effect 1 ====
Cire1 <- Matrix(as.matrix(Ci)[(ncol(fit$X) + 1):758, (ncol(fit$X) + 1):758])
## Random effect 2 ====
Cire2 <- Matrix(as.matrix(Ci)[759:1514, 759:1514])
## Random effect 3 ====
Cire3 <- Matrix(as.matrix(Ci)[1515:2004, 1515:2004])
traces3 <- Matrix(0, nrow(sigma3), ncol(sigma3), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigma3)) {
  for (j in i:ncol(sigma3)) {
    traces3[i, j] <- sum(Cire3[((i - 1) * nrow(Ki) + 1):(i * nrow(Ki)), ((j - 1) * nrow(Ki) + 1):(j * nrow(Ki))] * Ki)
  }
}
traces3 <- forceSymmetric(traces3)
all.equal(traces3, fit$traces[[3]]) # TRUE

Wu.check <- Matrix(0, nrow(fit$y), ncol(fit$Wu), sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit$derivs[[1]])) {
  Wu.check[, param] <- fit$Z[[1]] %*% fit$derivs[[1]][[i]] %*% solve(G1) %*% fit$u[1:756, 1, drop = FALSE]
  param <- param + 1
}
for (i in 1:length(fit$derivs[[2]])) {
  Wu.check[, param] <- fit$Z[[2]] %*% fit$derivs[[2]][[i]] %*% solve(G2) %*% fit$u[757:1512, 1, drop = FALSE]
  param <- param + 1
}
for (i in 1:length(fit$derivs[[3]])) {
  Wu.check[, param] <- fit$Z[[3]] %*% kronecker(fit$derivs[[3]][[i]], K) %*% solve(G3) %*% fit$u[1513:2002, 1, drop = FALSE]
  param <- param + 1
}
for (i in 1:length(fit$derivs[[4]])) {
  Wu.check[, param] <- fit$derivs[[4]][[i]] %*% Ri %*% fit$e
  param <- param + 1
}
all.equal(Wu.check, fit$Wu) # TRUE
scores.check <- Matrix(0, nrow(fit$scores), 1, sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit$derivs[[1]])) {
  s1 <- sum(diag(fit$derivs[[1]][[i]] %*% solve(G1)))
  s2 <- sum(diag(solve(G1) %*% fit$derivs[[1]][[i]] %*% solve(G1) %*% Cire1))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
for (i in 1:length(fit$derivs[[2]])) {
  s1 <- sum(diag(fit$derivs[[2]][[i]] %*% solve(G2)))
  s2 <- sum(diag(solve(G2) %*% fit$derivs[[2]][[i]] %*% solve(G2) %*% Cire2))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
for (i in 1:length(fit$derivs[[3]])) {
  s1 <- sum(diag(kronecker(fit$derivs[[3]][[i]], K) %*% solve(G3)))
  s2 <- sum(diag(solve(G3) %*% (kronecker(fit$derivs[[3]][[i]], K) %*% solve(G3)) %*% Cire3))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
for (i in 1:length(fit$derivs[[4]])) {
  s1 <- sum(diag(fit$derivs[[4]][[i]] %*% Ri))
  s2 <- sum(diag((t(W) %*% Ri %*% fit$derivs[[4]][[i]] %*% Ri %*% W) %*% Ci))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
all.equal(scores.check, fit$scores) # TRUE

# Single-stage with corg for 2 environments, fitted as dense ====
rm(list = ls())
d <- readRDS("~/repos/kmmr/checks/dataPlotsDROPS.rds")
K <- readRDS("~/repos/kmmr/checks/kinshipPlotsDROPS.rds")
Ki <- solve(K)
attr(Ki, "logDet") <- sum(log(eigen(K, only.values = TRUE)$values))
dss <- droplevels(d[d$Env %in% c("Gai12", "Cra12"),])

fit <- kmmr(fixed = GY ~ 1 + Env,
            random = ~ atEach("Env", ~ re(~ ar1(Row):ar1(Col), vars = ~ 1)) +
              re(~ corg(Env):lk(Gen, Ki), vars = ~ Env),
            residual = ~ re(~ id(Env):id(Gen):id(Rep), vars = ~ Env),
            data = dss, sparse = FALSE, maxit = 1, debugging = TRUE)

all.equal(c(rep("Gai12", 720), rep("Cra12", 721)), as.character(dss$Env))
levels(dss$Env) # Swapped around! Cra12 first, then Gai12!
R <- Diagonal(nrow(fit$y), c(rep(fit$k[11, 1], 720), rep(fit$k[10, 1], 721))) # First 11, then 10!
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), fit$logDetR) # TRUE
## Random effect 1 ====
vars1 <- outer(sqrt(c(fit$k[1:1, 1])), sqrt(c(fit$k[1:1, 1])))
all.equal(vars1, as.matrix(fit$vars[[1]])) # TRUE
sigma1 <- as.numeric(vars1) * kronecker(fit$k[2, 1]^fit$bases[[1]][[1]], fit$k[3, 1]^fit$bases[[1]][[2]])
sigma1 <- forceSymmetric(sigma1)
all.equal(sigma1, as(solve(fit$sigmai[[1]]), "unpackedMatrix")) # TRUE
all.equal(sum(log(eigen(sigma1, only.values = TRUE)$values)), fit$logDetvarGtotal[1]) # TRUE
G1 <- sigma1
## Random effect 2 ====
vars2 <- outer(sqrt(c(fit$k[4:4, 1])), sqrt(c(fit$k[4:4, 1])))
all.equal(vars2, as.matrix(fit$vars[[2]])) # TRUE
sigma2 <- as.numeric(vars2) * kronecker(fit$k[5, 1]^fit$bases[[2]][[1]], fit$k[6, 1]^fit$bases[[2]][[2]])
sigma2 <- forceSymmetric(sigma2)
all.equal(sigma2, as(solve(fit$sigmai[[2]]), "unpackedMatrix")) # TRUE
all.equal(sum(log(eigen(sigma2, only.values = TRUE)$values)), fit$logDetvarGtotal[2]) # TRUE
G2 <- sigma2
## Random effect 3 ====
vars3 <- outer(sqrt(c(fit$k[7:8, 1])), sqrt(c(fit$k[7:8, 1])))
all.equal(vars3, as.matrix(fit$vars[[3]])) # TRUE
sigma3 <- vars3 * Matrix(c(1, fit$k[9, 1], fit$k[9, 1], 1), 2, 2)
sigma3 <- forceSymmetric(sigma3)
all.equal(sigma3, as(solve(fit$sigmai[[3]]), "unpackedMatrix")) # TRUE
all.equal(sum(log(eigen(sigma3, only.values = TRUE)$values)), fit$logDetvarGtotal[3]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit$logDetnvarGtotal[3]) # TRUE
G3 <- kronecker(sigma3, K)

G <- as(bdiag(G1, G2, G3), "unpackedMatrix")
V <- cbind(fit$Z[[1]], fit$Z[[2]], fit$Z[[3]]) %*% G %*% t(cbind(fit$Z[[1]], fit$Z[[2]], fit$Z[[3]])) + R
Vi <- solve(V)
X <- fit$X
P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi
all.equal(as.numeric(t(fit$y) %*% P %*% fit$y), fit$yPy) # TRUE
Ri <- solve(R)
Gi <- solve(G)
W <- cbind(fit$X, fit$Z[[1]], fit$Z[[2]], fit$Z[[3]])
C <- as(t(W) %*% Ri %*% W, "unpackedMatrix")
C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] <- C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] + Gi
all.equal(sum(log(eigen(C, only.values = TRUE)$values)), fit$logDetC) # TRUE
Ci <- solve(C)
## Random effect 1 ====
Cire1 <- Matrix(as.matrix(Ci)[(ncol(fit$X) + 1):758, (ncol(fit$X) + 1):758])
## Random effect 2 ====
Cire2 <- Matrix(as.matrix(Ci)[759:1514, 759:1514])
## Random effect 3 ====
Cire3 <- Matrix(as.matrix(Ci)[1515:2004, 1515:2004])
traces3 <- Matrix(0, nrow(sigma3), ncol(sigma3), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigma3)) {
  for (j in i:ncol(sigma3)) {
    traces3[i, j] <- sum(Cire3[((i - 1) * nrow(Ki) + 1):(i * nrow(Ki)), ((j - 1) * nrow(Ki) + 1):(j * nrow(Ki))] * Ki)
  }
}
traces3 <- forceSymmetric(traces3)
all.equal(traces3, fit$traces[[3]]) # TRUE

Wu.check <- Matrix(0, nrow(fit$y), ncol(fit$Wu), sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit$derivs[[1]])) {
  Wu.check[, param] <- fit$Z[[1]] %*% fit$derivs[[1]][[i]] %*% solve(G1) %*% fit$u[1:756, 1, drop = FALSE]
  param <- param + 1
}
for (i in 1:length(fit$derivs[[2]])) {
  Wu.check[, param] <- fit$Z[[2]] %*% fit$derivs[[2]][[i]] %*% solve(G2) %*% fit$u[757:1512, 1, drop = FALSE]
  param <- param + 1
}
for (i in 1:length(fit$derivs[[3]])) {
  Wu.check[, param] <- fit$Z[[3]] %*% kronecker(fit$derivs[[3]][[i]], K) %*% solve(G3) %*% fit$u[1513:2002, 1, drop = FALSE]
  param <- param + 1
}
for (i in 1:length(fit$derivs[[4]])) {
  Wu.check[, param] <- fit$derivs[[4]][[i]] %*% Ri %*% fit$e
  param <- param + 1
}
all.equal(Wu.check, fit$Wu) # TRUE
scores.check <- Matrix(0, nrow(fit$scores), 1, sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit$derivs[[1]])) {
  s1 <- sum(diag(fit$derivs[[1]][[i]] %*% solve(G1)))
  s2 <- sum(diag(solve(G1) %*% fit$derivs[[1]][[i]] %*% solve(G1) %*% Cire1))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
for (i in 1:length(fit$derivs[[2]])) {
  s1 <- sum(diag(fit$derivs[[2]][[i]] %*% solve(G2)))
  s2 <- sum(diag(solve(G2) %*% fit$derivs[[2]][[i]] %*% solve(G2) %*% Cire2))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
for (i in 1:length(fit$derivs[[3]])) {
  s1 <- sum(diag(kronecker(fit$derivs[[3]][[i]], K) %*% solve(G3)))
  s2 <- sum(diag(solve(G3) %*% (kronecker(fit$derivs[[3]][[i]], K) %*% solve(G3)) %*% Cire3))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
for (i in 1:length(fit$derivs[[4]])) {
  s1 <- sum(diag(fit$derivs[[4]][[i]] %*% Ri))
  s2 <- sum(diag((t(W) %*% Ri %*% fit$derivs[[4]][[i]] %*% Ri %*% W) %*% Ci))
  s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
  scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  param <- param + 1
}
all.equal(scores.check, fit$scores) # TRUE

# Multi-trait with 2 x corg for 2 traits ====
rm(list = ls())
d <- kmmr::dataDROPS
K <- kmmr::kinshipDROPS
Ki <- solve(K)
attr(Ki, "logDet") <- sum(log(eigen(K, only.values = TRUE)$values))
set.seed(1997)
dss <- droplevels(d[d$Man == "R" & d$Env %in% levels(d$Env)[1:2],])
dss <- dss[sample(1:nrow(dss), nrow(dss)),] # Shuffling data.
dss$GY[sample(1:nrow(dss), floor(nrow(dss) / 4))] <- NA # Introducing missing values.

fit <- kmmr(fixed = GY ~ -1 + Env,
            random = ~ re(~ corg(Env):lk(Gen, Ki), vars = ~ Env),
            residual = ~ re(~ corg(Env):id(Gen), vars = ~ Env),
            data = dss, sparse = FALSE, maxit = 1, debugging = TRUE)

R <- kronecker(fit$vars[[2]] * Matrix(c(1, fit$k[6, 1], fit$k[6, 1], 1), 2, 2),
               fit$bases[[2]][[2]])
R <- fit$Z[[2]] %*% R %*% t(fit$Z[[2]])
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), as.numeric(fit$logDetR)) # TRUE
vars <- outer(sqrt(fit$k[1:2, 1]), sqrt(fit$k[1:2, 1]))
all.equal(vars, as.matrix(fit$vars[[1]])) # TRUE
sigma <- vars * Matrix(c(1, fit$k[3, 1], fit$k[3, 1], 1), 2, 2)
all.equal(sigma, solve(fit$sigmai[[1]])) # TRUE
all.equal(sum(log(eigen(sigma, only.values = TRUE)$values)), fit$logDetvarGtotal[1]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit$logDetnvarGtotal[1]) # TRUE
G <- kronecker(sigma, K)
all.equal(sum(log(eigen(G, only.values = TRUE)$values)), fit$logDetG) # TRUE
V <- fit$Z[[1]] %*% G %*% t(fit$Z[[1]]) + R
Vi <- solve(V)
X <- fit$X
P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi
all.equal(as.numeric(t(fit$y) %*% P %*% fit$y), fit$yPy) # TRUE
Ri <- solve(R)
Gi <- solve(G)
W <- cbind(fit$X, fit$Z[[1]])
C <- as(t(W) %*% Ri %*% W, "unpackedMatrix")
C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] <- C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] + Gi
all.equal(sum(log(eigen(C, only.values = TRUE)$values)), fit$logDetC) # TRUE
Ci <- solve(C)
Cire1 <- Ci[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)]
traces <- Matrix(0, nrow(sigma), ncol(sigma), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigma)) {
  for (j in i:ncol(sigma)) {
    traces[i, j] <- sum(Cire1[((i - 1) * nrow(Ki) + 1):(i * nrow(Ki)), ((j - 1) * nrow(Ki) + 1):(j * nrow(Ki))] * Ki)
  }
}
traces <- forceSymmetric(traces)
all.equal(traces, fit$traces[[1]]) # TRUE
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

# exfa for 10 environments with weights and intercept ====
rm(list = setdiff(ls(), c("d", "K", "Ki", "fit")))
set.seed(1997)
dss <- droplevels(d[d$Man == "R" & d$Env %in% levels(d$Env)[1:10],])
dss$weight.kmmr <- 1 / dss$weight
dss <- dss[sample(1:nrow(dss), nrow(dss)),] # Shuffling data.
dss$GY[sample(1:nrow(dss), floor(nrow(dss) / 4))] <- NA # Introducing missing values.

fit <- kmmr(fixed = GY ~ 1 + Env,
            random = ~ re(~ exfa(Env, 2):lk(Gen, Ki), vars = ~ Env),
            residual = ~ units,
            weights = "weight.kmmr",
            data = dss, sparse = TRUE, maxit = 1, debugging = TRUE, denseThreshold = 0.2,
            getStartingValues = TRUE)

old.var <- fit$Init[31]
new.var <- mean(na.omit(dss)$weight.kmmr)
inc <- new.var / old.var
fit$Init[31] <- 1
fit$Fixed[31] <- TRUE
fit$Init[1:10] <- fit$Init[1:10] * inc

fit <- kmmr(fixed = GY ~ 1 + Env,
            random = ~ re(~ exfa(Env, 2):lk(Gen, Ki), vars = ~ Env),
            residual = ~ units,
            weights = "weight.kmmr",
            init = fit$Init, isFixed = fit$Fixed,
            data = dss, sparse = TRUE, maxit = 1, debugging = TRUE, denseThreshold = 0.2)

R <- Diagonal(nrow(fit$y), na.omit(dss)$weight.kmmr)
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), fit$logDetR) # TRUE
vars <- outer(sqrt(c(fit$k[1:10, 1], 1, 1)), sqrt(c(fit$k[1:10, 1], 1, 1)))
all.equal(vars, as.matrix(fit$vars[[1]])) # TRUE
L <- Matrix(fit$k[11:30, 1], length(levels(dss$Env)), 2)
sigma <- vars * rbind(cbind(L %*% t(L) + diag(1 - diag(L %*% t(L))), -L), cbind(-t(L), Diagonal(2, 1)))
sigma <- forceSymmetric(sigma)
all.equal(sigma, as(solve(fit$sigmai[[1]]), "unpackedMatrix")) # TRUE
all.equal(sum(log(eigen(sigma, only.values = TRUE)$values)), fit$logDetvarGtotal[1]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit$logDetnvarGtotal[1]) # TRUE
G <- kronecker(sigma, K)
all.equal(sum(log(eigen(G, only.values = TRUE)$values)), fit$logDetG) # TRUE
V <- fit$Z[[1]] %*% G %*% t(fit$Z[[1]]) + R
Vi <- solve(V)
X <- fit$X
P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi
all.equal(as.numeric(t(fit$y) %*% P %*% fit$y), fit$yPy) # TRUE
Ri <- solve(R)
Gi <- solve(G)
W <- cbind(fit$X, fit$Z[[1]])
C <- as(t(W) %*% Ri %*% W, "unpackedMatrix")
C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] <- C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] + Gi
all.equal(sum(log(eigen(C, only.values = TRUE)$values)), as.numeric(fit$logDetC)) # TRUE
Ci <- solve(C)
Cire1 <- Ci[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)]
traces <- Matrix(0, nrow(sigma), ncol(sigma), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigma)) {
  for (j in i:ncol(sigma)) {
    traces[i, j] <- sum(Cire1[((i - 1) * nrow(Ki) + 1):(i * nrow(Ki)), ((j - 1) * nrow(Ki) + 1):(j * nrow(Ki))] * Ki)
  }
}
traces <- forceSymmetric(traces)
all.equal(traces[fit$traces[[1]] != 0], fit$traces[[1]][fit$traces[[1]] != 0]) # TRUE
Wu.check <- Matrix(0, nrow(fit$y), ncol(fit$Wu), sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit$derivs[[1]])) {
  if (!is.null(fit$derivs[[1]][[i]])) {
    Wu.check[, param] <- fit$Z[[1]] %*% kronecker(fit$derivs[[1]][[i]], K) %*% Gi %*% fit$u
  }
  param <- param + 1
}
for (i in 1:length(fit$derivs[[2]])) {
  if (!is.null(fit$derivs[[2]][[i]])) {
    Wu.check[, param] <- fit$derivs[[2]][[i]] %*% Ri %*% fit$e
  }
  param <- param + 1
}
all.equal(Wu.check, fit$Wu) # TRUE
scores.check <- Matrix(0, nrow(fit$scores), 1, sparse = FALSE, doDiag = FALSE)
param <- 1
for (i in 1:length(fit$derivs[[1]])) {
  if (!is.null(fit$derivs[[1]][[i]])) {
    s1 <- sum(diag(kronecker(fit$derivs[[1]][[i]], K) %*% Gi))
    s2 <- sum(diag(Gi %*% (kronecker(fit$derivs[[1]][[i]], K) %*% Gi) %*% Cire1))
    s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
    scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  }
  param <- param + 1
}
for (i in 1:length(fit$derivs[[2]])) {
  if (!is.null(fit$derivs[[2]][[i]])) {
    s1 <- sum(diag(fit$derivs[[2]][[i]] %*% Ri))
    s2 <- sum(diag((t(W) %*% Ri %*% fit$derivs[[2]][[i]] %*% Ri %*% W) %*% Ci))
    s3 <- as.numeric(t(fit$y) %*% P %*% fit$Wu[, param])
    scores.check[param, 1] <- -0.5 * (s1 - s2 - s3)
  }
  param <- param + 1
}
all.equal(scores.check, fit$scores) # TRUE

# gk + fa stacked for 10 environments ====
rm(list = setdiff(ls(), c("d", "K", "Ki")))
ED <- kmmr::distmatDROPS
set.seed(1997)
dss <- droplevels(d[d$Man == "R" & d$Env %in% levels(d$Env)[1:10],])
dss <- dss[sample(1:nrow(dss), nrow(dss)),] # Shuffling data.
dss$GY[sample(1:nrow(dss), floor(nrow(dss) / 4))] <- NA # Introducing missing values.
ED <- ED[levels(dss$Env), levels(dss$Env)]

fit <- kmmr(fixed = GY ~ -1 + Env,
            random = ~ re(~ stack(~ gk(Env, ED) + fa(Env, 1)):lk(Gen, Ki), vars = ~ Env),
            residual = ~ units,
            data = dss, sparse = FALSE, maxit = 2, debugging = TRUE,
            stepSize = c(0.25, 0.5)^2)

fit2 <- kmmr(fixed = GY ~ -1 + Env,
             random = ~ re(~ gk(Env, ED):lk(Gen, Ki), vars = ~ Env) +
                        re(~ fa(Env, 1):lk(Gen, Ki), vars = ~ Env),
             residual = ~ units,
             data = dss, sparse = FALSE, maxit = 2, debugging = TRUE,
             stepSize = c(0.25, 0.5)^2)

R <- Diagonal(nrow(fit$y), fit$k[32, 2])
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), fit$logDetR) # TRUE
all.equal(sum(log(eigen(R, only.values = TRUE)$values)), fit2$logDetR) # TRUE
vars1 <- outer(sqrt(fit$k[1:10, 2]), sqrt(fit$k[1:10, 2]))
vars2 <- outer(sqrt(fit$k[11:20, 2]), sqrt(fit$k[11:20, 2]))
all.equal(vars1, as.matrix(fit$vars[[1]][[1]])) # TRUE
all.equal(vars2, as.matrix(fit$vars[[1]][[2]])) # TRUE
L <- Matrix(fit$k[22:31, 2], 10, 1)
sigma1 <- forceSymmetric(vars1 * exp(-fit$k[21, 2] * ED))
sigma2 <- forceSymmetric(vars2 * (L %*% t(L) + diag(1 - diag(L %*% t(L)))))
sigma <- sigma1 + sigma2
all.equal(as.matrix(sigma), as.matrix(solve(fit$sigmai[[1]]))) # TRUE
all.equal(sum(log(eigen(sigma, only.values = TRUE)$values)), fit$logDetvarGtotal[1]) # TRUE
all.equal(sum(log(eigen(K, only.values = TRUE)$values)), fit$logDetnvarGtotal[1]) # TRUE
G <- kronecker(sigma, K)
all.equal(sum(log(eigen(G, only.values = TRUE)$values)), fit$logDetG) # TRUE
V <- fit$Z[[1]] %*% G %*% t(fit$Z[[1]]) + R
Vi <- solve(V)
X <- fit$X
P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi
all.equal(as.numeric(t(fit$y) %*% P %*% fit$y), fit$yPy) # TRUE
all.equal(as.numeric(t(fit$y) %*% P %*% fit$y), fit2$yPy) # TRUE
Ri <- solve(R)
Gi <- solve(G)
W <- cbind(fit$X, fit$Z[[1]])
C <- as(t(W) %*% Ri %*% W, "unpackedMatrix")
C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] <- C[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)] + Gi
all.equal(sum(log(eigen(C, only.values = TRUE)$values)), fit$logDetC) # TRUE
Ci <- solve(C)
Cire1 <- Ci[(ncol(fit$X) + 1):nrow(C), (ncol(fit$X) + 1):ncol(C)]
traces <- Matrix(0, nrow(sigma), ncol(sigma), sparse = FALSE, doDiag = FALSE)
for (i in 1:nrow(sigma)) {
  for (j in i:ncol(sigma)) {
    traces[i, j] <- sum(Cire1[((i - 1) * nrow(Ki) + 1):(i * nrow(Ki)), ((j - 1) * nrow(Ki) + 1):(j * nrow(Ki))] * Ki)
  }
}
traces <- forceSymmetric(traces)
all.equal(traces, fit$traces[[1]]) # TRUE
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
all.equal(Wu.check, fit2$Wu[, c(1:10, 12:21, 11, 22:32)]) # TRUE
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
all.equal(scores.check, fit2$scores[c(1:10, 12:21, 11, 22:32),, drop = FALSE]) # TRUE





