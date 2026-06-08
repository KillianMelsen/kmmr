library(asreml)
asreml.options(score = TRUE, wvr = TRUE)
K <- readRDS("dev/exampledata_K.rds")
Ki <- solve(K)
d <- readRDS("dev/exampledata.rds")
d <- droplevels(d[d$Man  == "R" & d$Env == "Cam12",])
fit <- asreml(fixed = GY ~ 1,
              random = ~ vm(Gen, K),
              residual = ~ idv(units),
              data = d)

k <- matrix(NA, 2, 20)
k[, 1] <- var(d$GY) / 20

X <- matrix(1, nrow(d), 1)
Z <- diag(nrow(d))
y <- matrix(d$GY)

# Option 1: Extended (y) formulation ====
# Option 1 uses the extended formulation where the coefficient matrix C is extended
# with the response variable. I suspect this will complicate the use of off-the-shelf
# sparse solvers if C is sparse.
Ri.1 <- diag(1 / k[2, 1], nrow(y))
W.1 <- cbind(X, Z)
Wy.1 <- cbind(W.1, y)
M.1 <- t(Wy.1) %*% Ri.1 %*% Wy.1
sigmai.1 <- matrix(1 / k[1, 1], 1, 1)
Gi.1 <- kronecker(sigmai.1, Ki)
M.1[2:247, 2:247] <- M.1[2:247, 2:247] + Gi.1
Mchol.1 <- chol(M.1)
yPy.1 <- Mchol.1[nrow(Mchol.1), ncol(Mchol.1)]^2
Mchol_XZ.1 <- Mchol.1[1:(nrow(Mchol.1) - 1), 1:(ncol(Mchol.1) - 1)]
My.1 <- Mchol.1[1:(nrow(Mchol.1) - 1), ncol(Mchol.1), drop = FALSE]
logDetC.1 <- 2 * sum(log(diag(Mchol_XZ.1)))
logDetG.1 <- -1 * sum(log(eigen(Gi.1)$values))
logDetR.1 <- -1 * sum(log(eigen(Ri.1)$values))
LL.1 <- -0.5 * (logDetC.1 + logDetR.1 + logDetG.1 + yPy.1)
bu.1 <- backsolve(Mchol_XZ.1, My.1)
b.1 <- bu.1[1, 1, drop = FALSE]
u.1 <- bu.1[2:nrow(bu.1), 1, drop = FALSE]
umat_sigmai.1 <- u.1 %*% sigmai.1
e.1 <- y - W.1 %*% bu.1
Wu.1 <- cbind(Z %*% umat_sigmai.1 %*% matrix(1, 1, 1),
              1 / k[2, 1] * e.1)
XWjxZWj.1 <- t(W.1) %*% Wu.1 * (1 / k[2, 1])
WiWj.1 <- t(Wu.1) %*% Wu.1 * (1 / k[2, 1])
buWu.1 <- backsolve(Mchol_XZ.1, forwardsolve(t(Mchol_XZ.1), XWjxZWj.1))
AI.1 <- (WiWj.1 - (t(buWu.1) %*% XWjxZWj.1)) / 2
AIi.1 <- solve(AI.1)
Ci.1 <- chol2inv(Mchol_XZ.1)

# Option 2: Split factorization ====
# This should be easier to use with off-the-shelf sparse solvers because we don't
# need to restrict any reordering due to y being in C.
Ri.2 <- diag(1 / k[2, 1], nrow(y))
W.2 <- cbind(X, Z)
M.2 <- t(W.2) %*% Ri.2 %*% W.2
cy.2 <- t(W.2) %*% Ri.2 %*% y
sigmai.2 <- matrix(1 / k[1, 1], 1, 1)
Gi.2 <- kronecker(sigmai.2, Ki)
M.2[2:247, 2:247] <- M.2[2:247, 2:247] + Gi.2
Mchol.2 <- chol(M.2)
z.2 <- forwardsolve(t(Mchol.2), cy.2)
yPy.2 <- t(y) %*% Ri.2 %*% y - t(z.2) %*% z.2
logDetC.2 <- 2 * sum(log(diag(Mchol.2)))
logDetG.2 <- -1 * sum(log(eigen(Gi.2)$values))
logDetR.2 <- -1 * sum(log(eigen(Ri.2)$values))
LL.2 <- -0.5 * (logDetC.2 + logDetR.2 + logDetG.2 + yPy.2)
bu.2 <- backsolve(Mchol.2, z.2)
b.2 <- bu.2[1, 1, drop = FALSE]
u.2 <- bu.2[2:nrow(bu.2), 1, drop = FALSE]
umat_sigmai.2 <- u.2 %*% sigmai.2
e.2 <- y - W.2 %*% bu.2
Wu.2 <- cbind(Z %*% umat_sigmai.2 %*% matrix(1, 1, 1),
              1 / k[2, 1] * e.2)
XWjxZWj.2 <- t(W.2) %*% Wu.2 * (1 / k[2, 1])
WiWj.2 <- t(Wu.2) %*% Wu.2 * (1 / k[2, 1])
buWu.2 <- backsolve(Mchol.2, forwardsolve(t(Mchol.2), XWjxZWj.2))
AI.2 <- (WiWj.2 - (t(buWu.2) %*% XWjxZWj.2)) / 2
AIi.2 <- solve(AI.2)
Ci.2 <- chol2inv(Mchol.2)














