library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = TRUE, score = TRUE, step.size = 1)

# Seed:
set.seed(1997)

# Loading example data (DROPS):
d <- readRDS("dev/exampledata.rds")[, c("Gen", "Man", "Env", "GY")]
ED <- readRDS("dev/exampledata_ED.rds")
EC <- readRDS("dev/exampledata_EC.rds")
GD <- readRDS("dev/exampledata_GD.rds")
K <- readRDS("dev/exampledata_K.rds")
Ki <- solve(K)
attr(Ki, "INVERSE") <- TRUE

# Adding a ManEnv factor for ASReml/cornfruit:
d$ManEnv <- factor(paste(as.character(d$Man), as.character(d$Env), sep = ":"),
                   levels = c(paste0("R:", levels(d$Env)),
                              paste0("W:", levels(d$Env))))

ds <- droplevels(d[d$Env %in% c("Bol12", "Kar12", "Cam12") & d$Man == "R",])
ED <- ED[levels(ds$Env), levels(ds$Env)]
EC1 <- exp(-0.05 * ED)
EC2 <- exp(-1.00 * ED)

vf <- function(order, kappa) {
  V <- kappa[1] * EC1 + kappa[2] * EC2
  dVk1 <- EC1
  dVk2 <- EC2
  return(list(V, dVk1, dVk2))
}

init <- c(0.7698575, 0.7698575)
type <- c("V", "V")
cons <- c("P", "P")
fit1 <- asreml(fixed = GY ~ -1 + Env,
               random = ~ own(Env, "vf", init, type, cons):vm(Gen, K),
               residual = ~ idv(units),
               data = ds,
               maxit = 1)

EC1K <- kronecker(EC1, K, make.dimnames = TRUE)
EC2K <- kronecker(EC2, K, make.dimnames = TRUE)
ds$EnvGen <- factor(paste(as.character(ds$Env), as.character(ds$Gen), sep = ":"),
                    paste(as.character(ds$Env), as.character(ds$Gen), sep = ":"),
                    paste(as.character(ds$Env), as.character(ds$Gen), sep = ":"))

fit2 <- asreml(fixed = GY ~ -1 + Env,
               random = ~ vm(EnvGen, EC1K) + vm(EnvGen, EC2K),
               residual = ~ idv(units),
               data = ds,
               maxit = 1)


summary(fit1)$varcomp
summary(fit2)$varcomp

# Running manually through the kmmr solver ====
# fit <- kmmr(fixed = GY ~ -1 + Env,
#             random = ~ re(~ lk(Gen, Ki)) + re(~ stack(~ gk(Env, ED) + corg(Env)):lk(Gen, Ki)),
#             residual = ~ units,
#             data = ds)

# Data preparation:
y <- matrix(ds$GY)
X <- model.matrix(~ -1 + Env, ds)
Z <- diag(nrow(ds))
nX <- ncol(X)

# Storage for LL and parameters:
LL <- numeric(10)
k <- matrix(NA, 3, 10)
k[, 1] <- as.numeric(fit1$trace[5:7, 1])

# Log|K|:
logDetnvarGtotal <- -1 * sum(log(eigen(Ki)$values))

# Non-variable part of Gi:
nvarGi <- Ki

## Iteration 1 ====
Ri <- 1 / k[nrow(k), 1] * diag(nrow(y))
W <- cbind(X, Z)
Wy <- cbind(W, y)
M0 <- t(Wy) %*% Wy
M <- M0 * (1 / k[nrow(k), 1])
sigmai <- solve(k[1, 1] * EC1 + k[2, 1] * EC2)
logDetvarGtotal <- -1 * sum(log(eigen(sigmai)$values))
Gi <- kronecker(sigmai, nvarGi)
M[4:741, 4:741] <- M[4:741, 4:741] + Gi
Mchol <- chol(M)
yPy <- Mchol[nrow(Mchol), ncol(Mchol)]^2 # R <- solve(Ri); G <- solve(Gi); Vi <- solve(G + R); P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi; all.equal(yPy, as.numeric(t(y) %*% P %*% y))
Mchol_XZ <- Mchol[1:(nrow(Mchol) - 1), 1:(ncol(Mchol) - 1)]
My <- Mchol[1:(nrow(Mchol) - 1), ncol(Mchol):ncol(Mchol), drop = FALSE]
logDetC <- 2 * sum(log(diag(Mchol_XZ))) # C <- rbind(cbind(t(X) %*% Ri %*% X, t(X) %*% Ri %*% Z), cbind(t(Z) %*% Ri %*% X, t(Z) %*% Ri %*% Z + Gi)); all.equal(C, M[1:nrow(C), 1:ncol(C)]); all.equal(logDetC, sum(log(eigen(C)$values)))
logDetG <- nrow(nvarGi) * logDetvarGtotal + nrow(sigmai) * logDetnvarGtotal # all.equal(logDetG, sum(log(eigen(kronecker(solve(sigmai), K))$values)))
logDetR <- nrow(y) * log(k[nrow(k), 1]) # all.equal(logDetR, sum(log(eigen(solve(Ri))$values)))
LL[1] <- -0.5 * (logDetC + logDetR + logDetG + yPy)
bu <- backsolve(Mchol_XZ, My)
b <- bu[1:3, 1, drop = FALSE]
u <- bu[4:nrow(bu), 1, drop = FALSE]
Wu <- matrix(NA, nrow(y), nrow(k))
umat <- matrix(u, nrow = nrow(nvarGi), ncol = nrow(sigmai))
umat_sigmai <- umat %*% sigmai
derivs <- vector("list", nrow(k) - 1)
derivs[[1]] <- EC1
derivs[[2]] <- EC2
Wu[, 1] <- Z %*% matrix(umat_sigmai %*% derivs[[1]]) # all.equal(Wu[, 1], as.numeric(Z %*% kronecker(derivs[[1]], K) %*% Gi %*% u))
Wu[, 2] <- Z %*% matrix(umat_sigmai %*% derivs[[2]]) # all.equal(Wu[, 2], as.numeric(Z %*% kronecker(derivs[[2]], K) %*% Gi %*% u))
e <- y - W %*% bu
Wu[, nrow(k)] <- (1 / k[nrow(k), 1]) * e # all.equal(Wu[, nrow(k)], as.numeric(diag(nrow(y)) %*% Ri %*% e))
XWjxZWj <- t(W) %*% Wu * (1 / k[nrow(k), 1])
WiWj <- t(Wu) %*% Wu * (1 / k[nrow(k), 1])
buWu <- backsolve(Mchol_XZ, forwardsolve(t(Mchol_XZ), XWjxZWj))
AI <- (WiWj - (t(buWu) %*% XWjxZWj)) / 2 # all.equal(AI, (t(Wu) %*% P %*% Wu) / 2)
AIi <- solve(AI)
Ci <- chol2inv(Mchol_XZ) # all.equal(Ci, solve(C), check.attributes = FALSE)
traces <- matrix(0, nrow(sigmai), ncol(sigmai))
for (i in 1:nrow(sigmai)) {
  for (j in 1:ncol(sigmai)) {
    if (j >= i) { # Only do the upper triangle
      rstart <- (nX + 1 + (i - 1) * nrow(nvarGi))
      rend <- rstart - 1 + nrow(nvarGi)
      cstart <- (nX + 1 + (j - 1) * ncol(nvarGi))
      cend <- cstart - 1 + ncol(nvarGi)
      traces[i, j] <- sum(nvarGi * Ci[rstart:rend, cstart:cend])
    }
  }
}
traces <- traces + t(traces) - diag(diag(traces))
basisScores <- -0.5 * (nrow(nvarGi) * sigmai - sigmai %*% traces %*% sigmai - t(umat_sigmai) %*% nvarGi %*% umat_sigmai)
scores <- matrix(0, nrow(k), 1)
scores[1, 1] <- sum(derivs[[1]] * basisScores) # all.equal(scores[1, 1], as.numeric(-0.5 * (sum(diag(kronecker(derivs[[1]], K) %*% Gi)) - sum(diag((Gi %*% kronecker(derivs[[1]], K) %*% Gi) %*% Ci[5:988, 5:988])) - (t(y) %*% P %*% Wu[, 1, drop = FALSE]))))
scores[2, 1] <- sum(derivs[[2]] * basisScores) # all.equal(scores[2, 1], as.numeric(-0.5 * (sum(diag(kronecker(derivs[[2]], K) %*% Gi)) - sum(diag((Gi %*% kronecker(derivs[[2]], K) %*% Gi) %*% Ci[5:988, 5:988])) - (t(y) %*% P %*% Wu[, 2, drop = FALSE]))))
WSW <- t(W) %*% diag(nrow(y)) %*% W
evari <- (1 / k[nrow(k), 1])
scores[nrow(scores), 1] <- as.numeric(-0.5 * (sum(diag(nrow(y)) * evari) - sum(evari^2 * Ci * WSW) - (evari^2 * t(e) %*% diag(nrow(y)) %*% e)))
# all.equal(scores[nrow(scores), 1], as.numeric(-0.5 * (sum(diag(diag(nrow(y)) %*% Ri)) - sum(diag((t(cbind(X, Z)) %*% Ri %*% diag(nrow(y)) %*% Ri %*% cbind(X, Z)) %*% Ci)) - (t(y) %*% P %*% Wu[, nrow(k), drop = FALSE]))))
delta <- AIi %*% scores
k[, 2] <- k[, 1] + delta

## Iteration 2 ====
Ri <- 1 / k[nrow(k), 2] * diag(nrow(y))
W <- cbind(X, Z)
Wy <- cbind(W, y)
M0 <- t(Wy) %*% Wy
M <- M0 * (1 / k[nrow(k), 2])
sigmai <- solve(k[1, 2] * EC1 + k[2, 2] * EC2)
logDetvarGtotal <- -1 * sum(log(eigen(sigmai)$values))
Gi <- kronecker(sigmai, nvarGi)
M[5:988, 5:988] <- M[5:988, 5:988] + Gi
Mchol <- chol(M)
yPy <- Mchol[nrow(Mchol), ncol(Mchol)]^2 # R <- solve(Ri); G <- solve(Gi); Vi <- solve(G + R); P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi; all.equal(yPy, as.numeric(t(y) %*% P %*% y))
Mchol_XZ <- Mchol[1:(nrow(Mchol) - 1), 1:(ncol(Mchol) - 1)]
My <- Mchol[1:(nrow(Mchol) - 1), ncol(Mchol):ncol(Mchol), drop = FALSE]
logDetC <- 2 * sum(log(diag(Mchol_XZ))) # C <- rbind(cbind(t(X) %*% Ri %*% X, t(X) %*% Ri %*% Z), cbind(t(Z) %*% Ri %*% X, t(Z) %*% Ri %*% Z + Gi)); all.equal(C, M[1:nrow(C), 1:ncol(C)]); all.equal(logDetC, sum(log(eigen(C)$values)))
logDetG <- nrow(nvarGi) * logDetvarGtotal + nrow(sigmai) * logDetnvarGtotal # all.equal(logDetG, sum(log(eigen(kronecker(solve(sigmai), K))$values)))
logDetR <- nrow(y) * log(k[nrow(k), 2]) # all.equal(logDetR, sum(log(eigen(solve(Ri))$values)))
LL[2] <- -0.5 * (logDetC + logDetR + logDetG + yPy)
bu <- backsolve(Mchol_XZ, My)
b <- bu[1:4, 1, drop = FALSE]
u <- bu[5:nrow(bu), 1, drop = FALSE]
Wu <- matrix(NA, nrow(y), nrow(k))
umat <- matrix(u, nrow = nrow(nvarGi), ncol = nrow(sigmai))
umat_sigmai <- umat %*% sigmai
derivs <- vector("list", nrow(k) - 1)
derivs[[1]] <- EC1
derivs[[2]] <- EC2
Wu[, 1] <- Z %*% matrix(umat_sigmai %*% derivs[[1]]) # all.equal(Wu[, 1], as.numeric(Z %*% kronecker(derivs[[1]], K) %*% Gi %*% u))
Wu[, 2] <- Z %*% matrix(umat_sigmai %*% derivs[[2]]) # all.equal(Wu[, 2], as.numeric(Z %*% kronecker(derivs[[2]], K) %*% Gi %*% u))
e <- y - W %*% bu
Wu[, nrow(k)] <- (1 / k[nrow(k), 2]) * e # all.equal(Wu[, nrow(k)], as.numeric(diag(nrow(y)) %*% Ri %*% e))
XWjxZWj <- t(W) %*% Wu * (1 / k[nrow(k), 2])
WiWj <- t(Wu) %*% Wu * (1 / k[nrow(k), 2])
buWu <- backsolve(Mchol_XZ, forwardsolve(t(Mchol_XZ), XWjxZWj))
AI <- (WiWj - (t(buWu) %*% XWjxZWj)) / 2 # all.equal(AI, (t(Wu) %*% P %*% Wu) / 2)
AIi <- solve(AI)
Ci <- chol2inv(Mchol_XZ) # all.equal(Ci, solve(C), check.attributes = FALSE)
traces <- matrix(0, nrow(sigmai), ncol(sigmai))
for (i in 1:nrow(sigmai)) {
  for (j in 1:ncol(sigmai)) {
    if (j >= i) { # Only do the upper triangle
      rstart <- (nX + 1 + (i - 1) * nrow(nvarGi))
      rend <- rstart - 1 + nrow(nvarGi)
      cstart <- (nX + 1 + (j - 1) * ncol(nvarGi))
      cend <- cstart - 1 + ncol(nvarGi)
      traces[i, j] <- sum(nvarGi * Ci[rstart:rend, cstart:cend])
    }
  }
}
traces <- traces + t(traces) - diag(diag(traces))
basisScores <- -0.5 * (nrow(nvarGi) * sigmai - sigmai %*% traces %*% sigmai - t(umat_sigmai) %*% nvarGi %*% umat_sigmai)
scores <- matrix(0, nrow(k), 1)
scores[1, 1] <- sum(derivs[[1]] * basisScores) # all.equal(scores[1, 1], as.numeric(-0.5 * (sum(diag(kronecker(derivs[[1]], K) %*% Gi)) - sum(diag((Gi %*% kronecker(derivs[[1]], K) %*% Gi) %*% Ci[5:988, 5:988])) - (t(y) %*% P %*% Wu[, 1, drop = FALSE]))))
scores[2, 1] <- sum(derivs[[2]] * basisScores) # all.equal(scores[2, 1], as.numeric(-0.5 * (sum(diag(kronecker(derivs[[2]], K) %*% Gi)) - sum(diag((Gi %*% kronecker(derivs[[2]], K) %*% Gi) %*% Ci[5:988, 5:988])) - (t(y) %*% P %*% Wu[, 2, drop = FALSE]))))
WSW <- t(W) %*% diag(nrow(y)) %*% W
evari <- (1 / k[nrow(k), 2])
scores[nrow(scores), 1] <- as.numeric(-0.5 * (sum(diag(nrow(y)) * evari) - sum(evari^2 * Ci * WSW) - (evari^2 * t(e) %*% diag(nrow(y)) %*% e)))
# all.equal(scores[nrow(scores), 1], as.numeric(-0.5 * (sum(diag(diag(nrow(y)) %*% Ri)) - sum(diag((t(cbind(X, Z)) %*% Ri %*% diag(nrow(y)) %*% Ri %*% cbind(X, Z)) %*% Ci)) - (t(y) %*% P %*% Wu[, nrow(k), drop = FALSE]))))
delta <- AIi %*% scores
k[, 3] <- k[, 2] + delta

## iteration 3 ====



