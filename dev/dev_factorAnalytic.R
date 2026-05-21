library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = TRUE, score = TRUE)

d <- readRDS("~/repos/kmmr/dev/weights_data.rds")
EC <- d$EC
ED <- d$ED
d <- d$ydata
K <- readRDS("~/repos/kmmr/dev/weights_K.rds")
Ki <- solve(K)
attr(Ki, "INVERSE") <- TRUE

# excl <- c("KIE_2017", "QLB_2016")
# d <- droplevels(d[d$Man == "LN" & !(d$Env %in% excl),])
d <- droplevels(d[d$Man == "LN" & d$Env %in% c("HAN_2015", "KAL_2015", "KIE_2015", "QLB_2015", "RHH_2015"),])

fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ facv(Env, 2):vm(Gen, Ki),
                  residual = ~ idv(units),
                  data = d,
                  maxit = 1)

# Manual:
p <- length(levels(d$Env))
init <- as.numeric(fit.asr$trace[5:20,1,drop=F])
isFixed <- init == 0
k <- matrix(NA, 16, 20)
LL <- numeric(20)
k[, 1] <- init
X <- model.matrix(GY ~ -1 + Env, d)
Z <- diag(nrow(d))
W <- cbind(X, Z)
y <- matrix(d$GY)

getLoadingDerivs <- function(rowL, colL, L) {
  deriv <- matrix(0, nrow(L), nrow(L))
  deriv[rowL,] <- deriv[,rowL] <- L[, colL]
  deriv[rowL, rowL] <- 2 * L[rowL, colL]
  return(deriv)
}

iter <- 1
R <- k[nrow(k), iter] * diag(nrow(d))
logDetR <- sum(log(eigen(R)$values))
Ri <- solve(R)
L <- matrix(k[1:10, iter], 5, 2)
PSI <- diag(k[11:15, iter])
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
LL[iter] <- -0.5 * (logDetC + logDetR + logDetG + t(y) %*% P %*% y)
RHS <- rbind(t(X) %*% Ri %*% y,
             t(Z) %*% Ri %*% y)
bu <- Ci %*% RHS
b <- bu[1:5, 1, drop = F]
u <- bu[6:1015, 1, drop = F]
e <- y - W %*% bu

dGdL11 <- kronecker(getLoadingDerivs(1, 1, L), K)
dGdL21 <- kronecker(getLoadingDerivs(2, 1, L), K)
dGdL31 <- kronecker(getLoadingDerivs(3, 1, L), K)
dGdL41 <- kronecker(getLoadingDerivs(4, 1, L), K)
dGdL51 <- kronecker(getLoadingDerivs(5, 1, L), K)

dGdL12 <- kronecker(getLoadingDerivs(1, 2, L), K)
dGdL22 <- kronecker(getLoadingDerivs(2, 2, L), K)
dGdL32 <- kronecker(getLoadingDerivs(3, 2, L), K)
dGdL42 <- kronecker(getLoadingDerivs(4, 2, L), K)
dGdL52 <- kronecker(getLoadingDerivs(5, 2, L), K)

dGdPSI1 <- matrix(0, nrow(L), nrow(L)); dGdPSI1[1, 1] <- 1; dGdPSI1 <- kronecker(dGdPSI1, K)
dGdPSI2 <- matrix(0, nrow(L), nrow(L)); dGdPSI2[2, 2] <- 1; dGdPSI2 <- kronecker(dGdPSI2, K)
dGdPSI3 <- matrix(0, nrow(L), nrow(L)); dGdPSI3[3, 3] <- 1; dGdPSI3 <- kronecker(dGdPSI3, K)
dGdPSI4 <- matrix(0, nrow(L), nrow(L)); dGdPSI4[4, 4] <- 1; dGdPSI4 <- kronecker(dGdPSI4, K)
dGdPSI5 <- matrix(0, nrow(L), nrow(L)); dGdPSI5[5, 5] <- 1; dGdPSI5 <- kronecker(dGdPSI5, K)

Wu <- matrix(0, nrow(d), nrow(k))
Wu[, 1] <- dGdL11 %*% Gi %*% u
Wu[, 2] <- dGdL21 %*% Gi %*% u
Wu[, 3] <- dGdL31 %*% Gi %*% u
Wu[, 4] <- dGdL41 %*% Gi %*% u
Wu[, 5] <- dGdL51 %*% Gi %*% u

Wu[, 6] <-  0
Wu[, 7] <-  dGdL22 %*% Gi %*% u
Wu[, 8] <-  dGdL32 %*% Gi %*% u
Wu[, 9] <-  dGdL42 %*% Gi %*% u
Wu[, 10] <- dGdL52 %*% Gi %*% u

Wu[, 11] <- dGdPSI1 %*% Gi %*% u
Wu[, 12] <- dGdPSI2 %*% Gi %*% u
Wu[, 13] <- dGdPSI3 %*% Gi %*% u
Wu[, 14] <- dGdPSI4 %*% Gi %*% u
Wu[, 15] <- dGdPSI5 %*% Gi %*% u

Wu[, 16] <- Ri %*% e

AI <- t(Wu) %*% P %*% Wu / 2
AIi <- matrix(0, nrow(AI), ncol(AI))
AIi[!isFixed, !isFixed] <- solve(AI[!isFixed, !isFixed])

scores <- matrix(0, nrow(k), 1)
scores[1, 1] <-  -0.5 * (sum(diag(dGdL11 %*% Gi)) - sum(diag((Gi %*% dGdL11 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 1])
scores[2, 1] <-  -0.5 * (sum(diag(dGdL21 %*% Gi)) - sum(diag((Gi %*% dGdL21 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 2])
scores[3, 1] <-  -0.5 * (sum(diag(dGdL31 %*% Gi)) - sum(diag((Gi %*% dGdL31 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 3])
scores[4, 1] <-  -0.5 * (sum(diag(dGdL41 %*% Gi)) - sum(diag((Gi %*% dGdL41 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 4])
scores[5, 1] <-  -0.5 * (sum(diag(dGdL51 %*% Gi)) - sum(diag((Gi %*% dGdL51 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 5])

scores[6, 1] <-  0
scores[7, 1] <-  -0.5 * (sum(diag(dGdL22 %*% Gi)) - sum(diag((Gi %*% dGdL22 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 7])
scores[8, 1] <-  -0.5 * (sum(diag(dGdL32 %*% Gi)) - sum(diag((Gi %*% dGdL32 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 8])
scores[9, 1] <-  -0.5 * (sum(diag(dGdL42 %*% Gi)) - sum(diag((Gi %*% dGdL42 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 9])
scores[10, 1] <- -0.5 * (sum(diag(dGdL52 %*% Gi)) - sum(diag((Gi %*% dGdL52 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 10])

scores[11, 1] <- -0.5 * (sum(diag(dGdPSI1 %*% Gi)) - sum(diag((Gi %*% dGdPSI1 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 11])
scores[12, 1] <- -0.5 * (sum(diag(dGdPSI2 %*% Gi)) - sum(diag((Gi %*% dGdPSI2 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 12])
scores[13, 1] <- -0.5 * (sum(diag(dGdPSI3 %*% Gi)) - sum(diag((Gi %*% dGdPSI3 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 13])
scores[14, 1] <- -0.5 * (sum(diag(dGdPSI4 %*% Gi)) - sum(diag((Gi %*% dGdPSI4 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 14])
scores[15, 1] <- -0.5 * (sum(diag(dGdPSI5 %*% Gi)) - sum(diag((Gi %*% dGdPSI5 %*% Gi) %*% Ci[6:1015, 6:1015])) - t(y) %*% P %*% Wu[, 15])

scores[16, 1] <- -0.5 * (sum(diag(Ri)) - sum(diag((t(W) %*% Ri %*% Ri %*% W) %*% Ci)) - t(y) %*% P %*% Wu[, 16])

delta <- fit.asr$trace[4,1] * AIi %*% scores

# Our delta does not match asreml's delta.
# The scores are identical, the AIi matrix is different, and the working variates also seem to be different. LL is also identical.
# Also, the stepsize seems to be squared.
AIi.asr <- as.matrix(fit.asr$ai[1:16, 1:16]); rownames(AIi.asr) <- colnames(AIi.asr) <- NULL
Wu.asr <- as.matrix(fit.asr$wrkv)

delta <- fit.asr$trace[4,1]^2 * AIi %*% scores
delta.asr <- fit.asr$trace[4,1]^2 * AIi.asr %*% scores
delta.asr.actual <- matrix(fit.asr$trace[5:20, 2] - fit.asr$trace[5:20, 1])




