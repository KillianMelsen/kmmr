# library(asreml)
# asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = TRUE, score = TRUE, step.size = 0.1)

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

ds <- droplevels(d[d$Env %in% c("Bol12", "Kar12", "Cam12", "Deb12"),])
ds$ManGen <- as.factor(paste(as.character(ds$Man), as.character(ds$Gen), sep = ":"))
combis <- c(paste(levels(ds$Env)[1], levels(ds$Gen), sep = ":"),
            paste(levels(ds$Env)[2], levels(ds$Gen), sep = ":"),
            paste(levels(ds$Env)[3], levels(ds$Gen), sep = ":"),
            paste(levels(ds$Env)[4], levels(ds$Gen), sep = ":"))
ds$EnvGen <- factor(combis, combis, combis)

# fit <- asreml(fixed = GY ~ 1,
#               random = ~ idv(Man):vm(Gen, Ki) + idv(Env):vm(Gen, Ki),
#               residual = ~ idv(units),
#               data = ds, maxit = 25)
# 
# summary(fit)$varcomp
# fit$trace

fit.kmmr <- kmmr(fixed = GY ~ 1,
                 random = ~ re(~ id(Man):lk(Gen, Ki), vars = ~ 1) + re(~ id(Env):lk(Gen, Ki), vars = ~ 1),
                 residual = ~ units,
                 data = ds,
                 maxit = 25)




# # Running manually through the kmmr solver ====
# 
# # Data preparation:
# y <- matrix(ds$GY)
# X <- model.matrix(~ 1, ds)
# Z1 <- model.matrix(~ -1 + ManGen, ds)
# Z2 <- model.matrix(~ -1 + EnvGen, ds)
# Z <- cbind(Z1, Z2)
# nX <- ncol(X)
# 
# # Storage for LL and parameters:
# LL <- numeric(10)
# k <- matrix(NA, 3, 10)
# k[, 1] <- as.numeric(fit$trace[5:7, 1])
# # k[, 1] <- summary(fit)$varcomp[1:3, 1]
# 
# # Log|K|:
# logDetnvarGtotal <- vector("list", 2)
# logDetnvarGtotal[[1]] <- logDetnvarGtotal[[2]] <- -1 * sum(log(eigen(Ki)$values))
# 
# # Non-variable part of Gi:
# nvarGi <- vector("list", 2)
# nvarGi[[1]] <- nvarGi[[2]] <- Ki
# 
# ## Iteration 1 ====
# iter <- 1
# Ri <- 1 / k[nrow(k), iter] * diag(nrow(y))
# W <- cbind(X, Z1, Z2)
# Wy <- cbind(W, y)
# M0 <- t(Wy) %*% Wy
# M <- M0 * (1 / k[nrow(k), iter])
# sigmai <- vector("list", 2)
# sigmai[[1]] <- 1 / k[1, iter] * diag(2)
# sigmai[[2]] <- 1 / k[2, iter] * diag(4)
# logDetvarGtotal <- vector("list", 2)
# logDetvarGtotal[[1]] <- -1 * sum(log(eigen(sigmai[[1]])$values))
# logDetvarGtotal[[2]] <- -1 * sum(log(eigen(sigmai[[2]])$values))
# Gi <- vector("list", 2)
# Gi[[1]] <- kronecker(sigmai[[1]], nvarGi[[1]])
# Gi[[2]] <- kronecker(sigmai[[2]], nvarGi[[2]])
# M[2:493, 2:493] <- M[2:493, 2:493] + Gi[[1]]
# M[494:1477, 494:1477] <- M[494:1477, 494:1477] + Gi[[2]]
# Mchol <- chol(M)
# yPy <- Mchol[nrow(Mchol), ncol(Mchol)]^2 # R <- solve(Ri); G1 <- solve(Gi[[1]]); G2 <- solve(Gi[[2]]); Vi <- solve(Z1 %*% G1 %*% t(Z1) + Z2 %*% G2 %*% t(Z2) + R); P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi; all.equal(yPy, as.numeric(t(y) %*% P %*% y))
# Mchol_XZ <- Mchol[1:(nrow(Mchol) - 1), 1:(ncol(Mchol) - 1)]
# My <- Mchol[1:(nrow(Mchol) - 1), ncol(Mchol):ncol(Mchol), drop = FALSE]
# logDetC <- 2 * sum(log(diag(Mchol_XZ)))
# # C <- rbind(cbind(t(X) %*% Ri %*% X,    t(X) %*% Ri %*% Z1,               t(X) %*% Ri %*% Z2),
# #            cbind(t(Z1) %*% Ri %*% X,   t(Z1) %*% Ri %*% Z1 + Gi[[1]],    t(Z1) %*% Ri %*% Z2),
# #            cbind(t(Z2) %*% Ri %*% X,   t(Z2) %*% Ri %*% Z1,              t(Z2) %*% Ri %*% Z2 + Gi[[2]]))
# # all.equal(C, M[1:nrow(C), 1:ncol(C)])
# # all.equal(logDetC, sum(log(eigen(C)$values)))
# 
# logDetG <- (nrow(nvarGi[[1]]) * logDetvarGtotal[[1]] + nrow(sigmai[[1]]) * logDetnvarGtotal[[1]]) + (nrow(nvarGi[[2]]) * logDetvarGtotal[[2]] + nrow(sigmai[[2]]) * logDetnvarGtotal[[2]])
# # all.equal(logDetG, sum(log(eigen(kronecker(solve(sigmai[[1]]), K))$values)) + sum(log(eigen(kronecker(solve(sigmai[[2]]), K))$values)))
# 
# logDetR <- nrow(y) * log(k[nrow(k), iter]) # all.equal(logDetR, sum(log(eigen(solve(Ri))$values)))
# LL[iter] <- -0.5 * (logDetC + logDetR + logDetG + yPy)
# bu <- backsolve(Mchol_XZ, My)
# b <- bu[1, 1, drop = FALSE]
# u <- vector("list", 2)
# u[[1]] <- bu[2:493, 1, drop = FALSE]
# u[[2]] <- bu[494:1477, 1, drop = FALSE]
# Wu <- matrix(NA, nrow(y), nrow(k))
# 
# umat_sigmai <- vector("list", 2)
# derivs <- vector("list", 2)
# # First random effect:
# umat <- matrix(u[[1]], nrow = nrow(nvarGi[[1]]), ncol = nrow(sigmai[[1]]))
# umat_sigmai[[1]] <- umat %*% sigmai[[1]]
# derivs[[1]] <- vector("list", 1)
# derivs[[1]][[1]] <- diag(2)
# Wu[, 1] <- Z1 %*% matrix(umat_sigmai[[1]] %*% derivs[[1]][[1]]) # all.equal(Wu[, 1], as.numeric(Z1 %*% kronecker(derivs[[1]][[1]], K) %*% Gi[[1]] %*% u[[1]]))
# 
# # Second random effect:
# umat <- matrix(u[[2]], nrow = nrow(nvarGi[[2]]), ncol = nrow(sigmai[[2]]))
# umat_sigmai[[2]] <- umat %*% sigmai[[2]]
# derivs[[2]] <- vector("list", 1)
# derivs[[2]][[1]] <- diag(4)
# Wu[, 2] <- Z2 %*% matrix(umat_sigmai[[2]] %*% derivs[[2]][[1]]) # all.equal(Wu[, 2], as.numeric(Z2 %*% kronecker(derivs[[2]][[1]], K) %*% Gi[[2]] %*% u[[2]]))
# 
# # Residuals:
# e <- y - W %*% bu
# Wu[, nrow(k)] <- (1 / k[nrow(k), iter]) * e # all.equal(Wu[, nrow(k)], as.numeric(diag(nrow(y)) %*% Ri %*% e))
# XWjxZWj <- t(W) %*% Wu * (1 / k[nrow(k), iter])
# WiWj <- t(Wu) %*% Wu * (1 / k[nrow(k), iter])
# buWu <- backsolve(Mchol_XZ, forwardsolve(t(Mchol_XZ), XWjxZWj))
# AI <- (WiWj - (t(buWu) %*% XWjxZWj)) / 2 # all.equal(AI, (t(Wu) %*% P %*% Wu) / 2)
# AIi <- solve(AI)
# Ci <- chol2inv(Mchol_XZ) # all.equal(Ci, solve(C), check.attributes = FALSE)
# 
# scores <- matrix(0, nrow(k), 1)
# # First random effect:
# traces <- matrix(0, nrow(sigmai[[1]]), ncol(sigmai[[1]]))
# for (i in 1:nrow(sigmai[[1]])) {
#   for (j in 1:ncol(sigmai[[1]])) {
#     if (j >= i) { # Only do the upper triangle
#       rstart <- (nX + 1 + (i - 1) * nrow(nvarGi[[1]]))
#       rend <- rstart - 1 + nrow(nvarGi[[1]])
#       cstart <- (nX + 1 + (j - 1) * ncol(nvarGi[[1]]))
#       cend <- cstart - 1 + ncol(nvarGi[[1]])
#       traces[i, j] <- sum(nvarGi[[1]] * Ci[rstart:rend, cstart:cend])
#     }
#   }
# }
# traces <- traces + t(traces) - diag(diag(traces))
# basisScores <- -0.5 * (nrow(nvarGi[[1]]) * sigmai[[1]] - sigmai[[1]] %*% traces %*% sigmai[[1]] - t(umat_sigmai[[1]]) %*% nvarGi[[1]] %*% umat_sigmai[[1]])
# scores[1, 1] <- sum(derivs[[1]][[1]] * basisScores) # all.equal(scores[1, 1], as.numeric(-0.5 * (sum(diag(kronecker(derivs[[1]][[1]], K) %*% Gi[[1]])) - sum(diag((Gi[[1]] %*% kronecker(derivs[[1]][[1]], K) %*% Gi[[1]]) %*% Ci[2:493, 2:493])) - (t(y) %*% P %*% Wu[, 1, drop = FALSE]))))
# 
# # Second random effect:
# traces <- matrix(0, nrow(sigmai[[2]]), ncol(sigmai[[2]]))
# for (i in 1:nrow(sigmai[[2]])) {
#   for (j in 1:ncol(sigmai[[2]])) {
#     if (j >= i) { # Only do the upper triangle
#       rstart <- (nX + 492 + 1 + (i - 1) * nrow(nvarGi[[2]]))
#       rend <- rstart - 1 + nrow(nvarGi[[2]])
#       cstart <- (nX + 492 + 1 + (j - 1) * ncol(nvarGi[[2]]))
#       cend <- cstart - 1 + ncol(nvarGi[[2]])
#       traces[i, j] <- sum(nvarGi[[2]] * Ci[rstart:rend, cstart:cend])
#     }
#   }
# }
# traces <- traces + t(traces) - diag(diag(traces))
# basisScores <- -0.5 * (nrow(nvarGi[[2]]) * sigmai[[2]] - sigmai[[2]] %*% traces %*% sigmai[[2]] - t(umat_sigmai[[2]]) %*% nvarGi[[2]] %*% umat_sigmai[[2]])
# scores[2, 1] <- sum(derivs[[2]][[1]] * basisScores) # all.equal(scores[2, 1], as.numeric(-0.5 * (sum(diag(kronecker(derivs[[2]][[1]], K) %*% Gi[[2]])) - sum(diag((Gi[[2]] %*% kronecker(derivs[[2]][[1]], K) %*% Gi[[2]]) %*% Ci[494:1477, 494:1477])) - (t(y) %*% P %*% Wu[, 2, drop = FALSE]))))
# 
# # Residual:
# WSW <- t(W) %*% diag(nrow(y)) %*% W
# evari <- (1 / k[nrow(k), iter])
# scores[nrow(scores), 1] <- as.numeric(-0.5 * (sum(diag(nrow(y)) * evari) - sum(evari^2 * Ci * WSW) - (evari^2 * t(e) %*% diag(nrow(y)) %*% e)))
# # all.equal(scores[nrow(scores), 1], as.numeric(-0.5 * (sum(diag(diag(nrow(y)) %*% Ri)) - sum(diag((t(cbind(X, Z1, Z2)) %*% Ri %*% diag(nrow(y)) %*% Ri %*% cbind(X, Z1, Z2)) %*% Ci)) - (t(y) %*% P %*% Wu[, nrow(k), drop = FALSE]))))
# delta <- AIi %*% scores * 0.3162278^2
# # delta - matrix(fit$trace[5:7, 2] - fit$trace[5:7, 1])
# k[, iter + 1] <- k[, iter] + delta
# 
# ## Iteration 2 ====
# iter <- 2
# Ri <- 1 / k[nrow(k), iter] * diag(nrow(y))
# W <- cbind(X, Z1, Z2)
# Wy <- cbind(W, y)
# M0 <- t(Wy) %*% Wy
# M <- M0 * (1 / k[nrow(k), iter])
# sigmai <- vector("list", 2)
# sigmai[[1]] <- k[1, iter] * diag(2)
# sigmai[[2]] <- k[2, iter] * diag(4)
# logDetvarGtotal <- vector("list", 2)
# logDetvarGtotal[[1]] <- -1 * sum(log(eigen(sigmai[[1]])$values))
# logDetvarGtotal[[2]] <- -1 * sum(log(eigen(sigmai[[2]])$values))
# Gi <- vector("list", 2)
# Gi[[1]] <- kronecker(sigmai[[1]], nvarGi[[1]])
# Gi[[2]] <- kronecker(sigmai[[2]], nvarGi[[2]])
# M[2:493, 2:493] <- M[2:493, 2:493] + Gi[[1]]
# M[494:1477, 494:1477] <- M[494:1477, 494:1477] + Gi[[2]]
# Mchol <- chol(M)
# yPy <- Mchol[nrow(Mchol), ncol(Mchol)]^2 # R <- solve(Ri); G1 <- solve(Gi[[1]]); G2 <- solve(Gi[[2]]); Vi <- solve(Z1 %*% G1 %*% t(Z1) + Z2 %*% G2 %*% t(Z2) + R); P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X) %*% t(X) %*% Vi; all.equal(yPy, as.numeric(t(y) %*% P %*% y))
# Mchol_XZ <- Mchol[1:(nrow(Mchol) - 1), 1:(ncol(Mchol) - 1)]
# My <- Mchol[1:(nrow(Mchol) - 1), ncol(Mchol):ncol(Mchol), drop = FALSE]
# logDetC <- 2 * sum(log(diag(Mchol_XZ)))
# # C <- rbind(cbind(t(X) %*% Ri %*% X,    t(X) %*% Ri %*% Z1,               t(X) %*% Ri %*% Z2),
# #            cbind(t(Z1) %*% Ri %*% X,   t(Z1) %*% Ri %*% Z1 + Gi[[1]],    t(Z1) %*% Ri %*% Z2),
# #            cbind(t(Z2) %*% Ri %*% X,   t(Z2) %*% Ri %*% Z1,              t(Z2) %*% Ri %*% Z2 + Gi[[2]]))
# # all.equal(C, M[1:nrow(C), 1:ncol(C)])
# # all.equal(logDetC, sum(log(eigen(C)$values)))
# 
# logDetG <- (nrow(nvarGi[[1]]) * logDetvarGtotal[[1]] + nrow(sigmai[[1]]) * logDetnvarGtotal[[1]]) + (nrow(nvarGi[[2]]) * logDetvarGtotal[[2]] + nrow(sigmai[[2]]) * logDetnvarGtotal[[2]])
# # all.equal(logDetG, sum(log(eigen(kronecker(solve(sigmai[[1]]), K))$values)) + sum(log(eigen(kronecker(solve(sigmai[[2]]), K))$values)))
# 
# logDetR <- nrow(y) * log(k[nrow(k), iter]) # all.equal(logDetR, sum(log(eigen(solve(Ri))$values)))
# LL[iter] <- -0.5 * (logDetC + logDetR + logDetG + yPy)
# bu <- backsolve(Mchol_XZ, My)
# b <- bu[1, 1, drop = FALSE]
# u <- vector("list", 2)
# u[[1]] <- bu[2:493, 1, drop = FALSE]
# u[[2]] <- bu[494:1477, 1, drop = FALSE]
# Wu <- matrix(NA, nrow(y), nrow(k))
# 
# umat_sigmai <- vector("list", 2)
# derivs <- vector("list", 2)
# # First random effect:
# umat <- matrix(u[[1]], nrow = nrow(nvarGi[[1]]), ncol = nrow(sigmai[[1]]))
# umat_sigmai[[1]] <- umat %*% sigmai[[1]]
# derivs[[1]] <- vector("list", 1)
# derivs[[1]][[1]] <- diag(2)
# Wu[, 1] <- Z1 %*% matrix(umat_sigmai[[1]] %*% derivs[[1]][[1]]) # all.equal(Wu[, 1], as.numeric(Z1 %*% kronecker(derivs[[1]][[1]], K) %*% Gi[[1]] %*% u[[1]]))
# 
# # Second random effect:
# umat <- matrix(u[[2]], nrow = nrow(nvarGi[[2]]), ncol = nrow(sigmai[[2]]))
# umat_sigmai[[2]] <- umat %*% sigmai[[2]]
# derivs[[2]] <- vector("list", 1)
# derivs[[2]][[1]] <- diag(4)
# Wu[, 2] <- Z2 %*% matrix(umat_sigmai[[2]] %*% derivs[[2]][[1]]) # all.equal(Wu[, 2], as.numeric(Z2 %*% kronecker(derivs[[2]][[1]], K) %*% Gi[[2]] %*% u[[2]]))
# 
# # Residuals:
# e <- y - W %*% bu
# Wu[, nrow(k)] <- (1 / k[nrow(k), iter]) * e # all.equal(Wu[, nrow(k)], as.numeric(diag(nrow(y)) %*% Ri %*% e))
# XWjxZWj <- t(W) %*% Wu * (1 / k[nrow(k), iter])
# WiWj <- t(Wu) %*% Wu * (1 / k[nrow(k), iter])
# buWu <- backsolve(Mchol_XZ, forwardsolve(t(Mchol_XZ), XWjxZWj))
# AI <- (WiWj - (t(buWu) %*% XWjxZWj)) / 2 # all.equal(AI, (t(Wu) %*% P %*% Wu) / 2)
# AIi <- solve(AI)
# Ci <- chol2inv(Mchol_XZ) # all.equal(Ci, solve(C), check.attributes = FALSE)
# 
# scores <- matrix(0, nrow(k), 1)
# # First random effect:
# traces <- matrix(0, nrow(sigmai[[1]]), ncol(sigmai[[1]]))
# for (i in 1:nrow(sigmai[[1]])) {
#   for (j in 1:ncol(sigmai[[1]])) {
#     if (j >= i) { # Only do the upper triangle
#       rstart <- (nX + 1 + (i - 1) * nrow(nvarGi[[1]]))
#       rend <- rstart - 1 + nrow(nvarGi[[1]])
#       cstart <- (nX + 1 + (j - 1) * ncol(nvarGi[[1]]))
#       cend <- cstart - 1 + ncol(nvarGi[[1]])
#       traces[i, j] <- sum(nvarGi[[1]] * Ci[rstart:rend, cstart:cend])
#     }
#   }
# }
# traces <- traces + t(traces) - diag(diag(traces))
# basisScores <- -0.5 * (nrow(nvarGi[[1]]) * sigmai[[1]] - sigmai[[1]] %*% traces %*% sigmai[[1]] - t(umat_sigmai[[1]]) %*% nvarGi[[1]] %*% umat_sigmai[[1]])
# scores[1, 1] <- sum(derivs[[1]][[1]] * basisScores) # all.equal(scores[1, 1], as.numeric(-0.5 * (sum(diag(kronecker(derivs[[1]][[1]], K) %*% Gi[[1]])) - sum(diag((Gi[[1]] %*% kronecker(derivs[[1]][[1]], K) %*% Gi[[1]]) %*% Ci[2:493, 2:493])) - (t(y) %*% P %*% Wu[, 1, drop = FALSE]))))
# 
# # Second random effect:
# traces <- matrix(0, nrow(sigmai[[2]]), ncol(sigmai[[2]]))
# for (i in 1:nrow(sigmai[[2]])) {
#   for (j in 1:ncol(sigmai[[2]])) {
#     if (j >= i) { # Only do the upper triangle
#       rstart <- (nX + 492 + 1 + (i - 1) * nrow(nvarGi[[2]]))
#       rend <- rstart - 1 + nrow(nvarGi[[2]])
#       cstart <- (nX + 492 + 1 + (j - 1) * ncol(nvarGi[[2]]))
#       cend <- cstart - 1 + ncol(nvarGi[[2]])
#       traces[i, j] <- sum(nvarGi[[2]] * Ci[rstart:rend, cstart:cend])
#     }
#   }
# }
# traces <- traces + t(traces) - diag(diag(traces))
# basisScores <- -0.5 * (nrow(nvarGi[[2]]) * sigmai[[2]] - sigmai[[2]] %*% traces %*% sigmai[[2]] - t(umat_sigmai[[2]]) %*% nvarGi[[2]] %*% umat_sigmai[[2]])
# scores[2, 1] <- sum(derivs[[2]][[1]] * basisScores) # all.equal(scores[2, 1], as.numeric(-0.5 * (sum(diag(kronecker(derivs[[2]][[1]], K) %*% Gi[[2]])) - sum(diag((Gi[[2]] %*% kronecker(derivs[[2]][[1]], K) %*% Gi[[2]]) %*% Ci[494:1477, 494:1477])) - (t(y) %*% P %*% Wu[, 2, drop = FALSE]))))
# 
# # Residual:
# WSW <- t(W) %*% diag(nrow(y)) %*% W
# evari <- (1 / k[nrow(k), iter])
# scores[nrow(scores), 1] <- as.numeric(-0.5 * (sum(diag(nrow(y)) * evari) - sum(evari^2 * Ci * WSW) - (evari^2 * t(e) %*% diag(nrow(y)) %*% e)))
# # all.equal(scores[nrow(scores), 1], as.numeric(-0.5 * (sum(diag(diag(nrow(y)) %*% Ri)) - sum(diag((t(cbind(X, Z1, Z2)) %*% Ri %*% diag(nrow(y)) %*% Ri %*% cbind(X, Z1, Z2)) %*% Ci)) - (t(y) %*% P %*% Wu[, nrow(k), drop = FALSE]))))
# delta <- AIi %*% scores * 0.1
# k[, iter + 1] <- k[, iter] + delta
# 
# ## iteration 3 ====



