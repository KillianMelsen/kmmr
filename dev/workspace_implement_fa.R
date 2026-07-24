d <- readRDS("dev/exampledata.rds")
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY")])
# d <- droplevels(d[, c("Man", "Env", "Gen", "GY")])
K <- readRDS("dev/exampledata_K.rds")
Ki <- Matrix(solve(K)) # Make sure any base matrices are in Matrix format!
attr(Ki, "logDet") <- sum(log(eigen(K)$values)) # We need the log determinant of K, NOT the inverse of K!!!
ED <- Matrix(readRDS("dev/exampledata_ED.rds"))
# EC <- readRDS("dev/exampledata_EC.rds")
# ECi <- solve(EC)

# Initial values for the fa(Env, ext = TRUE) parameters
vars <- diag(cov(unstack(d, GY ~ Env), use = "complete"))
fit.factanal <- factanal(covmat = cov(unstack(d, GY ~ Env), use = "complete"), factors = 3)
L <- fit.factanal$loadings[,]
L <- kmmr:::rotateInitL(L)
resvar <- mean(vars) / 20
genvar <- vars / 20
init <- c(genvar, matrix(L), resvar)
# init <- c(rep(genvar / 2, 2), 0.15, matrix(L), resvar)
isFixed <- rep(FALSE, length(init))
isFixed[init == 0] <- TRUE

system.time(
fit <- kmmr(fixed = GY ~ -1 + Env,
            random = ~ re(~ fa(Env, 3):lk(Gen, Ki), vars = ~ Env),
            residual = ~ units,
            data = d,
            init = init,
            isFixed = isFixed,
            sparse = FALSE)
)

Lhat.kmmr <- matrix(summary(fit)$varcomp$Estimate[15:56], 14, 3)
corhat.kmmr <- Lhat.kmmr %*% t(Lhat.kmmr); diag(corhat.kmmr) <- 1

library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = F, score = F, threads = -1)
Ki2 <- solve(K)
attr(Ki2, "INVERSE") <- TRUE
system.time(
fit2 <- asreml(fixed = GY ~ -1 + Env,
               random = ~ facv(Env, 3):vm(Gen, Ki2),
               residual = ~ units,
               data = d,
               maxit = 20)
)

Lhat.asr <- matrix(summary(fit2)$varcomp[1:42, "component"], 14, 3)
PSIhat.asr <- diag(summary(fit2)$varcomp[43:56, "component"])
corhat.asr <- cov2cor(Lhat.asr %*% t(Lhat.asr) + PSIhat.asr)

range(abs(corhat.asr - corhat.kmmr))

plot(fit2$coefficients$random, fit$u)
cor(fit2$coefficients$random, fit$u@x)
