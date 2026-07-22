d <- readRDS("dev/exampledata.rds")
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY")])
# d <- droplevels(d[, c("Man", "Env", "Gen", "GY")])
K <- readRDS("dev/exampledata_K.rds")
Ki <- solve(Matrix(K))
# ED <- readRDS("dev/exampledata_ED.rds")
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
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = TRUE, score = TRUE)
system.time(
fit2 <- asreml(fixed = GY ~ -1 + Env,
               random = ~ facv(Env, 3):vm(Gen, K),
               residual = ~ units,
               data = d)
)

Lhat.asr <- matrix(summary(fit2)$varcomp[1:42, "component"], 14, 3)
PSIhat.asr <- diag(summary(fit2)$varcomp[43:56, "component"])
corhat.asr <- cov2cor(Lhat.asr %*% t(Lhat.asr) + PSIhat.asr)

range(abs(corhat.asr - corhat.kmmr))

plot(fit2$coefficients$random, fit$u)

