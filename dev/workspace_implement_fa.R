library(devtools)
load_all()
d <- readRDS("dev/exampledata.rds")
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY")])
# d <- droplevels(d[, c("Man", "Env", "Gen", "GY")])
K <- readRDS("dev/exampledata_K.rds")
Ki <- Matrix(solve(K)) # Make sure any base matrices are in Matrix format!
attr(Ki, "logDet") <- sum(log(eigen(K)$values)) # We need the log determinant of K, NOT the inverse of K!!!
ED <- Matrix(readRDS("dev/exampledata_ED.rds"))
# EC <- Matrix(readRDS("dev/exampledata_EC.rds"))
# ECi <- solve(EC)

# Initial values for the fa(Env, ext = TRUE) parameters
# vars <- diag(cov(unstack(d, GY ~ Env), use = "complete"))
# fit.factanal <- factanal(covmat = cov(unstack(d, GY ~ Env), use = "complete"), factors = 3)
# L <- fit.factanal$loadings[,]
# L <- kmmr:::rotateInitL(L)
# resvar <- mean(vars) / 20
# genvar <- vars / 20
# init <- c(genvar, matrix(L), resvar)
# # init <- c(rep(genvar / 2, 2), matrix(L), resvar)
# # init <- c(rep(genvar / 2, 2), 0.15, matrix(L), resvar)
# isFixed <- rep(FALSE, length(init))
# isFixed[init == 0] <- TRUE

# Normal FA models:
system.time(
  fit <- kmmr(fixed = GY ~ -1 + Env,
              random = ~ re(~ fa(Env, 3):lk(Gen, Ki), vars = ~ Env),
              residual = ~ units,
              data = d,
              timesToBoundary = 7,
              getStartingValues = FALSE)
)

# system.time(
#   fit <- kmmr(fixed = GY ~ -1 + Env,
#               random = ~ re(~ stack(~ gk(Env, ED) + fa(Env, 2)):lk(Gen, Ki), vars = ~ Env),
#               residual = ~ units,
#               data = d,
#               timesToBoundary = 7,
#               stepSize = c(0.25, 0.35, 0.45, 0.5, 0.75, rep(1, 20)),
#               avInfReg = c(rep(1.05, 25)),
#               maxit = 25, getStartingValues = T)
# )

library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = F, score = F, threads = -1)
Ki2 <- solve(K)
attr(Ki2, "INVERSE") <- TRUE
system.time(
fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ facv(Env, 3):vm(Gen, Ki2),
                  residual = ~ idv(units),
                  data = d,
                  maxit = 25)
)

# system.time(
# fit <- kmmr(fixed = GY ~ -1 + Env,
#             random = ~ re(~ stack(~ gk(Env, ED) + fa(Env, 2)):lk(Gen, Ki), vars = ~ Env),
#             residual = ~ units,
#             data = d,
#             init = init,
#             isFixed = isFixed,
#             sparse = FALSE,
#             stepSize = c(0.25, 0.35, 0.45, 0.5, 0.75, rep(1, 20)),
#             avInfReg = c(rep(1.05, 25)),
#             maxit = 25)
# 
# fit2 <- kmmr(fixed = GY ~ -1 + Env,
#             random = ~ re(~ gk(Env, ED):lk(Gen, Ki), vars = ~ Env),
#             residual = ~ units,
#             data = d,
#             init = c(genvar, 0.15, resvar),
#             sparse = FALSE,
#             stepSize = c(0.25, 0.35, 0.45, 0.5, 0.75, rep(1, 20)),
#             avInfReg = c(rep(1.05, 25)),
#             maxit = 25)
# 
# fit3 <- kmmr(fixed = GY ~ -1 + Env,
#              random = ~ re(~ fa(Env, 2):lk(Gen, Ki), vars = ~ Env),
#              residual = ~ units,
#              data = d,
#              init = c(genvar, matrix(L), resvar),
#              isFixed = c(rep(F, 14), rep(F, 14), T, rep(F, 13), F),
#              sparse = FALSE,
#              stepSize = c(0.25, 0.35, 0.45, 0.5, 0.75, rep(1, 20)),
#              avInfReg = c(rep(1.05, 25)),
#              maxit = 25)
# )

# summary(fit)$logLik # stack of Gaussian kernel and FA
# summary(fit2)$logLik # Gaussian kernel only
# summary(fit3)$logLik # FA only
# 
# u.stack <- fit$u
# u.gk <- fit2$u
# u.fa <- fit3$u
# pairs(data.frame(u.stack@x, u.gk@x, u.fa@x))

fit <- kmmr(fixed = GY ~ -1 + Env,
            random = ~ re(~ stack(~ lk(Env, EC) + fa(Env, 2)):lk(Gen, Ki), vars = ~ Env),
            residual = ~ units,
            data = d,
            init = init,
            isFixed = isFixed,
            sparse = FALSE,
            stepSize = c(0.25, 0.35, 0.45, 0.5, 0.75, rep(1, 20)),
            avInfReg = c(rep(1.05, 25)),
            maxit = 25)

Lhat.kmmr <- matrix(summary(fit)$varcomp$Estimate[29:56], 14, 2)
corhat.kmmr <- Lhat.kmmr %*% t(Lhat.kmmr); diag(corhat.kmmr) <- 1
covhat.kmmr <- outer(sqrt(summary(fit)$varcomp$Estimate[15:28]), sqrt(summary(fit)$varcomp$Estimate[15:28])) * corhat.kmmr

library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = F, score = F, threads = -1)
Ki2 <- solve(K)
attr(Ki2, "INVERSE") <- TRUE
# system.time(
fit2 <- asreml(fixed = GY ~ -1 + Env,
               random = ~ corgh(Env):vm(Gen, Ki2) + facv(Env, 2):vm(Gen, Ki2),
               residual = ~ idv(units),
               data = d,
               maxit = 25,
               start.values = TRUE)

vc <- fit2$vparameters.table

vc$Value[1:91] <- EC[upper.tri(EC)]
vc$Constraint[1:91] <- "F"
fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ corgh(Env):vm(Gen, Ki2) + facv(Env, 2):vm(Gen, Ki2),
                  residual = ~ idv(units),
                  data = d,
                  maxit = 25,
                  G.param = vc)
# )

Lhat.asr <- matrix(summary(fit.asr)$varcomp[106:133, "component"], 14, 2)
PSIhat.asr <- diag(summary(fit.asr)$varcomp[134:147, "component"])
covhat.asr <- Lhat.asr %*% t(Lhat.asr) + PSIhat.asr
corhat.asr <- cov2cor(covhat.asr)

lkhat.asr <- outer(sqrt(summary(fit.asr)$varcomp[92:105, "component"]), sqrt(summary(fit.asr)$varcomp[92:105, "component"])) * EC
lkhat.kmmr <- outer(sqrt(summary(fit)$varcomp$Estimate[1:14]), sqrt(summary(fit)$varcomp$Estimate[1:14])) * EC

range(abs(corhat.asr - corhat.kmmr))
range(abs(lkhat.asr - lkhat.kmmr))
total.asr <- lkhat.asr + covhat.asr
total.kmmr <- lkhat.kmmr + covhat.kmmr
range(abs(cov2cor(total.asr) - cov2cor(total.kmmr)))
mean(abs(cov2cor(total.asr) - cov2cor(total.kmmr)))

u.asr <- fit.asr$coefficients$random[1:3444] + fit.asr$coefficients$random[3445:6888]
u.kmmr <- fit$u@x
plot(u.asr, u.kmmr)
cor(u.asr, u.kmmr)
