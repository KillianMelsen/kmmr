Ki <- solve(Matrix(readRDS("~/repos/kmmr/dev/exampledata_K.rds")))
Ki.asr <- as.matrix(Ki)
attr(Ki.asr, "INVERSE") <- TRUE
d <- readRDS("~/repos/kmmr/dev/exampledata.rds")
# d <- droplevels(d[d$Man == "R" & d$Env %in% c("Cam12", "Kar12"),])

fit.CPU <- kmmr(fixed = GY ~ -1 + ManEnv,
                random = ~ re(~ fa(ManEnv, 3):lk(Gen, Ki), vars = ~ ManEnv),
                residual = ~ units,
                data = d,
                sparse = FALSE,
                useGPU = FALSE)

fit.GPU <- kmmr(fixed = GY ~ -1 + ManEnv,
                random = ~ re(~ fa(ManEnv, 3):lk(Gen, Ki), vars = ~ ManEnv),
                residual = ~ units,
                data = d,
                sparse = FALSE,
                useGPU = TRUE)

fit.sparse <- kmmr(fixed = GY ~ -1 + ManEnv,
                   random = ~ re(~ exfa(ManEnv, 3):lk(Gen, Ki), vars = ~ ManEnv),
                   residual = ~ units,
                   data = d,
                   sparse = TRUE,
                   useGPU = FALSE)

library(asreml)
asreml.options(workspace = "5000mb")
d <- droplevels(d[d$Man == "R",])
fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ fa(Env, 3):vm(Gen, Ki.asr),
                  residual = ~ units,
                  data = d,
                  maxit = 20)

# This is really interesting, kmmr becomes much slower when I include an overall intercept.
# Seems like the Cholesky() reordering is screwed up because of the single dense column in X.
# When including an overall intercept, the CF has only 16 supernodes: 14 for the fixed effects, 1 smaller 246x246 one,
# and one absolutely massive supernode. There is a ton of fill-in. Not sure how to solve this, maybe by dividing C into
# a dense part and a sparse part?
fit.kmmr <- kmmr(fixed = GY ~ -1 + Env,
                 random = ~ re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
                 residual = ~ units,
                 data = d,
                 maxit = 20,
                 sparse = TRUE,
                 useGPU = FALSE)

# Trying to fix the fill-in issue ====
Ki <- solve(Matrix(readRDS("~/repos/kmmr/dev/exampledata_K.rds")))
d <- readRDS("~/repos/kmmr/dev/exampledata.rds")
set.seed(1997)
d$covariate <- rnorm(nrow(d), 5, 0.5)
fit.asr <- asreml(fixed = GY ~ 1 + ManEnv + covariate,
                  random = ~ fa(ManEnv, 3):vm(Gen, Ki.asr),
                  residual = ~ units,
                  data = d,
                  maxit = 20)

fit.sparse <- kmmr(fixed = GY ~ 1 + ManEnv + covariate,
                   random = ~ re(~ exfa(ManEnv, 3):lk(Gen, Ki), vars = ~ ManEnv),
                   residual = ~ units,
                   data = d,
                   sparse = TRUE,
                   useGPU = FALSE)

# Arguments:
fixed = GY ~ 1 + ManEnv + covariate
random = ~ re(~ exfa(ManEnv, 3):lk(Gen, Ki), vars = ~ ManEnv)
residual = ~ units
data = d
sparse = TRUE
useGPU = FALSE
maxit = 20
stepSize = NULL
avInfReg = NULL
weights = NULL
init = NULL
tolConvLL = 0.001
isFixed = NULL
verbose = TRUE
timesToBoundary = 5
getStartingValues = FALSE
super = TRUE




