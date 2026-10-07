d <- readRDS("~/repos/kmmr/dev/exampledata.weights.rds")
K <- readRDS("~/repos/kmmr/dev/exampledata_K.rds")[levels(d$Gen), levels(d$Gen)]
Ki.asr <- solve(K)
Ki.kmr <- Matrix(Ki.asr)
attr(Ki.asr, "INVERSE") <- TRUE
d$weight.kmr <- 1 / d$weight

fit.asr <- asreml(fixed = GY ~ -1 + ManEnv,
                  random = ~ fa(ManEnv, 1):vm(Gen, Ki.asr),
                  residual = ~ units,
                  data = d, maxit = 50,
                  weights = weight, family = asr_gaussian(dispersion = 1.0))

sv <- kmmr(fixed = GY ~ -1 + ManEnv,
           random = ~ re(~ exfa(ManEnv, 1):lk(Gen, Ki.kmr), vars = ~ ManEnv),
           residual = ~ units,
           data = d, maxit = 50, getStartingValues = TRUE)
sv$Init[nrow(sv)] <- 1
sv$Fixed[nrow(sv)] <- T

fit.kmr <- kmmr(fixed = GY ~ -1 + ManEnv,
                random = ~ re(~ exfa(ManEnv, 1):lk(Gen, Ki.kmr), vars = ~ ManEnv),
                residual = ~ units,
                data = d, maxit = 50, weights = "weight.kmr", init = sv$Init, isFixed = sv$Fixed)
