d <- readRDS("~/repos/kmmr/dev/weights_data.rds")
EC <- d$EC
ED <- d$ED
d <- d$ydata
K <- readRDS("~/repos/kmmr/dev/weights_K.rds")
Ki <- solve(K)
attr(Ki, "INVERSE") <- TRUE

d <- droplevels(d[d$Man == "LN" & d$Env %in% c("HAN_2015", "HAN_2016", "HAN_2017"),])

library(asreml)

fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ corgh(Env):vm(Gen, Ki),
                  residual = ~ idv(units),
                  data = d, maxit = 1)

fit.asr.wt <- asreml(fixed = GY ~ -1 + Env,
                     random = ~ corgh(Env):vm(Gen, Ki),
                     residual = ~ units,
                     data = d, maxit = 1,
                     weights = weight,
                     family = asr_gaussian(dispersion = 1.0))

fit.kmmr <- kmmr(fixed = GY ~ -1 + Env,
                 random = ~ re(~ corg(Env):lk(Gen, Ki), vars = ~ Env),
                 residual = ~ units,
                 data = d,
                 tolConvLL = 1e-8,
                 init = c(8, 8, 8, 0.5, 0.5, 0.5, 8))

fit.kmmr.wt <- kmmr(fixed = GY ~ -1 + Env,
                    random = ~ re(~ corg(Env):lk(Gen, Ki), vars = ~ Env),
                    residual = ~ units,
                    data = d,
                    weights = "weight",
                    init =    c(8, 8, 8, 0.5, 0.5, 0.5, 1),
                    isFixed = c(F, F, F, F,   F,   F,   TRUE),
                    tolConvLL = 1e-8)

summary(fit.asr)$varcomp
summary(fit.kmmr)$varcomp
summary(fit.asr.wt)$varcomp
summary(fit.kmmr.wt)$varcomp



