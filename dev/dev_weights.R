d <- readRDS("~/repos/kmmr/dev/weights_data.rds")
EC <- d$EC
ED <- d$ED
d <- d$ydata
K <- readRDS("~/repos/kmmr/dev/weights_K.rds")
Ki <- solve(K)
attr(Ki, "INVERSE") <- TRUE

d <- droplevels(d[d$Man == "LN" & d$Env == "HAN_2015",])

library(asreml)

fit.asr <- asreml(fixed = GY ~ 1,
                  random = ~ vm(Gen, Ki),
                  residual = ~ units,
                  data = d)

fit.asr.wt <- asreml(fixed = GY ~ 1,
                     random = ~ vm(Gen, Ki),
                     residual = ~ units,
                     data = d,
                     weights = weight,
                     family = asr_gaussian(dispersion = 1.0))

fit.kmmr <- kmmr(fixed = GY ~ 1,
                 random = ~ re(~ lk(Gen, Ki), vars = ~ 1),
                 residual = ~ units,
                 data = d)

fit.kmmr.wt <- kmmr(fixed = GY ~ 1,
                    random = ~ re(~ lk(Gen, Ki), vars = ~ 1),
                    residual = ~ units,
                    data = d,
                    weights = "weight",
                    init = c(10.8364028, 1.0),
                    isFixed = c(FALSE, TRUE))





