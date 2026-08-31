Ki <- solve(Matrix(readRDS("~/repos/kmmr/dev/exampledata_K.rds")))
Ki.asr <- as.matrix(Ki)
attr(Ki.asr, "INVERSE") <- TRUE
d <- readRDS("~/repos/kmmr/dev/exampledata.rds")
# d <- droplevels(d[d$Man == "R" & d$Env %in% c("Cam12", "Kar12"),])

fit.CPU <- kmmr(fixed = GY ~ 1 + ManEnv,
                random = ~ re(~ exfa(ManEnv, 3):lk(Gen, Ki), vars = ~ ManEnv),
                residual = ~ units,
                data = d,
                sparse = TRUE,
                useGPU = FALSE)

fit.GPU <- kmmr(fixed = GY ~ 1 + ManEnv,
                random = ~ re(~ fa(ManEnv, 3):lk(Gen, Ki), vars = ~ ManEnv),
                residual = ~ units,
                data = d,
                sparse = FALSE,
                useGPU = TRUE)

library(asreml)
asreml.options(workspace = "8000mb")
fit.asr <- asreml(fixed = GY ~ 1 + ManEnv,
                  random = ~ fa(ManEnv, 3):vm(Gen, Ki.asr),
                  residual = ~ units,
                  data = d,
                  maxit = 20)
