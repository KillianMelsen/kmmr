Ki <- solve(Matrix(readRDS("~/repos/kmmr/dev/exampledata_K.rds")))
d <- readRDS("~/repos/kmmr/dev/exampledata.rds")
# d <- droplevels(d[d$Man == "R" & d$Env %in% c("Cam12", "Kar12"),])

fit <- kmmr(fixed = GY ~ 1 + ManEnv,
            random = ~ re(~ fa(ManEnv, 3):lk(Gen, Ki), vars = ~ ManEnv),
            residual = ~ units,
            data = d,
            sparse = FALSE,
            useGPU = TRUE)

