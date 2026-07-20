d <- readRDS("dev/exampledata.rds")
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY")])
K <- readRDS("dev/exampledata_K.rds")
Ki <- solve(K)
ED <- readRDS("dev/exampledata_ED.rds")
EC <- readRDS("dev/exampledata_EC.rds")

fit <- kmmr(fixed = GY ~ -1 + Env,
            random = ~ re(~ fa(Env, 3):lk(Gen, Ki), vars = ~ Env),
            residual = ~ units,
            data = d)

