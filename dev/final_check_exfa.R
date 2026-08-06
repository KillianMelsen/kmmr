d <- readRDS("dev/exampledata.rds")
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY")])
Ki <- solve(readRDS("dev/exampledata_K.rds"))
KiMatrix <- Matrix(Ki)
attr(Ki, "INVERSE") <- TRUE

library(asreml)
asreml.options(workspace = "10000mb", pworkspace = "10000mb")

system.time({
fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ fa(Env, 3):vm(Gen, Ki),
                  residual = ~ idv(units),
                  data = d,
                  maxit = 50)
}) # 74.580 seconds, 34 iterations, 2.193529 s/iter, LL = -684.787

system.time({
fit.kmr <- kmmr(fixed = GY ~ -1 + Env,
                random = ~ re(~ exfa(Env, 3):lk(Gen, KiMatrix), vars = ~ Env),
                residual = ~ units,
                data = d,
                maxit = 25,
                sparse = TRUE,
                super = TRUE)
}) # 11.284 seconds, 13 iterations, 0.868 s/iter, LL = -684.7136

plot(fit.asr$coefficients$random[1:3444], fit.kmr$u@x[1:3444])
cor(fit.asr$coefficients$random[1:3444], fit.kmr$u@x[1:3444], method = "pearson")
cor(fit.asr$coefficients$random[1:3444], fit.kmr$u@x[1:3444], method = "spearman")

plot(fit.asr$coefficients$random[3445:4182], fit.kmr$u@x[3445:4182]) # Factor scores


