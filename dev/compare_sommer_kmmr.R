tolLL <- 5e-5
d <- readRDS("~/repos/kmmr/dev/exampledata.rds")
K <- readRDS("~/repos/kmmr/dev/exampledata_K.rds")
Ki <- Matrix(solve(K))
attr(Ki, "inverse") <- TRUE
system.time(
  fit.kmr <- kmmr(fixed = GY ~ -1 + ManEnv,
                  random = ~ re(~ exfa(ManEnv, 3):lk(Gen, Ki), vars = ~ ManEnv),
                  residual = ~ units,
                  data = d, maxit = 150, timing = TRUE, tolConvLL = tolLL)
) # 35.1 seconds, 28 iters, 1.25 s/iter

system.time(
  fit.smr <- mmes(fixed = GY ~ -1 + ManEnv,
                  random = ~ vsm(fam(ManEnv, 3), ism(Gen), Gu = Ki),
                  rcov = ~ units,
                  data = d, nIters = 150, tolParConvLL = tolLL)
) # 101.9 seconds, 87 iters, 1.17 s/iter

plot(fit.kmr$u@x[1:6888],
     fit.smr$u[,]); abline(0, 1, col = "red")
