library(asreml)
asreml.options(workspace = "8000mb", pworkspace = "8000mb")

Ki.asr <- solve(readRDS("~/repos/kmmr/dev/exampledata_K.rds"))
Ki.kmr <- Matrix(Ki.asr)
attr(Ki.asr, "INVERSE") <- TRUE
d <- readRDS("~/repos/kmmr/dev/exampledata.rds")
d <- droplevels(d[d$Man == "R",])

# Fitting models ====
system.time(
fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ facv(Env, 3):vm(Gen, Ki.asr),
                  residual = ~ units,
                  data = d,
                  maxit = 30)
)

system.time(
fit.kmr <- kmmr(fixed = GY ~ -1 + Env,
                random = ~ re(~ fa(Env, 3):lk(Gen, Ki.kmr), vars = ~ Env),
                residual = ~ units,
                data = d,
                maxit = 30, denseThreshold = 0.00001)
)

# Predicting ====
p.asr <- predict(fit.asr, classify = "Env:Gen")$pvals

pv <- Matrix(0, 14 + 3444, 6, sparse = FALSE)
pv[1, 1:6] <- 1 # BLUE for Gai12
pv[cbind(c(15:20), c(1:6))] <- 1 # BLUPs for the first 6 genotypes in Gai12.
pv <- as(pv, "sparseMatrix")
p.kmr <- predict(fit.kmr, predictionMatrix = pv)

head(p.asr)
p.kmr


