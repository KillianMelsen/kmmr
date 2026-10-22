suppressMessages(library(asreml))
suppressMessages(library(kmmr))
suppressMessages(library(Matrix))
asreml.options(workspace = "5000mb")

d <- readRDS("~/repos/kmmr-benchmarking/data/DROPS_weights.rds")
Ki.asr <- solve(readRDS("~/repos/kmmr-benchmarking/data/DROPS_K.rds")[levels(d$Gen), levels(d$Gen)])
Ki.smr <- Ki.kmr <- Matrix(Ki.asr)
attr(Ki.asr, "INVERSE") <- TRUE
attr(Ki.smr, "inverse") <- TRUE

fit.asr <- asreml(fixed = GY ~ -1 + ManEnv,
                  random = ~ fa(ManEnv, 1):vm(Gen, Ki.asr),
                  residual = ~ units,
                  weights = weight, family = asr_gaussian(dispersion = 1.0),
                  data = d, maxit = 150)
tolLL <- fit.asr$trace[1, ncol(fit.asr$trace) - 1] - fit.asr$trace[1, ncol(fit.asr$trace) - 2]

# mixed Cmode ====
sv <- kmmr(fixed = GY ~ -1 + ManEnv,
           random = ~ re(~ exfa(ManEnv, 1):lk(Gen, Ki.kmr), vars = ~ ManEnv),
           residual = ~ units,
           weights = "weight",
           data = d, maxit = 150, getStartingValues = TRUE)
sv$Init[nrow(sv)] <- 1.0
sv$Fixed[nrow(sv)] <- TRUE
fit.kmr <- kmmr(fixed = GY ~ -1 + ManEnv,
                random = ~ re(~ exfa(ManEnv, 1):lk(Gen, Ki.kmr), vars = ~ ManEnv),
                residual = ~ units,
                weights = "weight",
                data = d, maxit = 150,
                init = sv$Init, isFixed = sv$Fixed, tolConvLL = tolLL)

system.time(preds.asr <- predict(fit.asr, classify = "ManEnv:Gen", maxit = 1)$pvals)
k <- matrix(0, nrow(fit.kmr$b) + nrow(fit.kmr$u), 28 * 234)
for (i in 1:28) {
  k[i, ((i - 1) * 234 + 1):(i * 234)] <- 1 # Fixed effect for each ManEnv.
}
k[29:(ncol(k) + 28), 1:ncol(k)] <- diag(1, 28 * 234)
k <- Matrix(k)
system.time(preds.kmr <- predict(fit.kmr, predictionMatrix = k, chunks = TRUE, chunkSize = 500))
plot(preds.asr$predicted.value, preds.kmr$Prediction); abline(0, 1, col = "red")
plot(preds.asr$std.error, preds.kmr$SE); abline(0, 1, col = "red")
rm(preds.kmr)

# sparse Cmode ====
sv <- kmmr(fixed = GY ~ -1 + ManEnv,
           random = ~ re(~ exfa(ManEnv, 1):lk(Gen, Ki.kmr), vars = ~ ManEnv),
           residual = ~ units,
           weights = "weight",
           data = d, maxit = 150, getStartingValues = TRUE)
sv$Init[nrow(sv)] <- 1.0
sv$Fixed[nrow(sv)] <- TRUE
fit.kmr <- kmmr(fixed = GY ~ -1 + ManEnv,
                random = ~ re(~ exfa(ManEnv, 1):lk(Gen, Ki.kmr), vars = ~ ManEnv),
                residual = ~ units,
                weights = "weight",
                data = d, maxit = 150,
                init = sv$Init, isFixed = sv$Fixed, tolConvLL = tolLL, denseThreshold = 0.1)

k <- matrix(0, nrow(fit.kmr$b) + nrow(fit.kmr$u), 28 * 234)
for (i in 1:28) {
  k[i, ((i - 1) * 234 + 1):(i * 234)] <- 1 # Fixed effect for each ManEnv.
}
k[29:(ncol(k) + 28), 1:ncol(k)] <- diag(1, 28 * 234)
k <- Matrix(k)
system.time(preds.kmr <- predict(fit.kmr, predictionMatrix = k, chunks = TRUE, chunkSize = 500))
plot(preds.asr$predicted.value, preds.kmr$Prediction); abline(0, 1, col = "red")
plot(preds.asr$std.error, preds.kmr$SE); abline(0, 1, col = "red")
rm(preds.kmr)

# dense Cmode ====
sv <- kmmr(fixed = GY ~ -1 + ManEnv,
           random = ~ re(~ fa(ManEnv, 1):lk(Gen, Ki.kmr), vars = ~ ManEnv),
           residual = ~ units,
           weights = "weight",
           data = d, maxit = 150, getStartingValues = TRUE)
sv$Init[nrow(sv)] <- 1.0
sv$Fixed[nrow(sv)] <- TRUE
fit.kmr <- kmmr(fixed = GY ~ -1 + ManEnv,
                random = ~ re(~ fa(ManEnv, 1):lk(Gen, Ki.kmr), vars = ~ ManEnv),
                residual = ~ units,
                weights = "weight",
                data = d, maxit = 150,
                init = sv$Init, isFixed = sv$Fixed, tolConvLL = tolLL, denseThreshold = 0.0)

k <- matrix(0, nrow(fit.kmr$b) + nrow(fit.kmr$u), 28 * 234)
for (i in 1:28) {
  k[i, ((i - 1) * 234 + 1):(i * 234)] <- 1 # Fixed effect for each ManEnv.
}
k[29:(ncol(k) + 28), 1:ncol(k)] <- diag(1, 28 * 234)
k <- Matrix(k)
system.time(preds.kmr <- predict(fit.kmr, predictionMatrix = k, chunks = TRUE, chunkSize = 500))
plot(preds.asr$predicted.value, preds.kmr$Prediction); abline(0, 1, col = "red")
plot(preds.asr$std.error, preds.kmr$SE); abline(0, 1, col = "red")



