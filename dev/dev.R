# Seed:
set.seed(1997)

# Loading example data (DROPS):
d <- readRDS("dev/exampledata.rds")[, c("Gen", "Man", "Env", "GY")]
ED <- readRDS("dev/exampledata_ED.rds")
EC <- readRDS("dev/exampledata_EC.rds")
GD <- readRDS("dev/exampledata_GD.rds")
K <- readRDS("dev/exampledata_K.rds")
Ki <- solve(K)
attr(Ki, "INVERSE") <- TRUE

# Adding a ManEnv factor for ASReml/cornfruit:
d$ManEnv <- factor(paste(as.character(d$Man), as.character(d$Env), sep = ":"),
                   levels = c(paste0("R:", levels(d$Env)),
                              paste0("W:", levels(d$Env))))

# corg:cors:lk ====
ds <- droplevels(d[d$Env %in% c("Bol12", "Cam12", "Kar12"),])
fit.kmmr <- kmmr(fixed = GY ~ -1 + Man:Env,
                 random = ~ re(~ corg(Man):corg(Env):lk(Gen, Ki), vars = ~ Man),
                 residual = ~ units,
                 data = ds,
                 tolConvLL = 0.0001)

summary(fit.kmmr, round = 3)
# summary(fit.kmmr)$varcomp
# summary(fit.kmmr)$logLik

# Checking match with ASReml/cornfruit:
library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb")
fit.asr <- asreml(fixed = GY ~ -1 + ManEnv,
                  random = ~ corgh(Man):corg(Env):vm(Gen, Ki),
                  residual = ~ idv(units),
                  data = ds)
summary(fit.asr)$varcomp

# cors:lk ====
ds <- droplevels(ds[ds$Man == "R",])
fit.kmmr <- kmmr(fixed = GY ~ -1 + Env,
                 random = ~ re(~ cors(Env):lk(Gen, Ki), vars = ~ 1),
                 residual = ~ units,
                 data = ds,
                 tolConvLL = 0.01)

# Checking match with ASReml/cornfruit:
library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb")
fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ corv(Env):vm(Gen, Ki),
                  residual = ~ idv(units),
                  data = ds)

# Single variance Gaussian kernel ====
fit.kmmr <- kmmr(fixed = GY ~ -1 + Man:Env,
                 random = ~ re(~ corg(Man):gk(Env, base = ED):lk(Gen, base = Ki), vars = ~ Man),
                 residual = ~ units,
                 data = d,
                 tolConvLL = 0.01)

# Checking match with ASReml/cornfruit:
library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb")
vf <- cornfruit::svgk(D = ED[levels(d$Env), levels(d$Env)])
init <- rep(0.1, 2 + 1 + 1)
type <- c("V", "V", "R", "V")
cons <- c("P", "P", "U", "P")
fit.asr <- asreml(fixed = GY ~ -1 + ManEnv,
                  random = ~ own(ManEnv, "vf", init, type, cons):vm(Gen, Ki),
                  residual = ~ idv(units),
                  data = d)

summary(fit.asr)$varcomp[1:5, "component"] - fit.kmmr$k[, ncol(fit.kmmr$k)]
fit.asr$loglik - fit.kmmr$LL[length(fit.kmmr$LL)]

# Multiple variance Gaussian kernel ====
fit.kmmr <- kmmr(fixed = GY ~ -1 + Man:Env,
                 random = ~ re(~ corg(Man):gk(Env, base = ED):lk(Gen, base = Ki), vars = ~ Man:Env),
                 residual = ~ units,
                 data = d,
                 tolConvLL = 0.01)

# Checking match with ASReml/cornfruit:
library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb")
vf <- cornfruit::mvgk(D = ED[levels(d$Env), levels(d$Env)])
init <- rep(0.1, 28 + 1 + 1)
type <- c(rep("V", 28), "R", "V")
cons <- c(rep("P", 28), "U", "P")
fit.asr <- asreml(fixed = GY ~ -1 + ManEnv,
                  random = ~ own(ManEnv, "vf", init, type, cons):vm(Gen, Ki),
                  residual = ~ idv(units),
                  data = d)

range(summary(fit.asr)$varcomp[1:31, "component"] - fit.kmmr$k[, ncol(fit.kmmr$k)])
fit.asr$loglik - fit.kmmr$LL[length(fit.kmmr$LL)]

# Multiple variance linear kernel ====
fit.kmmr <- kmmr(fixed = GY ~ -1 + Man:Env,
                 random = ~ re(~ corg(Man):lk(Env, base = EC):lk(Gen, base = Ki), vars = ~ Man:Env),
                 residual = ~ units,
                 data = d,
                 tolConvLL = 0.01)

# Checking match with ASReml/cornfruit:
library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb")
vf <- cornfruit::mvlk(C = EC[levels(d$Env), levels(d$Env)])
init <- rep(0.1, 28 + 1)
type <- c(rep("V", 28), "R")
cons <- c(rep("P", 28), "U")
fit.asr <- asreml(fixed = GY ~ -1 + ManEnv,
                  random = ~ own(ManEnv, "vf", init, type, cons):vm(Gen, Ki),
                  residual = ~ idv(units),
                  data = d)

range(summary(fit.asr)$varcomp[1:30, "component"] - fit.kmmr$k[, ncol(fit.kmmr$k)])
fit.asr$loglik - fit.kmmr$LL[length(fit.kmmr$LL)]




# Gaussian kernel GBLUP ====
# This one looks fine
ds <- droplevels(d[d$Man == "R" & d$Env == "Bol12",])
fit.kmmr <- kmmr(fixed = GY ~ 1,
                 random = ~ re(~ gk(Gen, base = GD), vars = ~ 1),
                 residual = ~ units,
                 data = ds,
                 tolConvLL = 0.01)

# This one converges but takes a very strange path, needs fixing of the max changes in variances
ds <- droplevels(d[d$Man == "R" & d$Env == "Kar12",])
fit.kmmr <- kmmr(fixed = GY ~ 1,
                 random = ~ re(~ gk(Gen, base = GD), vars = ~ 1),
                 residual = ~ units,
                 data = ds,
                 tolConvLL = 0.01, maxit = 100)


