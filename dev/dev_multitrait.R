d <- readRDS("dev/exampledata.rds")[, c("Env", "Man", "Gen", "GY")]
d <- droplevels(d[d$Env == "Bol12",])
K <- readRDS("dev/exampledata_K.rds")
Ki.asr <- solve(K); attr(Ki.asr, "INVERSE") <- TRUE
Ki <- Matrix(solve(K))
# library(asreml)
# 
# microbenchmark::microbenchmark(fit.kmmr <- kmmr(fixed = GY ~ 1 + Man,
#                                                 random = ~ re(~ corg(Man):lk(Gen, Ki), vars = ~ Man),
#                                                 residual = ~ re(~ corg(Man):id(Gen), vars = ~ Man),
#                                                 data = d),
#                                fit.asr <- asreml(fixed = GY ~ 1 + Man,
#                                                  random = ~ corgh(Man):vm(Gen, Ki.asr),
#                                                  residual = ~ corgh(Man):id(units),
#                                                  asmv = Man,
#                                                  data = d),
#                                times = 20)
# 
# plot(fit.kmmr$u@x, fit.asr$coefficients$random[,])
# cor(fit.kmmr$u@x, fit.asr$coefficients$random[,])

# Checking scrambled data ====
set.seed(1997)
d2 <- d[sample(1:nrow(d), nrow(d)),]
# Ordering does not matter for kmmr:
fit.kmmr1 <- kmmr(fixed = GY ~ 1 + Man,
                  random = ~ re(~ corg(Man):lk(Gen, Ki), vars = ~ Man),
                  residual = ~ re(~ corg(Man):id(Gen), vars = ~ Man),
                  data = d)

fit.kmmr2 <- kmmr(fixed = GY ~ 1 + Man,
                  random = ~ re(~ corg(Man):lk(Gen, Ki), vars = ~ Man),
                  residual = ~ re(~ corg(Man):id(Gen), vars = ~ Man),
                  data = d2)

all.equal(fit.kmmr1$u, fit.kmmr2$u)
all.equal(fit.kmmr1$k, fit.kmmr2$k)
all.equal(fit.kmmr1$k[,1], fit.kmmr2$k[,1])

fit.asr <- asreml(fixed = GY ~ 1 + Man,
                  random = ~ corgh(Man):vm(Gen, Ki.asr),
                  residual = ~ corgh(Man):id(units),
                  asmv = Man,
                  data = d)

plot(fit.kmmr2$u@x, fit.asr$coefficients$random[,])
cor(fit.kmmr2$u@x, fit.asr$coefficients$random[,])

# Checking missing values ====
set.seed(1997)
setNA <- sample(1:nrow(d), nrow(d) * 0.25)
d2 <- d
d2$GY[setNA] <- NA

fit.kmmr <- kmmr(fixed = GY ~ 1 + Man,
                 random = ~ re(~ corg(Man):lk(Gen, Ki), vars = ~ Man),
                 residual = ~ re(~ corg(Man):id(Gen), vars = ~ Man),
                 data = d2)

fit.asr <- asreml(fixed = GY ~ 1 + Man,
                  random = ~ corgh(Man):vm(Gen, Ki.asr),
                  residual = ~ corgh(Man):id(units),
                  asmv = Man,
                  data = d2)

plot(fit.kmmr$u@x, fit.asr$coefficients$random[,])
cor(fit.kmmr$u@x, fit.asr$coefficients$random[,])

# Checking scrambled data and missing values ====
set.seed(1997)
d2 <- d2[sample(1:nrow(d2), nrow(d2)),]
fit.kmmr <- kmmr(fixed = GY ~ 1 + Man,
                 random = ~ re(~ corg(Man):lk(Gen, Ki), vars = ~ Man),
                 residual = ~ re(~ corg(Man):id(Gen), vars = ~ Man),
                 data = d2)

plot(fit.kmmr$u@x, fit.asr$coefficients$random[,])
cor(fit.kmmr$u@x, fit.asr$coefficients$random[,])

# Checking GxE data ====
rm(d, d2, fit.asr, fit.kmmr, setNA)
d <- readRDS("dev/exampledata.rds")[, c("Env", "Man", "Gen", "GY")]
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY")])

fit.kmmr <- kmmr(fixed = GY ~ -1 + Env,
                 random = ~ re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
                 residual = ~ units,
                 data = d,
                 sparse = TRUE)

fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ fa(Env, 3):vm(Gen, Ki.asr),
                  residual = ~ units,
                  data = d,
                  maxit = 20)

plot(fit.kmmr$u@x[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d),])
cor(fit.kmmr$u@x[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d),])

# Checking scrambled GxE data ====
set.seed(1997)
d2 <- d[sample(1:nrow(d), nrow(d)),]

fit.kmmr2 <- kmmr(fixed = GY ~ -1 + Env,
                 random = ~ re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
                 residual = ~ units,
                 data = d2,
                 sparse = TRUE)

fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ fa(Env, 3):vm(Gen, Ki.asr),
                  residual = ~ units,
                  data = d2,
                  maxit = 20)

plot(fit.kmmr2$u@x[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d),])
cor(fit.kmmr2$u@x[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d),])

# Checking scrambled GxE data with missing values ====
set.seed(1997)
d2$GY[sample(1:nrow(d), nrow(d) * 0.25)] <- NA

fit.kmmr2 <- kmmr(fixed = GY ~ -1 + Env,
                  random = ~ re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
                  residual = ~ units,
                  data = d2,
                  sparse = TRUE)

fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ fa(Env, 3):vm(Gen, Ki.asr),
                  residual = ~ units,
                  data = d2,
                  maxit = 25)

plot(fit.kmmr2$u@x[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d),])
cor(fit.kmmr2$u@x[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d),])

# Showing the equivalent of the shorthand residual = ~ units argument:
fit.kmmr3 <- kmmr(fixed = GY ~ -1 + Env,
                  random = ~ re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
                  residual = ~ re(~ id(Env):id(Gen), vars = ~ 1),
                  data = d2,
                  sparse = TRUE)

all.equal(fit.kmmr2$u, fit.kmmr3$u)
all.equal(fit.kmmr2$k, fit.kmmr3$k)

# Checking GxE data with weights ====
rm(list = ls())
# Weights in this dataframe are the reciprocals of variances, that's what asreml needs.
# kmmr needs the variances themselves.
d <- readRDS("dev/exampledata.weights.rds")[, c("Env", "Man", "Gen", "GY", "weight")]
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY", "weight")])
d$weight.kmmr <- 1 / d$weight
K <- readRDS("dev/exampledata_K.rds")[levels(d$Gen), levels(d$Gen)]
Ki.asr <- solve(K); attr(Ki.asr, "INVERSE") <- TRUE
Ki <- Matrix(solve(K))

fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ fa(Env, 3):vm(Gen, Ki.asr),
                  residual = ~ units,
                  data = d,
                  maxit = 25,
                  family = asr_gaussian(dispersion = 1.0),
                  weights = weight)

# First we get starting values:
init <- kmmr(fixed = GY ~ -1 + Env,
             random = ~ re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
             residual = ~ units,
             data = d,
             sparse = TRUE,
             getStartingValues = TRUE)

# Then we fix the residual variance at 1:
init$Init[init$Parameter == "residual(vars{1})"] <- 1.0
init$Fixed[init$Parameter == "residual(vars{1})"] <- TRUE

# Then we fit the model and supply weights:
fit.kmmr <- kmmr(fixed = GY ~ -1 + Env,
                 random = ~ re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
                 residual = ~ units,
                 data = d,
                 sparse = TRUE,
                 init = init$Init,
                 isFixed = init$Fixed,
                 weights = "weight.kmmr")

plot(fit.kmmr$u@x[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d),])
cor(fit.kmmr$u@x[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d),])

# Checking scrambled GxE data with weights and missing values ====
rm(list = ls())
# Weights in this dataframe are the reciprocals of variances, that's what asreml needs.
# kmmr needs the variances themselves.
d <- readRDS("dev/exampledata.weights.rds")[, c("Env", "Man", "Gen", "GY", "weight")]
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY", "weight")])
d$weight.kmmr <- 1 / d$weight
K <- readRDS("dev/exampledata_K.rds")[levels(d$Gen), levels(d$Gen)]
Ki.asr <- solve(K); attr(Ki.asr, "INVERSE") <- TRUE
Ki <- Matrix(solve(K))

# Interesting, at seed 1997 asreml gets a lower LL, at seed 2005 kmmr gets a lower LL...
set.seed(2005)
d2 <- d[sample(1:nrow(d), nrow(d)),]
d2$GY[sample(1:nrow(d), nrow(d) * 0.25)] <- NA

fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ fa(Env, 3):vm(Gen, Ki.asr),
                  residual = ~ units,
                  data = d2,
                  maxit = 50,
                  family = asr_gaussian(dispersion = 1.0),
                  weights = weight)

# First we get starting values:
init <- kmmr(fixed = GY ~ -1 + Env,
             random = ~ re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
             residual = ~ units,
             data = d2,
             sparse = TRUE,
             getStartingValues = TRUE)

# Then we fix the residual variance at 1:
init$Init[init$Parameter == "residual(vars{1})"] <- 1.0
init$Fixed[init$Parameter == "residual(vars{1})"] <- TRUE

# Then we fit the model and supply weights:
fit.kmmr <- kmmr(fixed = GY ~ -1 + Env,
                 random = ~ re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
                 residual = ~ units,
                 data = d2,
                 sparse = TRUE,
                 init = init$Init,
                 isFixed = init$Fixed,
                 weights = "weight.kmmr",
                 maxit = 50)

plot(fit.kmmr$u@x[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d),])
cor(fit.kmmr$u@x[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d),])

# quite a large difference in LL and BLUPs. Is this purely due to starting values?
init.asr <- fit.asr$trace[5:nrow(fit.asr$trace), 1]
L.init.asr <- matrix(init.asr[15:56], 14, 3)
PSI.init.asr <- diag(init.asr[1:14])
covmat.init.asr <- L.init.asr %*% t(L.init.asr) + PSI.init.asr
vars.init.asr <- diag(covmat.init.asr)
L.init.asr <- (1 / sqrt(vars.init.asr)) * L.init.asr
PSI.init.asr <- (1 / vars.init.asr) * PSI.init.asr
all.equal(L.init.asr %*% t(L.init.asr) + PSI.init.asr, cov2cor(covmat.init.asr)) # TRUE

# And fit the model again:
fit.kmmr <- kmmr(fixed = GY ~ -1 + Env,
                 random = ~ re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
                 residual = ~ units,
                 data = d2,
                 sparse = TRUE,
                 init = c(vars.init.asr, as.numeric(L.init.asr), 1.0), # Replace our initial values with those from asreml
                 isFixed = init$Fixed,
                 weights = "weight.kmmr")

# Doesn't look good, but the LL is the same, so this seems to purely be an issue with constraining updates:
all.equal(fit.asr$trace["LogLik", 1], fit.kmmr$LL[1]) # TRUE

# Try with asremls final estimates:
init.asr <- fit.asr$trace[5:nrow(fit.asr$trace), ncol(fit.asr$trace)]
L.init.asr <- matrix(init.asr[15:56], 14, 3)
PSI.init.asr <- diag(init.asr[1:14])
covmat.init.asr <- L.init.asr %*% t(L.init.asr) + PSI.init.asr
vars.init.asr <- diag(covmat.init.asr)
L.init.asr <- (1 / sqrt(vars.init.asr)) * L.init.asr
PSI.init.asr <- (1 / vars.init.asr) * PSI.init.asr
all.equal(L.init.asr %*% t(L.init.asr) + PSI.init.asr, cov2cor(covmat.init.asr)) # TRUE

# Strange that the LL goes down, then back up a litle again...
fit.kmmr <- kmmr(fixed = GY ~ -1 + Env,
                 random = ~ re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
                 residual = ~ units,
                 data = d2,
                 sparse = TRUE,
                 init = c(vars.init.asr, as.numeric(0.9999 * L.init.asr), 1.0), # We need the 0.9999 because asreml's fit has a single communality of 1.
                 isFixed = init$Fixed,
                 weights = "weight.kmmr",
                 maxit = 25)

# We do get essentially the same fit now, so it must be a case where the constraint rules of kmmr
# are not ideal.
plot(fit.kmmr$u@x[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d),])
cor(fit.kmmr$u@x[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d),])
