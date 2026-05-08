library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = TRUE, score = TRUE)

# Seed:
set.seed(1997)

# Loading example data (DROPS):
d <- readRDS("dev/exampledata.rds")[, c("Gen", "Man", "Env", "GY")]
ED <- readRDS("dev/exampledata_ED.rds")
EC <- readRDS("dev/exampledata_EC.rds")
K <- readRDS("dev/exampledata_K.rds")

d <- droplevels(d[d$Env %in% c("Bol12", "Kar12", "Cam12") & d$Man == "R",])
ED <- ED[levels(d$Env), levels(d$Env)]
EC <- EC[levels(d$Env), levels(d$Env)]
K <- K[levels(d$Gen), levels(d$Gen)]

ECi <- solve(EC)
attr(ECi, "INVERSE") <- TRUE
Ki <- solve(K)
attr(Ki, "INVERSE") <- TRUE

d$EnvGen <- factor(paste(as.character(d$Env), as.character(d$Gen), sep = ":"), paste(as.character(d$Env), as.character(d$Gen), sep = ":"), paste(as.character(d$Env), as.character(d$Gen), sep = ":"))
ECK <- kronecker(EC, K, make.dimnames = T)

# These two are the same, as expected:
fit.asr <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ corv(Env):vm(Gen, K) + vm(EnvGen, ECK),
                  residual = ~ idv(units),
                  data = d)

fit.kmmr.stack <- kmmr(fixed = GY ~ -1 + Env,
                       random = ~ re(~ stack(~ cors(Env) + lk(Env, EC)):lk(Gen, Ki), vars = ~ 1),
                       residual = ~ units,
                       data = d)#,
                       # init = fit.asr$trace[c(6, 7, 5, 8), 1])

plot(fit.asr$coefficients$random[1:738] + fit.asr$coefficients$random[739:1476], fit.kmmr.stack$u)

# But this one is different... Why?
# Fixed! It was a bug with the wrapper.
fit.kmmr <- kmmr(fixed = GY ~ -1 + Env,
                 random = ~ re(~ cors(Env):lk(Gen, Ki), vars = ~ 1) + re(~ lk(Env, ECi):lk(Gen, Ki), vars = ~ 1),
                 residual = ~ units,
                 data = d,
                 init = fit.asr$trace[c(6, 5, 7, 8), 1])

plot(fit.asr$coefficients$random[1:738] + fit.asr$coefficients$random[739:1476], fit.kmmr$u[1:738] + fit.kmmr$u[739:1476])
plot(fit.asr$coefficients$random, fit.kmmr$u)
