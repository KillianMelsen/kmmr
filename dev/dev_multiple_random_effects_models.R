library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = TRUE, score = TRUE)

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

ds <- droplevels(d[d$Env %in% c("Bol12", "Kar12", "Cam12", "Deb12"),])
ds$ManGen <- as.factor(paste(as.character(ds$Man), as.character(ds$Gen), sep = ":"))
combis <- c(paste(levels(ds$Env)[1], levels(ds$Gen), sep = ":"),
            paste(levels(ds$Env)[2], levels(ds$Gen), sep = ":"),
            paste(levels(ds$Env)[3], levels(ds$Gen), sep = ":"),
            paste(levels(ds$Env)[4], levels(ds$Gen), sep = ":"))
ds$EnvGen <- factor(combis, combis, combis)

fit <- asreml(fixed = GY ~ 1,
              random = ~ corgv(Man):vm(Gen, Ki) + idh(Env):vm(Gen, Ki),
              residual = ~ idv(units),
              data = ds, maxit = 20)

# summary(fit)$varcomp
# fit$trace
init <- fit$trace[c(6, 5, 7:10, 11), 1]
fit.kmmr <- kmmr(fixed = GY ~ 1,
                 random = ~ re(~ corg(Man):lk(Gen, Ki), vars = ~ 1) + re(~ id(Env):lk(Gen, Ki), vars = ~ Env),
                 residual = ~ units,
                 data = ds,
                 maxit = 20,
                 # init = as.numeric(init),
                 tolConvLL = 0.0001)



