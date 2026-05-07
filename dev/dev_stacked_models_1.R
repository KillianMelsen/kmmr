library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = TRUE, score = TRUE)

# Seed:
set.seed(1997)

# Loading example data (DROPS):
d <- readRDS("dev/exampledata.rds")[, c("Gen", "Man", "Env", "GY")]
GD <- readRDS("dev/exampledata_GD.rds")
K <- readRDS("dev/exampledata_K.rds")
Ki <- solve(K)
attr(Ki, "INVERSE") <- TRUE

d <- droplevels(d[d$Man == "R" & d$Env == "Bol12",])
K <- K[levels(d$Gen), levels(d$Gen)]
GD <- GD[levels(d$Gen), levels(d$Gen)]
Ki <- Ki[levels(d$Gen), levels(d$Gen)]

fit1 <- kmmr(fixed = GY ~ 1,
             random = ~ re(~ stack(~ gk(Gen, GD) + lk(Gen, K)), vars = ~ 1),
             residual = ~ units,
             data = d,
             tolConvLL = 1e-9)

fit2 <- kmmr(fixed = GY ~ 1,
             random = ~ re(~ stack(~ lk(Gen, K) + gk(Gen, GD)), vars = ~ 1),
             residual = ~ units,
             data = d,
             tolConvLL = 1e-9)

summary(fit1)
summary(fit2)

fit3 <- kmmr(fixed = GY ~ 1,
             random = ~ re(~ gk(Gen, GD), vars = ~ 1) + re(~ lk(Gen, Ki), vars = ~ 1),
             residual = ~ units,
             data = d,
             tolConvLL = 1e-9)

summary(fit1)
summary(fit2)
summary(fit3)
