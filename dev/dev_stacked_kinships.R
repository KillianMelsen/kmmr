library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = TRUE, score = TRUE)

# Seed:
set.seed(1997)

# Loading example data (DROPS):
d <- readRDS("dev/exampledata.rds")[, c("Gen", "Man", "Env", "GY")]
d <- droplevels(d[d$Man == "R" & d$Env == "Bol12",])
GD <- readRDS("dev/exampledata_GD.rds")[levels(d$Gen), levels(d$Gen)]
K1 <- readRDS("dev/exampledata_K.rds")[levels(d$Gen), levels(d$Gen)]
K2 <- exp(-0.1 * GD)
K1i <- solve(K1)
K2i <- solve(K2)
attr(K1i, "INVERSE") <- TRUE
attr(K2i, "INVERSE") <- TRUE

fit.asr <- asreml(fixed = GY ~ 1,
                  random = ~ vm(Gen, K1) + vm(Gen, K2),
                  residual = ~ idv(units),
                  data = d)

# This doesn't work for whatever reason. It should work...
# vf <- function(order, kappa) {
#   V <- kappa[1] * K1 + kappa[2] * K2
#   return(list(V, K1, K2))
# }
# init <- c(0.2, 0.2)
# type <- c("V", "V")
# cons <- c("P", "P")
# fit2 <- asreml(fixed = GY ~ 1,
#                random = ~ own(Gen, "vf", init, type, cons),
#                residual = ~ idv(units),
#                data = d)

summary(fit.asr)$varcomp

fit <- kmmr(fixed = GY ~ 1,
            random = ~ re(~ stack(~ lk(Gen, K1i) + lk(Gen, K2i)), vars = ~ 1),
            residual = ~ units,
            data = d)


