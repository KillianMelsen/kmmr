# Implement the ar1 correlation model
K <- readRDS("dev/exampledata_K.rds")
Ki.asr <- solve(K); attr(Ki.asr, "INVERSE") <- TRUE
Ki.kmmr <- solve(Matrix(K))
d <- readRDS("~/repos/kmmr/dev/exampledata_singleStage.rds")
d <- droplevels(d[d$type == "CrossValidation",])
d <- droplevels(d[d$Variety_ID %in% rownames(K),])
d <- droplevels(d[d$Experiment != "Deb13R",])
d$Env <- substr(d$Experiment, 1, 5)
d$Man <- substr(d$Experiment, 6, 6)
d$Gen <- d$Variety_ID
d$GY <- d$grain.yield
d$Rep <- d$Replicate
d$Block <- d$block
d$Plot <- d$plotID
d$Row <- d$Row
d$Col <- d$Column

d <- d[, c("Plot", "Row", "Col", "Block", "Rep", "Man", "Env", "Gen", "GY")]
d <- na.omit(d)
d$Rep <- factor(d$Rep, sort(unique(d$Rep)), sort(unique(d$Rep)))
d$Row <- factor(d$Row, sort(unique(d$Row)), sort(unique(d$Row))) # Make sure it's converted to a factor correctly!
d$Col <- factor(d$Col, sort(unique(d$Col)), sort(unique(d$Col))) # Make sure it's converted to a factor correctly!
d$Block <- factor(d$Block, sort(unique(d$Block)), sort(unique(d$Block)))
d$Gen <- as.factor(d$Gen)
d$Env <- as.factor(d$Env)

# dss <- droplevels(d[d$Man == "R" & d$Env == "Gai12",]) # Or Cra12, both fit fine in asreml and kmmr
# ggplot2::ggplot(dss, ggplot2::aes(x = as.numeric(as.character(Col)), y = as.numeric(as.character(Row)), shape = Rep, color = GY)) +
#   ggplot2::geom_point(size = 5) +
#   ggplot2::scale_color_continuous(palette = "viridis")
# 
# fit.asr <- asreml(GY ~ 1 + Gen,
#                   random = ~ ar1v(Row):ar1(Col),
#                   residual = ~ idv(units),
#                   data = dss,
#                   maxit = 25)
# 
# fit.kmmr <- kmmr(GY ~ 1 + Gen,
#                 random = ~ re(~ ar1(Row):ar1(Col), vars = ~ 1),
#                 residual = ~ units,
#                 data = dss,
#                 sparse = T,
#                 stepSize = c(0.1, 0.3, 0.5))
# 
# cor(fit.kmmr$u@x, fit.asr$coefficients$random[,])

# Now work on trial/environment specific random effects for single-stage analyses:
dss <- droplevels(d[d$Man == "R" & d$Env %in% c("Gai12", "Cra12"),])
ggplot2::ggplot(dss, ggplot2::aes(x = as.numeric(as.character(Col)), y = as.numeric(as.character(Row)), shape = Rep, color = GY)) +
  ggplot2::facet_wrap(ggplot2::vars(Env)) +
  ggplot2::geom_point(size = 5) +
  ggplot2::scale_color_continuous(palette = "viridis")

system.time(
fit.asr <- asreml(GY ~ -1 + Env,
                  random = ~ at(Env):ar1v(Row):ar1(Col) +
                             corgv(Env):vm(Gen, Ki.asr),
                  residual = ~ idv(units),
                  data = dss,
                  maxit = 30)
)

system.time(
fit.kmmr <- kmmr(fixed = GY ~ -1 + Env,
                 random = ~ atEach("Env", ~ re(~ ar1(Row):ar1(Col), vars = ~ 1)) +
                            re(~ corg(Env):lk(Gen, Ki.kmmr), vars = ~ 1),
                 residual = ~ units,
                 data = dss,
                 stepSize = c(0.1, 0.3, 0.5),
                 init = c(0.1193462, 0.1, 0.1, 0.1193462, 0.1, 0.1, 0.1078008, 0.5, 0.1078008),
                 isFixed = c(F, F, F, F, F, F, F, F, F),
                 sparse = TRUE,
                 maxit = 30)
)

# The last 492 BLUPs should be for the GxE part:
u.asr <- as.numeric(fit.asr$coefficients$random[(nrow(fit.asr$coefficients$random) - 492 + 1):nrow(fit.asr$coefficients$random),])
u.kmmr <- as.numeric(fit.kmmr$u@x[(length(fit.kmmr$u@x) - 492 + 1):length(fit.kmmr$u@x)])
plot(u.asr, u.kmmr)
cor(u.asr, u.kmmr)

fixed = GY ~ -1 + Env:Gen
random = ~ atEach("Env", ~ re(~ ar1(Row):ar1(Col), vars = ~ 1))
residual = ~ units
data = dss
maxit = 20
stepSize = NULL
avInfReg = NULL
weights = NULL
init = NULL
tolConvLL = 0.001
isFixed = NULL
verbose = TRUE
sparse = FALSE
timesToBoundary = 5
getStartingValues = FALSE
super = TRUE







