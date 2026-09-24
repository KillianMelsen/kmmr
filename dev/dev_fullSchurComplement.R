d <- kmmr::dataDROPS
K <- kmmr::kinshipDROPS
ED <- kmmr::distmatDROPS
Ki <- solve(K)
attr(Ki, "logDet") <- sum(log(eigen(K, only.values = TRUE)$values))
d <- droplevels(d[d$Man == "R" & d$Env %in% levels(d$Env)[1:10],])
ED <- ED[levels(d$Env), levels(d$Env)]

# fit <- kmmr(fixed = GY ~ 1 + Env,
#             random = ~ re(~ exfa(Env, 1):lk(Gen, Ki), vars = ~ Env) + re(~ gk(Env, ED):lk(Gen, Ki), vars = ~ Env),
#             residual = ~ units,
#             data = d, sparse = TRUE, maxit = 20, avInfReg = rep(1.1, 20))

fixed = GY ~ 1 + Env
random = ~ re(~ exfa(Env, 1):lk(Gen, Ki), vars = ~ Env) + re(~ gk(Env, ED):lk(Gen, Ki), vars = ~ Env)
residual = ~ units
data = d
maxit = 20
avInfReg = rep(1.1, 20)
stepSize = NULL
weights = NULL
init = NULL
tolConvLL = 0.001
isFixed = NULL
verbose = TRUE
timesToBoundary = 5
getStartingValues = FALSE
super = TRUE
useGPU = FALSE
denseThreshold = 0.2
debugging = FALSE
