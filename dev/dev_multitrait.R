d <- readRDS("dev/exampledata.rds")
d <- droplevels(d[d$Env == "Bol12",])
K <- readRDS("dev/exampledata_K.rds")
Ki.asr <- solve(K); attr(Ki.asr, "INVERSE") <- TRUE
Ki <- Matrix(solve(K))
# library(asreml)
# 
# microbenchmark::microbenchmark(fit.kmmr <- kmmr(fixed = GY ~ 1 + Man,
#                                                 random = ~ re(~ corg(Man):lk(Gen, Ki), vars = ~ Man),
#                                                 residual = ~ re(~ corg(Man):id(units), vars = ~ 1),
#                                                 data = d),
#                                fit.asr <- asreml(fixed = GY ~ 1 + Man,
#                                                  random = ~ corgh(Man):vm(Gen, Ki.asr),
#                                                  residual = ~ corgv(Man):id(units),
#                                                  asmv = Man,
#                                                  data = d),
#                                times = 20)
# 
# plot(fit.kmmr$u@x, fit.asr$coefficients$random[,])
# cor(fit.kmmr$u@x, fit.asr$coefficients$random[,])

set.seed(1997)
d2 <- d[sample(1:nrow(d), nrow(d)),]

# Need to figure out why the starting values are the same when the dataset is scrambled but LL and estimates
# are not identical...
# Scrambling doesn't matter when I use residual = ~ re(~ id(units), vars = ~ 1) and make sure
# data$units is just a factor like factor(sprintf("U%d", 1:nrow(data)), sprintf("U%d", 1:nrow(data)), sprintf("U%d", 1:nrow(data)))
# Not sure yet what goes wrong when scrambling and doing residual correlations.
fit.kmmr1 <- kmmr(fixed = GY ~ 1 + Man,
                  random = ~ re(~ corg(Man):lk(Gen, Ki), vars = ~ Man),
                  residual = ~ re(~ corg(Man):id(units), vars = ~ Man),
                  data = d)

fit.kmmr2 <- kmmr(fixed = GY ~ 1 + Man,
                  random = ~ re(~ corg(Man):lk(Gen, Ki), vars = ~ Man),
                  residual = ~ re(~ corg(Man):id(units), vars = ~ Man),
                  data = d2)

all.equal(fit.kmmr1$u, fit.kmmr2$u)
all.equal(fit.kmmr1$k, fit.kmmr2$k)
all.equal(fit.kmmr1$k[,1], fit.kmmr2$k[,1])

fit.asr <- asreml(fixed = GY ~ 1 + Man,
                  random = ~ corgh(Man):vm(Gen, Ki.asr),
                  residual = ~ units,
                  data = d2)

fixed = GY ~ -1 + Man
random = ~ re(~ corg(Man):lk(Gen, Ki), vars = ~ Man)
residual = ~ re(~ corg(Man):id(units), vars = ~ Man)
data = d
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