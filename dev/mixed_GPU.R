library(Matrix)
{
D <- Matrix(readRDS("~/repos/kmmr/dev/exampledata_ED.rds"))
K <- Matrix(readRDS("~/repos/kmmr/dev/exampledata_K.rds"))
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
dss <- droplevels(d[d$Man == "R",])
colnames(dss)[colnames(dss) == "Rep"] <- "compBlock"
colnames(dss)[colnames(dss) == "Block"] <- "incompBlock"
dss$Rep <- ave(numeric(nrow(dss)), dss$Env, dss$Gen, FUN = seq_along)
dss$Rep <- factor(dss$Rep, 1:max(dss$Rep))
dss$compBlock <- as.numeric(dss$compBlock)
dss$compBlock[dss$compBlock == 4] <- 2
dss$compBlock[dss$compBlock == 5] <- 1
dss$compBlock <- as.factor(dss$compBlock)
D <- D[levels(dss$Env), levels(dss$Env)]
K <- K[levels(dss$Gen), levels(dss$Gen)]
Ki <- solve(K)
attr(Ki, "logDet") <- -1 * sum(log(eigen(Ki)$values))
}

fit.kmmr <- kmmr(fixed = GY ~ -1 + Env,
                 random = ~ atEach("Env", ~ re(~ id(compBlock), vars = ~ 1)) +
                            atEach("Env", ~ re(~ id(compBlock):id(incompBlock), vars = ~ 1)) +
                            atEach("Env", ~ re(~ ar1(Row):ar1(Col), vars = ~ 1)) +
                            re(~ gk(Env, D):lk(Gen, Ki), vars = ~ Env),
                 residual = ~ re(~ id(Env):id(Gen):id(Rep), vars = ~ Env),
                 data = dss,
                 stepSize = c(0.1, 0.2),#, 0.3, 0.4, 0.5),
                 avInfReg = rep(1.1, 2),
                 maxit = 2, timing = TRUE, useGPU = TRUE)

