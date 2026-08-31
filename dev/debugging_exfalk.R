library(asreml)
asreml.options(workspace = "8000mb")

# DROPS ====
d <- readRDS("~/repos/kmmr/dev/exampledata.rds")
d <- droplevels(d[d$Man == "R",])
K <- readRDS("~/repos/kmmr/dev/exampledata_K.rds")
Ki <- Matrix(solve(K))
m <- 100
setup <- exfalkSetup(K, m = m)

fit1 <- kmmr(fixed = GY ~ -1 + Env,
             random = ~ re(~ exfa(Env, 3):exfalk(Gen, setup), vars = ~ Env),
             residual = ~ units,
             data = d,
             sparse = TRUE)

falki <- solve(tcrossprod(setup$L) + solve(setup$PSIi))
fit2 <- kmmr(fixed = GY ~ -1 + Env,
             random = ~ re(~ exfa(Env, 3):lk(Gen, falki), vars = ~ Env),
             residual = ~ units,
             data = d,
             sparse = TRUE)

system.time(
fit3 <- kmmr(fixed = GY ~ -1 + Env,
             random = ~ re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
             residual = ~ units,
             data = d,
             sparse = TRUE)
)

isZero <- kronecker(c(rep(F, length(levels(d$Env))), rep(T, 3)),
                    c(rep(F, nrow(K)), rep(T, m)),
                    FUN = "|")

plot(fit1$u[which(!isZero),], fit2$u[1:3444,])
plot(fit1$u[which(!isZero),], fit3$u[1:3444,])

cor(fit1$u[which(!isZero),], fit2$u[1:3444,])
cor(fit1$u[which(!isZero),], fit3$u[1:3444,])

# CIMMYT ====
d <- read.csv("dev/Krause_et_al_2018_Yield_BLUEs.csv")
d <- droplevels(d[d$Breeding.Cycle == "2016-17",])
d$GID <- paste0("G", d$GID)
d$Man <- factor(d$Managed_Treatment, unique(d$Managed_Treatment), unique(d$Managed_Treatment))
d$Gen <- factor(d$GID, unique(d$GID), unique(d$GID))
d <- d[, c("Man", "Gen", "Grain_Yield_BLUE")]
colnames(d)[3] <- "GY"
K <- Matrix(readRDS("dev/Krause_K.rds"))
rownames(K) <- colnames(K) <- paste0("G", rownames(K))
K <- K[levels(d$Gen), levels(d$Gen)]
Ki.asr <- solve(as.matrix(K)); attr(Ki.asr, "INVERSE") <- TRUE
Ki.kmr <- solve(K)

# system.time(
#   fit.asr <- asreml(fixed = GY ~ -1 + Man,
#                     random = ~ fa(Man, 2):vm(Gen, Ki.asr),
#                     residual = ~ units,
#                     data = d,
#                     maxit = 25)
# ) # 181.533 seconds, 14 iterations, 12.97 seconds/iter
#
# system.time(
#   fit.kmr <- kmmr(fixed = GY ~ -1 + Man,
#                   random = ~ re(~ exfa(Man, 2):lk(Gen, Ki.kmr), vars = ~ Man),
#                   residual = ~ units,
#                   data = d,
#                   maxit = 25,
#                   stepSize = c(0.15, 0.25, 0.5, 0.75, rep(1, 21)),
#                   avInfReg = c(1.1, 1.1, 1.05, 1.05, rep(1, 21)),
#                   sparse = TRUE)
# )  # 96.837 seconds, 15 iterations, 6.46 seconds/iter
# plot(fit.kmr$u[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d)])
# cor(fit.kmr$u[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d)])
# 
# system.time(
#   fit.kmr2 <- kmmr(fixed = GY ~ -1 + Man,
#                   random = ~ re(~ exfa(Man, 2):lk(Gen, Ki.kmr), vars = ~ Man),
#                   residual = ~ units,
#                   data = d,
#                   init = summary(fit.kmr)$varcomp$Estimate,
#                   isFixed = summary(fit.kmr)$varcomp$Estimate == 0,
#                   maxit = 25,
#                   stepSize = c(0.15, 0.25, 0.5, 0.75, rep(1, 21)),
#                   avInfReg = c(1.1, 1.1, 1.05, 1.05, rep(1, 21)),
#                   sparse = TRUE)
# ) # 69.555 seconds, 11 iterations, 6.32 seconds/iter
# plot(fit.kmr2$u[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d)])
# cor(fit.kmr2$u[1:nrow(d)], fit.asr$coefficients$random[1:nrow(d)])
# 
# # You cannot directly compare the variance components of asreml and kmmr as they
# # are estimated on different scales!
# # You have to compare the FA implied covariance matrices:
# covmat.asr <- matrix(summary(fit.asr)$varcomp[6:15, 1], 5, 2) %*% t(matrix(summary(fit.asr)$varcomp[6:15, 1], 5, 2)) + diag(summary(fit.asr)$varcomp[1:5, 1])
# cormat.asr <- cov2cor(covmat.asr)
# cormat.kmr2 <- matrix(summary(fit.kmr2)$varcomp[6:15, 2], 5, 2) %*% t(matrix(summary(fit.kmr2)$varcomp[6:15, 2], 5, 2))
# diag(cormat.kmr2) <- 1
# covmat.kmr2 <- outer(sqrt(summary(fit.kmr2)$varcomp[1:5, 2]), sqrt(summary(fit.kmr2)$varcomp[1:5, 2])) * cormat.kmr2
# 
# m <- 100
# setup <- exfalkSetup(K, m = m)
# system.time(
#   fit.kmr.exfalk <- kmmr(fixed = GY ~ -1 + Man,
#                          random = ~ re(~ exfa(Man, 2):exfalk(Gen, setup), vars = ~ Man),
#                          residual = ~ units,
#                          data = d,
#                          maxit = 25,
#                          stepSize = c(0.15, 0.25, 0.5, 0.75, rep(1, 21)),
#                          avInfReg = c(1.1, 1.1, 1.05, 1.05, rep(1, 21)),
#                          sparse = TRUE)
# )
# 
# isZero <- kronecker(c(rep(F, length(levels(d$Man))), rep(T, 2)),
#                     c(rep(F, nrow(K)), rep(T, m)),
#                     FUN = "|")
# 
# plot(fit.kmr.exfalk$u[which(!isZero),], fit.kmr2$u[1:nrow(d)])
# cor(fit.kmr.exfalk$u[which(!isZero),], fit.kmr2$u[1:nrow(d)])



runs <- 50
results <- expand.grid(Run = 1:runs,
                       Model = c("Full", "exfalk", "asreml"),
                       Man = levels(d$Man),
                       cor = NA,
                       LL = NA,
                       time = NA)
setup <- exfalkSetup(K, m = 100, method = "fad")

## Sampling training/test sets ====
test.set <- vector("list", runs)
set.seed(1997)
for (run in 1:runs) {
  test.set[[run]] <- sample(rownames(K), 200)
}

for (run in 1:runs) {
  cat(sprintf("Run = %d\n\n", run))
  d.train <- d
  d.train[d.train$Gen %in% test.set[[run]], "GY"] <- NA
  
  ## Fitting models ====
  # start <- proc.time()[3]
  # fit.asr <- asreml(fixed = GY ~ -1 + Man,
  #                   random = ~ fa(Man, 2):vm(Gen, Ki.asr),
  #                   residual = ~ units,
  #                   data = d.train,
  #                   maxit = 25)
  # stop <- proc.time()[3]
  # time.asr <- stop - start
  # LL.asr <- fit.asr$loglik
  
  start <- proc.time()[3]
  fit.kmr <- kmmr(fixed = GY ~ -1 + Man,
                  random = ~ re(~ exfa(Man, 2):lk(Gen, Ki.kmr), vars = ~ Man),
                  residual = ~ units,
                  data = d.train,
                  maxit = 25,
                  stepSize = c(0.15, 0.25, 0.5, 0.75, rep(1, 21)),
                  avInfReg = c(1.1, 1.1, 1.05, 1.05, rep(1, 21)),
                  sparse = TRUE)
  stop <- proc.time()[3]
  time.kmr <- stop - start
  LL.kmr <- fit.kmr$LL[length(fit.kmr$LL)]

  start <- proc.time()[3]
  fit.kmr.exfalk <- kmmr(fixed = GY ~ -1 + Man,
                         random = ~ re(~ exfa(Man, 2):exfalk(Gen, setup), vars = ~ Man),
                         residual = ~ units,
                         data = d.train,
                         maxit = 25,
                         stepSize = c(0.15, 0.25, 0.5, 0.75, rep(1, 21)),
                         avInfReg = c(1.1, 1.1, 1.05, 1.05, rep(1, 21)),
                         sparse = TRUE)
  stop <- proc.time()[3]
  time.kmr.exfalk <- stop - start
  LL.kmr.exfalk <- fit.kmr.exfalk$LL[length(fit.kmr.exfalk$LL)]

  ## Gathering predictions ====
  isZero <- kronecker(c(rep(F, length(levels(d$Man))), rep(T, 2)), c(rep(F, nrow(K)), rep(T, setup$m)), FUN = "|")

  # preds.asr <- fit.asr$coefficients$random[1:nrow(d),]
  # preds.asr <- data.frame(Man = gsub("fa\\(Man, 2\\)_(.*):vm\\(Gen, Ki.asr\\)_(.*)", "\\1", names(preds.asr)),
  #                         Gen = gsub("fa\\(Man, 2\\)_(.*):vm\\(Gen, Ki.asr\\)_(.*)", "\\2", names(preds.asr)),
  #                         pred = as.numeric(preds.asr))
  
  preds.kmr <- fit.kmr$u[1:nrow(d),, drop = FALSE]
  preds.kmr <- data.frame(Man = gsub("(.*):(.*)", "\\1", rownames(preds.kmr)),
                          Gen = gsub("(.*):(.*)", "\\2", rownames(preds.kmr)),
                          pred = preds.kmr@x)

  preds.kmr.exfalk <- fit.kmr.exfalk$u[which(!isZero),, drop = FALSE]
  preds.kmr.exfalk <- data.frame(Man = gsub("(.*):(.*)", "\\1", rownames(preds.kmr.exfalk)),
                                 Gen = gsub("(.*):(.*)", "\\2", rownames(preds.kmr.exfalk)),
                                 pred = preds.kmr.exfalk@x)

  # preds.asr <- droplevels(preds.asr[preds.asr$Gen %in% test.set[[run]],])
  preds.kmr <- droplevels(preds.kmr[preds.kmr$Gen %in% test.set[[run]],])
  preds.kmr.exfalk <- droplevels(preds.kmr.exfalk[preds.kmr.exfalk$Gen %in% test.set[[run]],])
  
  ## Computing accuracies ====
  target <- droplevels(d[d$Gen %in% test.set[[run]],])
  for (i in 1:length(levels(d$Man))) {
    target.ss <- droplevels(target[target$Man == levels(target$Man)[i],])
    # preds.asr.ss <- droplevels(preds.asr[preds.asr$Man == levels(target$Man)[i],])
    preds.kmr.ss <- droplevels(preds.kmr[preds.kmr$Man == levels(target$Man)[i],])
    preds.kmr.exfalk.ss <- droplevels(preds.kmr.exfalk[preds.kmr.exfalk$Man == levels(target$Man)[i],])
    
    # preds.asr.ss <- preds.asr.ss[match(as.character(target.ss$Gen), preds.asr.ss$Gen),]
    preds.kmr.ss <- preds.kmr.ss[match(as.character(target.ss$Gen), preds.kmr.ss$Gen),]
    preds.kmr.exfalk.ss <- preds.kmr.exfalk.ss[match(as.character(target.ss$Gen), preds.kmr.exfalk.ss$Gen),]

    # results[results$Run == run & results$Model == "asreml" & results$Man == levels(target$Man)[i], "cor"] <- cor(preds.asr.ss$pred, target.ss$GY)
    results[results$Run == run & results$Model == "Full" & results$Man == levels(target$Man)[i], "cor"] <- cor(preds.kmr.ss$pred, target.ss$GY)
    results[results$Run == run & results$Model == "exfalk" & results$Man == levels(target$Man)[i], "cor"] <- cor(preds.kmr.exfalk.ss$pred, target.ss$GY)
    
    # results[results$Run == run & results$Model == "asreml" & results$Man == levels(target$Man)[i], "LL"] <- LL.asr
    results[results$Run == run & results$Model == "Full" & results$Man == levels(target$Man)[i], "LL"] <- LL.kmr
    results[results$Run == run & results$Model == "exfalk" & results$Man == levels(target$Man)[i], "LL"] <- LL.kmr.exfalk
    
    # results[results$Run == run & results$Model == "asreml" & results$Man == levels(target$Man)[i], "time"] <- time.asr
    results[results$Run == run & results$Model == "Full" & results$Man == levels(target$Man)[i], "time"] <- time.kmr
    results[results$Run == run & results$Model == "exfalk" & results$Man == levels(target$Man)[i], "time"] <- time.kmr.exfalk
  }
}
saveRDS(results, "dev/results.exfalk.CIMMYT.fad.rds")
aggregate(results, cor ~ Model, FUN = mean)
aggregate(results, cor ~ Model, FUN = median)

# AGT ====
d <- read.csv("dev/Grain.yield_adjusted.csv")
d <- d[, c(1, 4)]
# M <- read.csv("dev/genetic_marker_matrix_M.csv")
# rownames(M) <- d$Genotype
# gdata <- statgenGWAS::createGData(geno = M)
# gdata <- statgenGWAS::codeMarkers(gdata, removeDuplicates = TRUE)
# K <- statgenGWAS::kinship(gdata$markers, method = "IBS")
# K <- Matrix(K)
# saveRDS(K, "dev/K.AGT.rds")
# K <- readRDS("dev/K.AGT.rds")
colnames(d) <- c("Gen", "GY")
d$Gen <- factor(d$Gen, d$Gen, d$Gen)
# K <- K[levels(d$Gen), levels(d$Gen)]
# Ki <- solve(K)
# rm(K); gc()
Ki <- readRDS("dev/Ki.AGT.rds")

# 30 seconds per iteration
fit2 <- kmmr(fixed = GY ~ 1,
            random = ~ re(~ lk(Gen, Ki), vars = ~ 1),
            residual = ~ units,
            data = d,
            sparse = FALSE,
            useGPU = TRUE)

# scaling <- diag(K)
# Kcor <- cov2cor(as.matrix(K))
# Kcor.ss <- Kcor[1:200, 1:200]

# 10 factors:
#       p    t10    t20    t100
# 1   200   0.10
# 2   400   0.25
# 3   800   1.40
# 4  1600   6.20
# 5  3200  30.00
# 6  6400  78.00 
# 7 10375 281.00 606.00 2839.00
# system.time(fit <- fad::fad(factors = 40, covmat = Kcor.ss))

system.time(setup <- exfalkSetup(K, 200, "fad"))
fit2 <- kmmr(fixed = GY ~ 1,
             random = ~ re(~ exfalk(Gen, setup), vars = ~ 1),
             residual = ~ units,
             data = d,
             sparse = TRUE)






