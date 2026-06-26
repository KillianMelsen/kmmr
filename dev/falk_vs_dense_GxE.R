library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", maxit = 25)
rotateInitL <- function(L) {
  L.rotated <- L
  m <- ncol(L.rotated)
  p <- nrow(L.rotated)
  coords <- which(upper.tri(matrix(0, p, m)), arr.ind = TRUE)
  givens <- vector("list", nrow(coords))
  for (i in 1:nrow(coords)) {
    row <- coords[i, "row"]
    col <- coords[i, "col"]
    givens[[i]] <- diag(m)
    h <- sqrt(sum(c(L.rotated[row, row], L.rotated[row, col])^2))
    c <- L.rotated[row, row] / h
    s <- L.rotated[row, col] / h
    diag(givens[[i]])[c(row, col)] <- c
    givens[[i]][row, col] <- -s
    givens[[i]][col, row] <- s
    L.rotated <- L.rotated %*% givens[[i]]
  }
  L.rotated[upper.tri(L.rotated)] <- 0
  rotmat <- Reduce("%*%", givens)
  return(list(L.rotated = L.rotated, rotmat = rotmat))
}
createCOO <- function(m, K) {
  sds <- diag(sqrt(diag(K)))
  FA <- factanal(covmat = K, factors = m)
  L <- rotateInitL(FA$loadings[,, drop = F])$L.rotated
  PSI <- diag(FA$uniquenesses)
  L <- sds %*% L
  PSI <- sds %*% PSI %*% sds
  FAKi <- rbind(cbind(solve(PSI), solve(PSI) %*% L),
                cbind(t(L) %*% solve(PSI), diag(m) + t(L) %*% solve(PSI) %*% L))
  rownames(FAKi) <- colnames(FAKi) <- c(rownames(K), paste0("GF_", 1:m))
  COO <- as.data.frame(which(lower.tri(FAKi, diag = TRUE), arr.ind = TRUE))
  COO$value <- FAKi[cbind(COO$row, COO$col)]
  colnames(COO) <- c("rowid", "colid", "vals")
  COO <- COO[!(COO$vals == 0),]
  COO <- COO[order(COO$rowid, COO$colid),]
  attr(COO, "INVERSE") <- TRUE
  attr(COO, "rowNames") <- rownames(FAKi)
  return(COO)
}
d <- readRDS("~/repos/kernels_2026/DROPS/data/ydata.SE.rds")
d <- droplevels(d[d$Man == "R",])
K <- readRDS("~/repos/kernels_2026/DROPS/data/K.SE.rds")
Ki <- solve(K)
attr(Ki, "INVERSE") <- TRUE
runs <- 20
trace <- TRUE
results1 <- expand.grid(Run = 1:runs,
                        Model = "Dense",
                        m = NA,
                        Env = levels(d$Env),
                        cor_pearson = NA,
                        RMSE = NA)

# Sample training and test sets:
sets <- vector("list", runs)
set.seed(1997)
for (run in 1:runs) {
  train.set <- sample(rownames(Ki), 150)
  test.set <- setdiff(rownames(Ki), train.set)
  sets[[run]] <- list(train.set = train.set, test.set = test.set)
}

# First do the standard FA-3 model with dense Ki ====
isFactor <- kronecker(c(rep(FALSE,  14), rep(TRUE, 3)),
                      c(rep(FALSE, 202), rep(TRUE, 0)),
                      FUN = "|")
run <- 1
for (run in 1:runs) {
  d.train <- d
  d.train[d.train$Gen %in% sets[[run]]$test.set, "GY"] <- NA
  target <- droplevels(d[d$Gen %in% sets[[run]]$test.set,])
  fit <- asreml(fixed = GY ~ -1 + Env,
                random = ~ fa(Env, 3):vm(Gen, Ki),
                residual = ~ units,
                data = d.train,
                weights = weight,
                family = asr_gaussian(dispersion = 1.0),
                trace = trace)
  
  preds <- as.data.frame(fit$coefficients$random[1:nrow(d.train),])
  names(preds) <- "predicted.value"
  preds$Env <- as.factor(gsub("fa\\(Env, 3\\)_(.*):vm\\(Gen, Ki\\)_(.*)", "\\1", rownames(preds)))
  preds$Gen <- as.factor(gsub("fa\\(Env, 3\\)_(.*):vm\\(Gen, Ki\\)_(.*)", "\\2", rownames(preds)))
  preds$EnvGen <- as.factor(paste(preds$Env, preds$Gen, sep = ":"))
  preds <- preds[match(target$EnvGen, preds$EnvGen),]
  
  for (e in levels(d.train$Env)) {
    preds2 <- droplevels(preds[preds$Env == e, c("Gen", "predicted.value")])
    target2 <- droplevels(target[target$Env == e, c("Gen", "GY")])
    target2 <- target2[match(preds2$Gen, target2$Gen),]
    
    # Accuracy measures:
    target2$GY <- scale(target2$GY, scale = FALSE)
    cor_pearson <- cor(preds2$predicted.value, target2$GY)
    RMSE <- sqrt(mean((preds2$predicted.value - target2$GY)^2))
    
    result.row <- which(results1$Run == run & results1$Model == "Dense" & results1$Env == e)
    results1[result.row, c("cor_pearson", "RMSE")] <- c(cor_pearson, RMSE)
  }
}

# Now the models with approximated kinships ====
factors <- c(30, 40)
results3 <- expand.grid(Run = 1:runs,
                        Model = "FALK",
                        m = factors,
                        Env = levels(d$Env),
                        cor_pearson = NA,
                        RMSE = NA)

for (m in factors) {
  COO <- createCOO(m, K)
  isFactor <- kronecker(c(rep(FALSE,  14), rep(TRUE, 3)),
                        c(rep(FALSE, 234), rep(TRUE, m)),
                        FUN = "|")
  for (run in 1:runs) {
    d.train <- d
    d.train[d.train$Gen %in% sets[[run]]$test.set, "GY"] <- NA
    target <- droplevels(d[d$Gen %in% sets[[run]]$test.set,])
    fit <- asreml(fixed = GY ~ -1 + Env,
                  random = ~ fa(Env, 3):vm(Gen, COO),
                  residual = ~ units,
                  data = d.train,
                  weights = weight,
                  family = asr_gaussian(dispersion = 1.0),
                  trace = trace)
    
    preds <- as.data.frame(fit$coefficients$random[!isFactor,])
    names(preds) <- "predicted.value"
    preds$Env <- as.factor(gsub("fa\\(Env, 3\\)_(.*):vm\\(Gen, COO\\)_(.*)", "\\1", rownames(preds)))
    preds$Gen <- as.factor(gsub("fa\\(Env, 3\\)_(.*):vm\\(Gen, COO\\)_(.*)", "\\2", rownames(preds)))
    preds$EnvGen <- as.factor(paste(preds$Env, preds$Gen, sep = ":"))
    preds <- preds[match(target$EnvGen, preds$EnvGen),]
    
    for (e in levels(d.train$Env)) {
      preds2 <- droplevels(preds[preds$Env == e, c("Gen", "predicted.value")])
      target2 <- droplevels(target[target$Env == e, c("Gen", "GY")])
      target2 <- target2[match(preds2$Gen, target2$Gen),]
      
      # Accuracy measures:
      target2$GY <- scale(target2$GY, scale = FALSE)
      cor_pearson <- cor(preds2$predicted.value, target2$GY)
      RMSE <- sqrt(mean((preds2$predicted.value - target2$GY)^2))
      
      result.row <- which(results3$Run == run & results3$Model == "FALK" & results3$Env == e & results3$m == m)
      results3[result.row, c("cor_pearson", "RMSE")] <- c(cor_pearson, RMSE)
    }
  }
}

results1$m <- 0
results <- rbind(results1,
                 results2,
                 results3)

x <- aggregate(results, cbind(cor_pearson, RMSE) ~ Model + m + Env, FUN = mean)
x <- droplevels(x[x$m %in% c(0, 50),])
plot(x[x$Model == "Dense", "cor_pearson"],
     x[x$Model == "FALK", "cor_pearson"])
abline(0, 1, col = "red")

















