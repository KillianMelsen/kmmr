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

# Sample training and test sets:
sets <- vector("list", runs)
set.seed(1997)
for (run in 1:runs) {
  train.set <- sample(rownames(Ki), 150)
  test.set <- setdiff(rownames(Ki), train.set)
  sets[[run]] <- list(train.set = train.set, test.set = test.set)
}

# First do the standard FA-3 model with dense Ki:
isFactor <- kronecker(c(rep(F, 14), rep(T, 3)),
                      c(rep(F, 234)),
                      FUN = "|")
for (run in 1:runs) {
  d.train <- d
  d.train[d.train$Gen %in% sets[[run]]$test.set, "GY"] <- NA
  fit <- asreml(fixed = GY ~ -1 + Env,
                random = ~ fa(Env, 3):vm(Gen, Ki),
                residual = ~ units,
                data = d.train,
                weights = weight,
                family = asr_gaussian(dispersion = 1.0),
                trace = trace)
  preds <- as.data.frame(fit$coefficients$random[1:(length(levels(d.train$Env)) * nrow(Ki)),])
}






