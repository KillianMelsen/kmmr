library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", maxit = 25)
set.seed(1997)

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

d <- readRDS("dev/exampledata.rds")
# d <- readRDS("~/repos/kernels_2026/DROPS/data/ydata.SE.rds")
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY")])
K <- readRDS("dev/exampledata_K.rds")
# diag(K) <- diag(K) + 1e-3
# K <- readRDS("~/repos/kernels_2026/DROPS/data/K.SE.rds")
Ki <- solve(K)
attr(Ki, "INVERSE") <- TRUE
factors <- c(5, 10, 20, 50)
runs <- 20
results <- expand.grid(m = factors,
                       Run = 1:runs,
                       Model = c("Dense", "FALK"),
                       Env = levels(d$Env),
                       cor_pearson = NA,
                       RMSE = NA)
trace <- FALSE

start <- Sys.time()
for (m in factors) {
  COO <- createCOO(m, K)
  ei <- 1
  for (e in levels(d$Env)) {
    d.e <- droplevels(d[d$Env == e,])
    for (run in 1:runs) {
      cat(sprintf("Running m = %d, Env = %d, run = %d\n", m, ei, run))
      train.set <- sample(levels(d.e$Gen), 150)
      test.set <- setdiff(levels(d.e$Gen), train.set)
      d.e.train <- d.e
      d.e.train[d.e.train$Gen %in% test.set, "GY"] <- NA
      target <- droplevels(d.e[d.e$Gen %in% test.set,])
      target$GY <- scale(target$GY, center = TRUE, scale = FALSE)
      
      if (m == factors[1]) {
        try({
          fit.dense <- asreml(fixed = GY ~ 1,
                              random = ~ vm(Gen, Ki),
                              residual = ~ units,
                              data = d.e.train,
                              trace = trace)#,
                              # weights = weight,
                              # family = asr_gaussian(dispersion = 1.0))
          
          preds <- fit.dense$coefficients$random[1:246,]
          names(preds) <- gsub("vm\\(Gen, Ki\\)_(.*)", "\\1", names(preds))
          preds <- preds[match(as.character(target$Gen), names(preds))]
          results[results$m == m & results$Run == run & results$Env == e & results$Model == "Dense", "cor_pearson"] <- cor(preds, target$GY)
          results[results$m == m & results$Run == run & results$Env == e & results$Model == "Dense", "RMSE"] <- sqrt(mean((preds - target$GY)^2))
        })
      }
      
      try({
        fit.FALK <- asreml(fixed = GY ~ 1,
                           random = ~ vm(Gen, COO),
                           residual = ~ units,
                           data = d.e.train,
                           trace = trace)#,
                           # weights = weight,
                           # family = asr_gaussian(dispersion = 1.0))
        
        preds <- fit.FALK$coefficients$random[1:246,]
        names(preds) <- gsub("vm\\(Gen, COO\\)_(.*)", "\\1", names(preds))
        preds <- preds[match(as.character(target$Gen), names(preds))]
        results[results$m == m & results$Run == run & results$Env == e & results$Model == "FALK", "cor_pearson"] <- cor(preds, target$GY)
        results[results$m == m & results$Run == run & results$Env == e & results$Model == "FALK", "RMSE"] <- sqrt(mean((preds - target$GY)^2))
      })
      
    }
    ei <- ei + 1
  }
}
end <- Sys.time()
end - start

# aggregate(results, cbind(cor_pearson, RMSE) ~ Model + m, FUN = mean)





