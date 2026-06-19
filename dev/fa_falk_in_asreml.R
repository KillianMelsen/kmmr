library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = TRUE, score = TRUE)
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

m <- 10
d <- readRDS("dev/exampledata.rds")
d <- droplevels(d[d$Man == "R", c("Env", "Gen", "GY")])
K <- readRDS("dev/exampledata_K.rds")
sds <- diag(sqrt(diag(K)))
FA <- factanal(covmat = K, factors = m)
L <- rotateInitL(FA$loadings[,, drop = F])$L.rotated
# L <- FA$loadings[,, drop = F]
PSI <- diag(FA$uniquenesses)
L <- sds %*% L
PSI <- sds %*% PSI %*% sds
FAKi <- rbind(cbind(solve(PSI), solve(PSI) %*% L),
              cbind(t(L) %*% solve(PSI), diag(m) + t(L) %*% solve(PSI) %*% L))

rownames(FAKi) <- colnames(FAKi) <- c(rownames(K), paste0("GF_", 1:m))
FAK <- L %*% t(L) + PSI
rownames(FAK) <- colnames(FAK) <- rownames(K)
Ki <- solve(K)
attr(Ki, "INVERSE") <- TRUE

# Store lower triangle and diagonal of FAKi in COO triplet format:
COO <- as.data.frame(which(lower.tri(FAKi, diag = TRUE), arr.ind = TRUE))
COO$value <- FAKi[cbind(COO$row, COO$col)]
colnames(COO) <- c("rowid", "colid", "vals")
COO <- COO[!(COO$vals == 0),]
COO <- COO[order(COO$rowid, COO$colid),]
attr(COO, "INVERSE") <- TRUE
attr(COO, "rowNames") <- rownames(FAKi)

isZero <- kronecker(c(rep(1, 14), rep(2, 2)),
                    c(rep(1, 246), rep(2, m))) > 1

system.time(fit1 <- asreml(fixed = GY ~ -1 + Env,
                           random = ~ fa(Env, 2):vm(Gen, Ki),
                           residual = ~ idv(units),
                           data = d, maxit = 25))

system.time(fit2 <- asreml(fixed = GY ~ -1 + Env,
                           random = ~ fa(Env, 2):vm(Gen, COO),
                           residual = ~ idv(units),
                           data = d, maxit = 25))

u1 <- matrix(fit1$coefficients$random[1:3444])
u2 <- matrix(fit2$coefficients$random[!isZero])
cor(u1, u2, method = "pearson")
cor(u1, u2, method = "spearman")
plot(u1, u2)
