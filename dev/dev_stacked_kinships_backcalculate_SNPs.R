library(rrBLUP)
library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = TRUE, score = TRUE)
library(ggplot2)
K <- readRDS("dev/DROPS_K.rds")
M <- readRDS("dev/DROPS_M.rds")
d <- readRDS("dev/exampledata.rds")[, c("Gen", "Man", "Env", "GY")]
d <- droplevels(d[d$Man == "R" & d$Env == "Bol12",])

K <- K[levels(d$Gen), levels(d$Gen)]
M <- M[levels(d$Gen),]

fit.snpblup <- rrBLUP::mixed.solve(y = d$GY, Z = M, X = matrix(1, 246, 1))
fit.gblup <- rrBLUP::kin.blup(data = d, geno = "Gen", pheno = "GY", K = K)
u.snp <- fit.snpblup$u
u.g <- fit.gblup$g

u.g.fromsnp <- M %*% matrix(u.snp)
u.snp.fromg <- t(M) %*% solve(K) %*% matrix(u.g)

plot(u.g.fromsnp, u.g); abline(0, 1, col = "red")
plot(u.snp.fromg, u.snp); abline(0, 1, col = "red")

# Now try with the group-regularized GBLUP / multi-GBLUP:
M.info <- read.csv("~/repos/kernels_2026/DROPS/raw_data/7b-InfoSNP_50K_41722.csv")
M.info$SNP.names <- gsub("-", ".", M.info$SNP.names)
M.info <- droplevels(M.info[M.info$SNP.names %in% colnames(M),])
SNPs <- vector("list", 10)
for (chr in 1:10) {
  SNPs[[chr]] <- M[, c(M.info[M.info$Chromosome == chr, "SNP.names"])]
}
sum(unlist(lapply(SNPs, ncol))) == ncol(M)
gdatas <- Ks <- vector("list", 10)
for (chr in 1:10) {
  Ks[[chr]] <- statgenGWAS::kinship(SNPs[[chr]], "IBS")[levels(d$Gen), levels(d$Gen)]
}
# Ks for chr 6 and 10 are singular:
matrix(unlist(lapply(Ks, function(K) min(eigen(K)$values))))
Ks[[6]] <- Ks[[6]] + diag(1e-4, 246, 246)
Ks[[10]] <- Ks[[10]] + diag(1e-4, 246, 246)

fit.asr <- asreml(fixed = GY ~ 1,
                  random = ~ vm(Gen, Ks[[1]]) + vm(Gen, Ks[[2]]) + vm(Gen, Ks[[3]]) + vm(Gen, Ks[[4]]) +
                    vm(Gen, Ks[[5]]) + vm(Gen, Ks[[6]]) + vm(Gen, Ks[[7]]) + vm(Gen, Ks[[8]]) + vm(Gen, Ks[[9]]) + vm(Gen, Ks[[10]]),
                  residual = ~ idv(units),
                  data = d)

u.g.grgblup <- vector("list", 10)
i <- 1
for (chr in 1:10) {
  u.g.grgblup[[chr]] <- fit.asr$coefficients$random[i:(i - 1 + 246), 1, drop = FALSE]
  rownames(u.g.grgblup[[chr]]) <- levels(d$Gen)
  i <- i + 246
}
plot(unlist(lapply(u.g.grgblup, as.numeric))) # These are chromosome specific / group-regularized GEBVs

# Now convert these to chromosome specific / group-regularized SNP effects:
Ksi <- lapply(Ks, solve)
u.snp.fromgrgblup <- lapply(1:10, function(chr) t(SNPs[[chr]]) %*% Ksi[[chr]] %*% u.g.grgblup[[chr]])
plot(unlist(lapply(u.snp.fromgrgblup, as.numeric))) # These are chromosome specific / group-regularized SNP effects
grsnpeffects <- unlist(lapply(u.snp.fromgrgblup, function(u) u[,]))
effects <- data.frame(effect = as.numeric(grsnpeffects), snp = M.info$SNP.names, chr = as.factor(M.info$Chromosome), pos = M.info$Position)
ggplot(effects, aes(x = 1:nrow(effects), y = abs(effect), color = chr)) + geom_point()

ggplot(effects, aes(x = pos, y = abs(effect), color = abs(effect))) +
  facet_wrap(vars(chr)) + geom_point() +
  scale_color_continuous(palette = "viridis")

pca <- prcomp(K)
summary(pca)
covariates <- pca$x[, 1:5]
covariates <- covariates[match(levels(d$Gen), rownames(covariates)),]
d2 <- cbind(d, covariates)

fit.asr2 <- asreml(fixed = GY ~ 1 + PC1 + PC2 + PC3 + PC4 + PC5,
                   random = ~ vm(Gen, Ks[[1]]) + vm(Gen, Ks[[2]]) + vm(Gen, Ks[[3]]) + vm(Gen, Ks[[4]]) +
                     vm(Gen, Ks[[5]]) + vm(Gen, Ks[[6]]) + vm(Gen, Ks[[7]]) + vm(Gen, Ks[[8]]) + vm(Gen, Ks[[9]]) + vm(Gen, Ks[[10]]),
                   residual = ~ idv(units),
                   data = d2)

u.g.grgblup2 <- vector("list", 10)
i <- 1
for (chr in 1:10) {
  u.g.grgblup2[[chr]] <- fit.asr2$coefficients$random[i:(i - 1 + 246), 1, drop = FALSE]
  rownames(u.g.grgblup2[[chr]]) <- levels(d2$Gen)
  i <- i + 246
}
plot(unlist(lapply(u.g.grgblup2, as.numeric))) # These are chromosome specific / group-regularized GEBVs

# Now convert these to chromosome specific / group-regularized SNP effects:
u.snp.fromgrgblup2 <- lapply(1:10, function(chr) t(SNPs[[chr]]) %*% Ksi[[chr]] %*% u.g.grgblup2[[chr]])
plot(unlist(lapply(u.snp.fromgrgblup2, as.numeric))) # These are chromosome specific / group-regularized SNP effects
grsnpeffects2 <- unlist(lapply(u.snp.fromgrgblup2, function(u) u[,]))
effects2 <- data.frame(effect = as.numeric(grsnpeffects2), snp = M.info$SNP.names, chr = as.factor(M.info$Chromosome), pos = M.info$Position)
ggplot(effects2, aes(x = 1:nrow(effects2), y = abs(effect), color = chr)) + geom_point()

ggplot(effects2, aes(x = pos, y = abs(effect), color = abs(effect))) +
  facet_wrap(vars(chr)) + geom_point() +
  scale_color_continuous(palette = "viridis")

ggsave(filename = "dev/effects2.png")
