library(asreml)
asreml.options(workspace = "5000mb", pworkspace = "5000mb", wvr = TRUE, score = TRUE)

# Seed:
set.seed(1997)

# Loading example data (DROPS):
d <- readRDS("dev/exampledata.rds")[, c("Gen", "Man", "Env", "GY")]
d <- droplevels(d[d$Man == "R" & d$Env == "Bol12",])
# GD <- readRDS("dev/exampledata_GD.rds")[levels(d$Gen), levels(d$Gen)]
# K1 <- readRDS("dev/exampledata_K.rds")[levels(d$Gen), levels(d$Gen)]
# K2 <- exp(-0.1 * GD)
# K1i <- solve(K1)
# K2i <- solve(K2)
# attr(K1i, "INVERSE") <- TRUE
# attr(K2i, "INVERSE") <- TRUE

# fit.asr <- asreml(fixed = GY ~ 1,
#                   random = ~ vm(Gen, K1) + vm(Gen, K2),
#                   residual = ~ idv(units),
#                   data = d)

# This doesn't work for whatever reason. It should work...
# vf <- function(order, kappa) {
#   V <- kappa[1] * K1 + kappa[2] * K2
#   return(list(V, K1, K2))
# }
# init <- c(0.2, 0.2)
# type <- c("V", "V")
# cons <- c("P", "P")
# fit2 <- asreml(fixed = GY ~ 1,
#                random = ~ own(Gen, "vf", init, type, cons),
#                residual = ~ idv(units),
#                data = d)

# summary(fit.asr)$varcomp

# Works!
# fit <- kmmr(fixed = GY ~ 1,
#             random = ~ re(~ stack(~ lk(Gen, K1) + lk(Gen, K2)), vars = ~ 1),
#             residual = ~ units,
#             data = d)

# Trying with 10 chromosome specific kinship matrices
M <- read.csv("~/repos/kernels_2026/DROPS/raw_data/7a-Genotyping_50K_41722.csv")
M.info <- read.csv("~/repos/kernels_2026/DROPS/raw_data/7b-InfoSNP_50K_41722.csv")
M.info$SNP.names <- gsub("-", ".", M.info$SNP.names)
SNPs <- vector("list", 10)
for (chr in 1:10) {
  SNPs[[chr]] <- M[, c("Ind", M.info[M.info$Chromosome == chr, "SNP.names"])]
}
sum(unlist(lapply(SNPs, ncol))) - 10 == ncol(M) - 1
gdatas <- Ks <- vector("list", 10)
for (chr in 1:10) {
  rownames(SNPs[[chr]]) <- SNPs[[chr]]$Ind
  SNPs[[chr]] <- SNPs[[chr]][, -1]
  gdatas[[chr]] <- statgenGWAS::createGData(geno = SNPs[[chr]])
  gdatas[[chr]] <- statgenGWAS::codeMarkers(gdatas[[chr]], removeDuplicates = T, MAF = 0.05)
  Ks[[chr]] <- statgenGWAS::kinship(gdatas[[chr]]$markers, "IBS")[levels(d$Gen), levels(d$Gen)]
}
# We filtered out some markers everywhere:
matrix(unlist(lapply(SNPs, function(SNP) ncol(SNP))) - unlist(lapply(gdatas, function(gdata) ncol(gdata$markers))))
# But still have plenty left:
matrix(unlist(lapply(gdatas, function(gdata) ncol(gdata$markers))))
# Ks for chr 6 and 10 are singular:
matrix(unlist(lapply(Ks, function(K) min(eigen(K)$values))))
Ks[[6]] <- Ks[[6]] + diag(1e-4, 246, 246)
Ks[[10]] <- Ks[[10]] + diag(1e-4, 246, 246)

system.time(
fit.asr <- asreml(fixed = GY ~ 1,
                  random = ~ vm(Gen, Ks[[1]]) + vm(Gen, Ks[[2]]) + vm(Gen, Ks[[3]]) + vm(Gen, Ks[[4]]) +
                    vm(Gen, Ks[[5]]) + vm(Gen, Ks[[6]]) + vm(Gen, Ks[[7]]) + vm(Gen, Ks[[8]]) + vm(Gen, Ks[[9]]) + vm(Gen, Ks[[10]]),
                  residual = ~ idv(units),
                  data = d))
summary(fit.asr)$varcomp

system.time(
fit.kmmr <- kmmr(fixed = GY ~ 1,
                 random = ~ re(~ stack(~ lk(Gen, Ks[[1]]) + lk(Gen, Ks[[2]]) + lk(Gen, Ks[[3]]) + 
                                         lk(Gen, Ks[[4]]) + lk(Gen, Ks[[5]]) + lk(Gen, Ks[[6]]) +
                                         lk(Gen, Ks[[7]]) + lk(Gen, Ks[[8]]) + lk(Gen, Ks[[9]]) +
                                         lk(Gen, Ks[[10]])), vars = ~ 1),
                 residual = ~ units,
                 data = d,
                 tolConvLL = 1e-9))

u.asr <- matrix(fit.asr$coefficients$random, nrow = 246)
u.asr <- rowSums(u.asr)
u.kmmr <- fit.kmmr$u
plot(u.kmmr, u.asr); abline(0, 1, col = "red")

plot(summary(fit.asr)$varcomp$component[1:11],
     summary(fit.kmmr)$varcomp$Estimate)
abline(0, 1, col = "red")

summary(fit.asr)$varcomp$component[1:11] - summary(fit.kmmr)$varcomp$Estimate


fit.1K <- asreml(fixed = GY ~ 1,
                 random = ~ vm(Gen, K1),
                 residual = ~ idv(units),
                 data = d)
u.1K <- fit.1K$coefficients$random
plot(u.kmmr, u.1K); abline(0, 1, col = "red")

# Taking 5 kinships
# system.time(
#   fit.asr <- asreml(fixed = GY ~ 1,
#                     random = ~ vm(Gen, Ks[[1]]) + vm(Gen, Ks[[3]]) + vm(Gen, Ks[[5]]) + vm(Gen, Ks[[7]]) +
#                       vm(Gen, Ks[[9]]),
#                     residual = ~ idv(units),
#                     data = d))
# summary(fit.asr)$varcomp
# 
# system.time(
#   fit.kmmr <- kmmr(fixed = GY ~ 1,
#                    random = ~ re(~ stack(~ lk(Gen, Ks[[1]]) + lk(Gen, Ks[[3]]) + lk(Gen, Ks[[5]]) + 
#                                            lk(Gen, Ks[[7]]) + lk(Gen, Ks[[9]])), vars = ~ 1),
#                    residual = ~ units,
#                    data = d))

