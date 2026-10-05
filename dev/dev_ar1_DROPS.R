# Implement the ar1 correlation model
K <- readRDS("dev/exampledata_K.rds")
Ki.asr <- solve(K); attr(Ki.asr, "INVERSE") <- TRUE
Ki.kmmr <- solve(Matrix(K))
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

# dss <- droplevels(d[d$Man == "R" & d$Env == "Gai12",]) # Or Cra12, both fit fine in asreml and kmmr
# ggplot2::ggplot(dss, ggplot2::aes(x = as.numeric(as.character(Col)), y = as.numeric(as.character(Row)), shape = Rep, color = GY)) +
#   ggplot2::geom_point(size = 5) +
#   ggplot2::scale_color_continuous(palette = "viridis")
# 
# fit.asr <- asreml(GY ~ 1 + Gen,
#                   random = ~ ar1v(Row):ar1(Col),
#                   residual = ~ idv(units),
#                   data = dss,
#                   maxit = 25)
# 
# fit.kmmr <- kmmr(GY ~ 1 + Gen,
#                 random = ~ re(~ ar1(Row):ar1(Col), vars = ~ 1),
#                 residual = ~ units,
#                 data = dss,
#                 sparse = T,
#                 stepSize = c(0.1, 0.3, 0.5))
# 
# cor(fit.kmmr$u@x, fit.asr$coefficients$random[,])

# Now work on trial/environment specific random effects for single-stage analyses:
dss <- droplevels(d[d$Man == "R" & d$Env %in% c("Gai12", "Cra12", "Gai13", "Gra13"),]) # , "Cam12"
ggplot2::ggplot(dss, ggplot2::aes(x = as.numeric(as.character(Col)), y = as.numeric(as.character(Row)), shape = Rep, color = GY)) +
  ggplot2::facet_wrap(ggplot2::vars(Env)) +
  ggplot2::geom_point(size = 5) +
  ggplot2::scale_color_continuous(palette = "viridis")

# Rep is not really an appropriate name here as it doesn't identify a unique row for B73:
# dss[dss$Env == "Gai12" & dss$Gen == "B73" & dss$Rep == 1,]

# Instead, call it a complete block and call the other blocks incomplete:
colnames(dss)[colnames(dss) == "Rep"] <- "compBlock"
colnames(dss)[colnames(dss) == "Block"] <- "incompBlock"
# Then add an actual replicate factor that, for each Env and Gen, uniquely identifies a single record:
dss$Rep <- ave(numeric(nrow(dss)), dss$Env, dss$Gen, FUN = seq_along)
dss$Rep <- factor(dss$Rep, 1:max(dss$Rep))

# Complete blocks are named inconveniently:
dss$compBlock <- as.numeric(dss$compBlock)
dss$compBlock[dss$compBlock == 4] <- 2
dss$compBlock[dss$compBlock == 5] <- 1
dss$compBlock <- as.factor(dss$compBlock)

system.time(
fit.asr <- asreml(GY ~ 1 + Env,
                  random = ~ at(Env):idv(compBlock) +
                             at(Env):idv(compBlock):id(incompBlock) +
                             at(Env):ar1v(Row):ar1(Col) +
                             fa(Env, 1):vm(Gen, Ki.asr),
                  residual = ~ dsum(~ units | Env),
                  data = dss,
                  maxit = 30)
)

system.time(
  
fit.kmmr <- kmmr(fixed = GY ~ 1 + Env,
                 random = ~ atEach("Env", ~ re(~ id(compBlock), vars = ~ 1)) +
                            atEach("Env", ~ re(~ id(compBlock):id(incompBlock), vars = ~ 1)) +
                            atEach("Env", ~ re(~ ar1(Row):ar1(Col), vars = ~ 1)) +
                            re(~ exfa(Env, 1):lk(Gen, Ki.kmmr), vars = ~ Env),
                 residual = ~ re(~ id(Env):id(Gen):id(Rep), vars = ~ Env),
                 data = dss,
                 stepSize = c(0.1, 0.2, 0.3, 0.4, 0.5),
                 maxit = 30)

)

# The last (nEnv + m) x 246 BLUPs should be for the GxE part:
u.asr <- fit.asr$coefficients$random[(nrow(fit.asr$coefficients$random) - (length(levels(dss$Env)) * 246) + 1):(nrow(fit.asr$coefficients$random)), , drop = FALSE]
u.kmmr <- fit.kmmr$u[(nrow(fit.kmmr$u) - (length(levels(dss$Env)) * 246) + 1):nrow(fit.kmmr$u), , drop = FALSE]
rownames(u.asr) <- gsub("Env_(.*):vm\\(Gen, Ki.asr\\)_(.*)", "\\1:\\2", rownames(u.asr))
plot(u.asr, u.kmmr)
cor(u.asr[,], u.kmmr@x)

u.asr <- u.asr[order(u.asr[,], decreasing = TRUE), , drop = FALSE]
u.kmmr <- u.kmmr[order(u.kmmr@x, decreasing = TRUE), , drop = FALSE]
top.asr <- u.asr[1:floor(0.1 * nrow(u.asr)), , drop = FALSE]
top.kmmr <- u.kmmr[1:floor(0.1 * nrow(u.kmmr)), , drop = FALSE]
plot(top.asr[,], top.kmmr@x)




