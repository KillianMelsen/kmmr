# # Implement the ar1 correlation model
# K <- readRDS("dev/exampledata_K.rds")
# d <- readRDS("~/repos/kmmr/dev/exampledata_singleStage.rds")
# d <- droplevels(d[d$type == "CrossValidation",])
# d <- droplevels(d[d$Variety_ID %in% rownames(K),])
# d <- droplevels(d[d$Experiment != "Deb13R",])
# d$Env <- substr(d$Experiment, 1, 5)
# d$Man <- substr(d$Experiment, 6, 6)
# d$Gen <- d$Variety_ID
# d$GY <- d$grain.yield
# d$Rep <- d$Replicate
# d$Block <- d$block
# d$Plot <- d$plotID
# d$Row <- d$Row
# d$Col <- d$Column
# 
# d <- d[, c("Plot", "Row", "Col", "Block", "Rep", "Man", "Env", "Gen", "GY")]
# d <- na.omit(d)
# d$Rep <- factor(d$Rep, sort(unique(d$Rep)), sort(unique(d$Rep)))
# d$Row <- factor(d$Row, sort(unique(d$Row)), sort(unique(d$Row)))
# d$Col <- factor(d$Col, sort(unique(d$Col)), sort(unique(d$Col)))
# d$Block <- factor(d$Block, sort(unique(d$Block)), sort(unique(d$Block)))
# d$Gen <- as.factor(d$Gen)
# 
# dss <- droplevels(d[d$Man == "R" & d$Env == "Bol12",])

data(shf)
fit.asr1 <- asreml(yield ~ -1 + Variety,
                  random = ~ ar1v(Row):ar1(Column),
                  residual = ~ idv(units),
                  data = shf,
                  maxit = 25)

fit.kmmr1 <- kmmr(yield ~ -1 + Variety,
                 random = ~ re(~ ar1(Row):ar1(Column), vars = ~ 1),
                 residual = ~ units,
                 data = shf,
                 sparse = TRUE)
cor(fit.kmmr1$u@x, fit.asr1$coefficients$random[,])
# ^^ Pretty much the same, so good.

# Now we set all data from row 4 to NA:
shf2 <- shf
shf2$yield[shf2$Row == 4] <- NA
fit.asr2 <- asreml(yield ~ -1 + Variety,
                   random = ~ ar1v(Row):ar1(Column),
                   residual = ~ idv(units),
                   data = shf2,
                   maxit = 25)

fit.kmmr2 <- kmmr(yield ~ -1 + Variety,
                  random = ~ re(~ ar1(Row):ar1(Column), vars = ~ 1),
                  residual = ~ units,
                  data = shf2,
                  sparse = TRUE)
cor(fit.kmmr2$u@x, fit.asr2$coefficients$random[,])
# ^^ Pretty much the same, so good.

# This is where asreml seems to have an issue. When we drop unused levels on the data
# we should lose the level corresponding to Row == 4 because there is no data from row 4.
# But completely dropping that level from the data should not result in a different fit.
# It should result in an identical fit to the #2 models above.
shf3 <- droplevels(na.omit(shf2))
fit.asr3 <- asreml(yield ~ -1 + Variety,
                   random = ~ ar1v(Row):ar1(Column),
                   residual = ~ idv(units),
                   data = shf3,
                   maxit = 25)

fit.kmmr3 <- kmmr(yield ~ -1 + Variety,
                  random = ~ re(~ ar1(Row):ar1(Column), vars = ~ 1),
                  residual = ~ units,
                  data = shf3,
                  sparse = TRUE)
cor(fit.kmmr3$u@x, fit.asr3$coefficients$random[,])
all.equal(fit.kmmr2$k, fit.kmmr3$k) # TRUE.

# Asreml seems to construct the distance matrix simply using 1:length(levels(d$Row)),
# ignoring the fact that there is now a gap between rows...
# We can't seem to fix this either with asreml by adding back the level for Row == 4:
shf4 <- shf3
shf4$Row <- factor(as.character(shf4$Row), as.character(min(as.numeric(levels(shf4$Row))):max(as.numeric(levels(shf4$Row)))))
fit.asr4 <- asreml(yield ~ -1 + Variety,
                   random = ~ ar1v(Row):ar1(Column),
                   residual = ~ idv(units),
                   data = shf4,
                   maxit = 25)

fit.kmmr4 <- kmmr(yield ~ -1 + Variety,
                  random = ~ re(~ ar1(Row):ar1(Column), vars = ~ 1),
                  residual = ~ units,
                  data = shf4,
                  sparse = TRUE)
cor(fit.kmmr4$u@x, fit.asr4$coefficients$random[,]) # Error, incompatible dimensions. Seems like asreml just discards levels for which there is no data.
all.equal(fit.kmmr2$k, fit.kmmr4$k) # TRUE.

# Let's try adding a single record with Row == 4 and missing GY:
shf5 <- rbind(shf4,
              data.frame(Rep = 1,
                         RowBlk = 1,
                         ColBlk = 1,
                         Row = 4,
                         Column = 1,
                         Variety = 1,
                         yield = NA))
fit.asr5 <- asreml(yield ~ -1 + Variety,
                   random = ~ ar1v(Row):ar1(Column),
                   residual = ~ idv(units),
                   data = shf5,
                   maxit = 25)

fit.kmmr5 <- kmmr(yield ~ -1 + Variety,
                  random = ~ re(~ ar1(Row):ar1(Column), vars = ~ 1),
                  residual = ~ units,
                  data = shf5,
                  sparse = TRUE)
cor(fit.kmmr5$u@x, fit.asr5$coefficients$random[,]) 
all.equal(fit.kmmr2$k, fit.kmmr5$k) # TRUE.
# Correlation is 1 again, so asreml needs to have at least a single record in the dataframe
# for each level of the Row or Column. Otherwise it thinks the distance between Row 3 and Row 5 is a single unit instead of 2.
# Try to verify that:
shf6 <- shf3
# The line below changes rows to 1:9 with levels 1:9. Row 4 completely disappears.
shf6$Row <- factor(as.numeric(shf6$Row), unique(as.numeric(shf6$Row)))
fit.kmmr6 <- kmmr(yield ~ -1 + Variety,
                  random = ~ re(~ ar1(Row):ar1(Column), vars = ~ 1),
                  residual = ~ units,
                  data = shf6,
                  sparse = TRUE)
cor(fit.kmmr6$u@x, fit.asr3$coefficients$random[,]) # Correlation of 1 with the wrong asreml fit, so indeed that's what happens.





