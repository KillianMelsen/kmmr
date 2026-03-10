# Seed:
set.seed(1997)

# Loading example data (DROPS):
d <- readRDS("dev/exampledata.rds")[, c("Gen", "Man", "Env", "GY")]
ED <- readRDS("dev/exampledata_ED.rds")
EC <- readRDS("dev/exampledata_EC.rds")
GD <- readRDS("dev/exampledata_GD.rds")
K <- readRDS("dev/exampledata_K.rds")
Ki <- solve(K)
attr(Ki, "INVERSE") <- TRUE

# Adding a ManEnv factor for ASReml/cornfruit:
d$ManEnv <- factor(paste(as.character(d$Man), as.character(d$Env), sep = ":"),
                   levels = c(paste0("R:", levels(d$Env)),
                              paste0("W:", levels(d$Env))))

# Single variance Gaussian kernel:
fit <- kmmr(fixed = GY ~ -1 + Man:Env,
            random = ~ re(~ corg(Man):gk(Env, base = ED):lk(Gen, base = Ki), vars = ~ Man),
            residual = ~ units,
            data = d,
            tolConvLL = 0.01)

# Checking match with ASReml/cornfruit:
vf <- cornfruit::svgk(D = ED[levels(d$Env), levels(d$Env)])
init <- rep(0.1, 4)
type <- c("V", "V", "R", "V")
cons <- c("P", "P", "U", "P")
fit <- asreml::asreml(fixed = GY ~ -1 + ManEnv,
                      random = ~ own(ManEnv, "vf", init, type, cons):vm(Gen, Ki),
                      residual = ~ units,
                      data = d)