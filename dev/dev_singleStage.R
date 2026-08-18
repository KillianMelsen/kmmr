fit <- kmmr(fixed = GY ~ -1 + Env,
            random = ~ atEach("Env", ~ re(~ id(Rep), vars = ~ 1), levels = levels(d$Env)) +
                       atEach("Env", ~ re(~ ar1(Row):ar1(Col), vars = ~ 1), levels = levels(d$Env)) +
                       re(~ exfa(Env, 3):lk(Gen, Ki), vars = ~ Env),
            residual = ~ re(~ id(Env):id(Gen):id(Rep), vars = ~ 1),
            data = d,
            sparse = TRUE)