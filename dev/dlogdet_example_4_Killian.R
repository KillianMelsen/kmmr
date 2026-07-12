library(spam)
library(LMMsolver)
library(microbenchmark)

trace <- function(A) {sum(diag(A))}

# splines example, just to get a sparse matrix C:
nseg <- 500
deg <- 3
n <- 1000
x <- seq(0, 1, length = n)
lambda <- 0.5

knots <- LMMsolver:::PsplinesKnots(0, 1, 3, nseg)
B <- LMMsolver:::Bsplines(knots, x)
q <- ncol(B)

D <- diff(diag.spam(q), diff = 2)
BtB <- crossprod(B)
DtD <- crossprod(D)

# define the matrix C, plus the partial derivatives:
my_own_def <- function(theta)
{
  C <- theta[1] * BtB + theta[2] * DtD
  list(
    C  = C,
    dC = list(BtB, DtD)
  )
}

theta <- c(1, lambda)
test <- my_own_def(theta)
display(test$C[1:20, 1:20])
dim(test$C)

# standard expression using the inverse of C to calculate the
# partial derivatives of logdet-C:
dlogdet_fun <- function(user_def, theta) {
  model_eval <- user_def(theta)
  Cinv <- solve(model_eval$C)
  dlogdet <- sapply(model_eval$dC, function(X) trace(Cinv %*% X))
}

dlogdet1 <- dlogdet_fun(user_def = my_own_def, theta = theta)

# using LMMsolver internal functions using Automated Differentation:
# 1) constructor of object
ADcholC_nl <- LMMsolver:::ADchol_nl(user_def = my_own_def, theta0 = theta)
# 2) calculates the partial derivatives, plus logdet as attr.
dlogdet2 <- LMMsolver:::dlogdet(ADcholC_nl, theta)

# compare results
dlogdet1
as.numeric(dlogdet2)

microbenchmark(
  standard = dlogdet_fun(my_own_def, theta),
  ADchol   = LMMsolver:::dlogdet(ADcholC_nl, theta),
  times = 100
)

