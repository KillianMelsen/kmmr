library(spam)
library(LMMsolver)

# small example with four supernodes:
n1 <- 2
n2 <- 5
n <- n1*n2
D1 <- diff.spam(diag.spam(n1), diff = 1)
D2 <- diff.spam(diag.spam(n2), diff = 1)

A <- diag.spam(1, n)
B1 <- diag.spam(1, n1) + t(D1)%*%D1
B2 <- diag.spam(1, n2) + t(D2)%*%D2
B <- kronecker(B1, B2)

# simple linear example, returns matrix of class spam:
f <- function(theta) { return(theta[1]*A + theta[2]*B)}

# partial derivatives, in this case linear
dC1 <- A
dC2 <- B

# init C0
theta0 <- c(1.0, 1.0)
C0 <- f(theta0)
class(C0)
display(C0)
display(t(chol(C0)))

# construct the object, using symbolic factorization
obj <- LMMsolver:::SparseCholesky(C0)
class(obj)

#
# update the object (in practice repeated many times),
# with internally the F and L objects as in Smith1995
#
theta <- c(1, 2)
C <- f(theta)
obj <- LMMsolver:::update(obj, C)

g <- vector(mode = "numeric", length = 2)
g[1] <- LMMsolver:::dlogdet(obj, dC1)
g[2] <- LMMsolver:::dlogdet(obj, dC2)
g

# a simple check, should be equal to n:
sum(theta*g) - n

# some other functions,
# all (for moment) internal:
set.seed(1234)
b <- rnorm(n)
x1 <- solve(C, b)
x2 <- LMMsolver:::solve(obj, b)
all.equal(x1, x2)

LMMsolver:::logdet(obj)
as.numeric(determinant(C)$modulus)



