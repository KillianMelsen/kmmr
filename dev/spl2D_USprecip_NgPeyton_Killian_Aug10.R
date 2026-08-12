library(LMMsolver)
library(spam)

# Get precipitation data from spam
data(USprecip)
dat <- data.frame(USprecip)

# only use observed data
dat <- subset(dat, infill==1)
nrow(dat) # 5906 true records, as in SAP2014 paper.

y <- dat$anomaly
x1 <- dat$lon
x2 <- dat$lat

# set parameters:
innerKnots <- 50 # See SAP 2014 paper
knots <- nseg <- c(innerKnots+1, innerKnots+1)

x1lim <- c(min(x1)-0.01, max(x1)+0.01)
x2lim <- c(min(x2)-0.01, max(x2)+0.01)

#
# use spatial option in LMMsolve:
#
# s <- proc.time()[3]
# obj1 <- LMMsolve(fixed = anomaly~1,
#                  spline = ~spl2D(x1 = lon, x2 = lat, nseg = nseg,
#                                   x1lim=x1lim, x2lim=x2lim),
#                  data = dat)
# e <- proc.time()[3]
# cat("Time LMMsolve ", e-s, " seconds\n")
# EDtbl <- summary(obj1)

knots <- list()
knots[[1]] <- LMMsolver:::PsplinesKnots(x1lim[1], x1lim[2], degree = 3, nseg = nseg[1])
knots[[2]] <- LMMsolver:::PsplinesKnots(x2lim[1], x2lim[2], degree = 3, nseg = nseg[2])

x1 <- dat$lon
x2 <- dat$lat
B1 <- LMMsolver:::Bsplines(knots[[1]], x1)
B2 <- LMMsolver:::Bsplines(knots[[2]], x2)
q <- c(ncol(B1), ncol(B2))

B <- LMMsolver:::RowKronecker(B1, B2)

dx1 <- attr(knots[[1]], which = 'dx')
dx2 <- attr(knots[[2]], which = 'dx')

P1 <- LMMsolver:::constructPenalty(q[1],2) %x% diag(1,q[2])
P2 <- diag(1,q[1]) %x% LMMsolver:::constructPenalty(q[2],2)


lC <- list()
lC[[1]] <- as.spam(crossprod(B))
lC[[2]] <- as.spam(P1)
lC[[3]] <- as.spam(P2)

C <- Reduce('+', lC)
# display(C)
# display(t(chol(C)))
summary(C)
summary(chol(C))

obj <- LMMsolver:::SparseCholesky(C)
CMatrix <- Matrix(as.matrix(C))

microbenchmark::microbenchmark(obj1 <- LMMsolver:::SparseCholesky(C),
                               obj <- LMMsolver:::update(obj, C),
                               sis <- SIS(Cholesky(CMatrix, LDL = FALSE, perm = TRUE, super = TRUE)),
                               times = 10)

Fmat <- chol(C) # This is purely to get a dummy object into which we can put obj@ADentries to make the F matrix
Fmat@entries <- obj@ADentries
Fmat <- as.matrix(Fmat)[obj@invpivot + 1, obj@invpivot + 1] # Converting to a standard matrix and permuting back

sum(Fmat != 0)
length(sis@x)

# Converting to a normal matrix:
sis <- as.matrix(sis)

# "Folding over" the F matrix:
Fmat[upper.tri(Fmat)] <- Fmat[upper.tri(Fmat, diag = FALSE)] / 2
Fmat <- Fmat + t(Fmat) - diag(diag(Fmat))

# sis has some more infill:
prod(dim(Fmat))
sum(Fmat != 0)
sum(sis != 0)

# Overlapping nonzeros:
nz1 <- which(Fmat != 0)
nz2 <- which(sis != 0)
nz <- intersect(nz1, nz2)
all.equal(Fmat[nz], sis[nz])
all.equal(solve(as.matrix(C))[nz1], Fmat[nz1]) # Not TRUE, don't know why
all.equal(solve(as.matrix(C))[nz2], sis[nz2]) # TRUE






