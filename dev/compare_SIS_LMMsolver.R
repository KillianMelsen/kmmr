library(spam)
library(LMMsolver)

# DROPS coefficient matrix of an expanded FA structure for 14 environments and
# dense kinship for 246 genotypes:
C <- readRDS("C_E14_M2_G246.rds")

system.time(Ci.naive <- solve(as.matrix(C))) # Naive inverse, ~ 2.8 seconds

# Code below uses the standard Matrix::Cholesky function for the supernodal factorization.
# Then it uses a slightly modified cholmod_extra C program to compute the sparse inverse subset.
# The C code already takes care of permuting it back.
system.time(Cinv <- kmmr::SIS(Matrix::Cholesky(C, LDL = FALSE, super = TRUE, perm = TRUE))) # ~ 0.16 seconds

# Converting to spam for LMMsolver
Cspam <- spam::as.spam(as.matrix(C))
# class(Cspam)
# display(Cspam)
# display(t(chol(Cspam)))
system.time(obj <- LMMsolver:::SparseCholesky(Cspam)) # Initial factorization, ~ 2.6 seconds
system.time(obj <- LMMsolver:::update(obj, Cspam)) # What would be used if the variance components changed, ~ 2.2 seconds


Fmat <- chol(Cspam) # This is purely to get a dummy object into which we can put obj@ADentries to make the F matrix
Fmat@entries <- obj@ADentries
Fmat <- as.matrix(Fmat)[obj@invpivot + 1, obj@invpivot + 1] # Converting to a standard matrix and permuting back

sum(Fmat != 0) # 2184002
length(Cinv@x) # 2184002, interestingly, if I set perm = TRUE on line 13 I get slightly more nonzeros (2190644)...

# Converting to a normal matrix:
Cinv <- as.matrix(Cinv)

# "Folding over" the F matrix:
Fmat[upper.tri(Fmat)] <- Fmat[upper.tri(Fmat, diag = FALSE)] / 2
Fmat <- Fmat + t(Fmat) - diag(diag(Fmat))

all.equal(Cinv, Fmat) # TRUE




