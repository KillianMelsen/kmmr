library(Matrix)

# Observation-level residuals ====
if ("weights" %in% names(mod.args)) { # If we're using weighting
  
  R_diag_vals <- 1 / data_obs[[mod.args$weights]] # Weights are reciprocals of variances.
  
} else if (grepl("dsum\\(", as.character(mod.args$residual)[2])) { # If we fit heterogeneous residual variances... 
  
  vc$term <- sub(paste0("^", trial, "_(.*?)!R$"), "\\1", rownames(vc)) # Trial names.
  resid_comp <- vc$component[match(data_obs[[trial]], vc$term)] # Putting residual variances in the correct order and expanding to the full dataset.
  R_diag_vals <- resid_comp
  
} else { # If we fit homog. residual variance
  
  resid_comp <- vc$component[grepl("units!R", rownames(vc))]
  R_diag_vals <- rep(resid_comp, N) # Expanding to the full dataset.
  
}

R_full <- Diagonal(n = N, x = R_diag_vals)
R_inv <- Diagonal(x = 1 / diag(R_full)) # solve(R_full)

# Construct C_full (sparse) and inverse ====
p <- ncol(X_full)
q <- ncol(Z_full)

XtRinv <- t(X_full) %*% R_inv # p × N
ZtRinv <- t(Z_full) %*% R_inv # q × N
XtRinvX <- XtRinv %*% X_full  # p × p <-- C11
# Alternatively:
# XtRinvX <- crossprod(X_full, R_inv %*% X_full)

k_val <- kappa(as.matrix(XtRinvX), exact = FALSE)
rank_val <- qr(as.matrix(XtRinvX))$rank
if(rank_val < p) {
  warning("C11 is rank deficient! The ECs are perfectly collinear. The inverse will be unstable.")
}
if(k_val > 1e8) {
  message("High ill-conditioning detected. Applying ridge penalty to stabilize...")
  diag(XtRinvX) <- diag(XtRinvX) + (mean(diag(XtRinvX)) * 1e-4)
}

XtRinvZ <- XtRinv %*% Z_full # p × q <-- C12
ZtRinvX <- t(XtRinvZ)        # q × p <-- C21

# Computing Zt * Rinv * Z (C22) is inherently very costly;
# the code below should work (?) provided Rinv is diagonal!!
# Currently, abort if not. Consider fitting different model, e.g. w/o covariance between traits?

if(!all(R_inv == diag(diag(R_inv)))) stop("Non-diagonal R_inv matrix not implemented for computing coefficient matrix.")

w <- diag(R_inv)
Z_scaled <- Z_full
Z_scaled@x <- Z_scaled@x * w[Z_scaled@i + 1]
ZtRinvZ <- crossprod(Z_scaled, Z_full)

# Now get C and invert it
top <- cbind(XtRinvX, XtRinvZ)
bottom <- cbind(ZtRinvX, ZtRinvZ + G_inv_full)
C_full <- rbind(top, bottom)

# If C not (near-)symmetrical, something's wrong.
if(!isSymmetric(C_full, tol = tolCsymm)) {
  stop(paste0("Coefficient matrix (C) not symmetric at ",tolCsymm," tolerance. Something wrong with construction of C?"))
}

Cchol <- Cholesky(C_full, LDL = FALSE) # Is C_full sparse?
C_inv <- solve(Cchol, Diagonal(nrow(C_full)))

# If the goal is to make predicted values and associated prediction error variances:
# Every possible prediction is a linear combination of the BLUEs and BLUPs.
# Assume that bu is a column vector of [BLUEs^T BLUPs^T]^T
# Then any prediction can be written as k^T %*% bu, where k is a prediction vector
# containing the coefficients of the linear combination.
# Say that bu = [E1
#                E2
#                G1
#                G2
#                E1G1
#                E1G2
#                E2G1
#                E2G2],
# where Ei is the BLUE for the main effect of environment i
#       Gj is the BLUP for the main effect of genotype j
#       EiGj is the BLUP for the interaction of environment i and genotype j
#
# Say we wanted the predicted value for genotype 2 in environment 1.
# Then k = [1
#           0
#           0
#           1
#           0
#           1
#           0
#           0]
# I believe this should work the same if covariates are involved in the fixed or random effects.
pval <- crossprod(k, bu)

# Now we need var(pval) = var(k^T %*% bu)
# Note that var(k^T %*% bu) = k^T %*% var(bu) %*% k = k^T %*% C^-1 %*% k,
# see equation (313), The Matrix Cookbook, Version: November 15, 2012
#
# We already have the Cholesky factor of C, so all we need to do is solve the system
# C %*% x = k,
# which can be done very quickly using:
x <- solve(Cchol, k)
# Since C %*% x = k,
# x = C^-1 %*% k, so k^T %*% C^-1 %*% k = k^T %*% x:
pev <- crossprod(k, x)







