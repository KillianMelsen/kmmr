resvar <- 1.5 # Residual variance
pn <- c(2, 3, 4, 5) # Plot numbers to be used as relative weights
Omega <- diag(resvar / pn) # Covariance matrix
wt <- diag(solve(Omega)) # Absolute weights

# Situation 1, using relative weights, and assuming we estimate a residual variance of 1.5:
Ri <- solve(1.5 * diag(4))
W <- diag(sqrt(pn))
solve(W %*% Ri %*% W)

# Situation 2, using absolute weights and a fixed residual variance of 1.0:
Ri <- solve(1.0 * diag(4))
W <- diag(sqrt(wt))
solve(W %*% Ri %*% W)

Omega





resvar <- 1.5 # Residual variance
pn <- c(2, 4, 2, 4) # Plot numbers to be used as relative weights
Omega <- diag(sqrt(resvar / pn)) %*% kronecker(matrix(c(1.0, 0.2, 0.2, 1.0), 2, 2), diag(2)) %*% diag(sqrt(resvar / pn)) # Covariance matrix
wt <- diag(solve(Omega)) # Absolute weights

# Situation 1, using relative weights, and assuming we estimate a residual variance of 1.5:
Ri <- solve(1.5 * kronecker(matrix(c(1.0, 0.2, 0.2, 1.0), 2, 2), diag(2)))
W <- diag(sqrt(pn))
solve(W %*% Ri %*% W)

# Situation 2, using absolute weights and a fixed residual variance of 1.0:
Ri <- solve(1.0 * kronecker(matrix(c(1.0, 0.2, 0.2, 1.0), 2, 2), diag(2)))
W <- diag(sqrt(wt))
solve(W %*% Ri %*% W)

Omega



