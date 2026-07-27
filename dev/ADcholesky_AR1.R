library(spam)

# Suppose we have dim plots in a single row, a column variance of 3, an autocorrelation of 0.99, and a residual variance of 2.
dim <- 1000
X <- matrix(1, dim, 1) # Fixed effect design matrix
Z <- diag(dim) # Random effect design matrix
W <- cbind(X, Z)
D <- abs(outer(1:dim, 1:dim, FUN = "-")) # Distance matrix of plots
theta <- c(3, 0.95, 2)

my_own_def <- function(theta) {
  # Making the inverse residual covariance matrix:
  R.inverse <- diag(theta[3]^-1, dim)
  
  # Making the inverse column AR1 covariance matrix:
  AR1.inverse <- matrix(0, dim, dim)
  diag(AR1.inverse) <- 1 + theta[2]^2 # Setting the main diagonal
  AR1.inverse[cbind(1:(dim - 1), 2:dim)] <- AR1.inverse[cbind(2:dim, 1:(dim - 1))] <- -theta[2] # Setting the super- and subdiagonals
  AR1.inverse[1, 1] <- AR1.inverse[dim, dim] <- 1 # Setting the corners
  AR1.inverse <- 1 / (1 - theta[2]^2) * AR1.inverse # Scaling the entire thing
  G.inverse <- theta[1]^-1 * AR1.inverse # Adding the inverse column variance
  
  # Making C:
  C <- t(W) %*% R.inverse %*% W
  C[2:(dim + 1), 2:(dim + 1)] <- C[2:(dim + 1), 2:(dim + 1)] + G.inverse
  
  # Derivative of C with respect to the residual variance:
  dCdtheta3 <- t(W) %*% diag(-theta[3]^-2, dim) %*% W
  
  # Derivative of C with respect to the column variance:
  dCdtheta1 <- rbind(cbind(matrix(0, 1, 1), matrix(0, 1, dim)),
                     cbind(matrix(0, dim, 1), -theta[1]^-2 * AR1.inverse))
  
  # Derivative of C with respect to the autocorrelation:
  # tmp is the explicit derivative of G^-1 with respect to the autocorrelation r.
  # If we compute dG^-1/dr as -G^-1 %*% dG/dr %*% G^-1 we get elements that are ~1e-16 which are not seen as zero by spam.
  # Explicit inverse can be checked using:
  # all.equal(tmp, -G.inverse %*% (theta[1] * D * theta[2]^(D - 1)) %*% G.inverse)
  tmp <- matrix(0, dim, dim)
  diag(tmp) <- (4 * theta[2]) / (theta[2]^2 - 1)^2
  tmp[cbind(1:(dim - 1), 2:dim)] <- tmp[cbind(2:dim, 1:(dim - 1))] <- (-theta[2]^2 - 1) / (theta[2]^2 - 1)^2
  tmp[1, 1] <- tmp[dim, dim] <- (2 * theta[2]) / (theta[2]^2 - 1)^2
  tmp <- tmp * theta[1]^-1
  dCdtheta2 <- rbind(cbind(matrix(0, 1, 1), matrix(0, 1, dim)),
                     cbind(matrix(0, dim, 1), tmp))
  
  return(list(C = as.spam(C),
              dC = list(as.spam(dCdtheta1),
                        as.spam(dCdtheta2),
                        as.spam(dCdtheta3))))
}

dlogdet_fun <- function(user_def, theta) {
  model_eval <- user_def(theta)
  Cinv <- solve(model_eval$C)
  dlogdet <- sapply(model_eval$dC, function(X) sum(diag(Cinv %*% X)))
}

dlogdet1 <- dlogdet_fun(user_def = my_own_def, theta = theta)

ADcholC_nl <- LMMsolver:::ADchol_nl(user_def = my_own_def, theta0 = theta)
dlogdet2 <- LMMsolver:::dlogdet(ADcholC_nl, theta)

all.equal(dlogdet1, as.numeric(dlogdet2)) # TRUE

microbenchmark::microbenchmark(
  standard = dlogdet_fun(my_own_def, theta),
  ADchol   = LMMsolver:::dlogdet(ADcholC_nl, theta),
  times = 10
)

test <- my_own_def(theta)
display(test$C[1:50, 1:50])
