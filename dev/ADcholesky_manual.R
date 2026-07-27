# C <- matrix(c(x^3,       2x + x^2,
#               2x + x^2,  x^4),
#             2, 2)
# x = 2.15:
C <- matrix(c(2.15^3,             2 * 2.15 + 2.15^2,
              2 * 2.15 + 2.15^2,  2.15^4),
            2, 2)

# dCdx <- matrix(c(3x^2,    2 + 2x,
#                  2 + 2x,  4x^3),
#                2, 2)
dCdx <- matrix(c(3 * 2.15^2,    2 + 2 * 2.15,
                 2 + 2 * 2.15,  4 * 2.15^3),
               2, 2)

# Naive:
dlogdetCdx.naive <- sum(dCdx * solve(C))

# Cholesky algorithm:
L <- matrix(0, nrow(C), ncol(C))
L[lower.tri(L, diag = TRUE)] <- C[lower.tri(C, diag = TRUE)]
N <- ncol(L)
for (k in 1:N) {
  # Pivot:
  L[k, k] <- sqrt(L[k, k])
  
  if (k < N) {
    # Lead column:
    for (j in (k + 1):N) {
      L[j, k] <- L[j, k] / L[k, k]
    }
    
    # Row operations:
    for (j in (k + 1):N) {
      for (i in j:N) {
        L[i, j] <- L[i, j] - L[i, k] * L[j, k]
      }
    }
  }
}

# Backward differentiation:
WS <- matrix(0, nrow(L), ncol(L)) # Workspace
diag(WS) <- 2 / diag(L) # d f(L) / d L, f(L) = sum(2 * log(diag(L)))
for (k in N:1) {
  if (k < N) {
    # Row operations:
    for (j in (k + 1):N) {
      for (i in j:N) {
        WS[i, k] <- WS[i, k] - WS[i, j] * L[j, k]
        WS[j, k] <- WS[j, k] - WS[i, j] * L[i, k]
      }
    }
    
    # Lead column:
    for (j in (k + 1):N) {
      WS[j, k] <- WS[j, k] / L[k, k]
      WS[k, k] <- WS[k, k] - L[j, k] * WS[j, k]
    }
  }
  
  # Pivot:
  WS[k, k] <- 0.5 * WS[k, k] / L[k, k]
}

dlogdetCdx.AD <- sum(dCdx * WS)

all.equal(dlogdetCdx.naive, dlogdetCdx.AD) # TRUE

