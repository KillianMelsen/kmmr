# load("~/repos/kmmr/dev/exfa_manual_check.RData")
# Then run the solver prelims
# Then run a single iteration up to and including delta
# Now check against manually computed quantities and matrices:

ymat <- as.matrix(y)
Xmat <- as.matrix(X)
Zmat <- as.matrix(Z[[1]])
Wmat <- cbind(Xmat, Zmat)
Kmat <- as.matrix(solve(Ki))
Rimat <- as.matrix(Ri)
Cmat <- t(Wmat) %*% Rimat %*% Wmat
Rmat <- solve(Rimat)
Lmat <- matrix(init[15:56], 14, 3)
PSImat <- diag(1 - diag(Lmat %*% t(Lmat)))
cormat <- rbind(cbind(Lmat %*% t(Lmat) + PSImat,  -Lmat),
                cbind(-t(Lmat),                    diag(ncol(Lmat))))
varsmat <- outer(sqrt(c(init[1:nrow(Lmat)], rep(1, ncol(Lmat)))),
                 sqrt(c(init[1:nrow(Lmat)], rep(1, ncol(Lmat)))))
covmat <- varsmat * cormat
covmati <- solve(covmat)
Gmat <- kronecker(covmat, Kmat)
Gimat <- solve(Gmat)
Cmat[15:nrow(Cmat), 15:ncol(Cmat)] <- Cmat[15:nrow(Cmat), 15:ncol(Cmat)] + Gimat
Vmat <- Rmat + Zmat %*% Gmat %*% t(Zmat)
Vimat <- solve(Vmat)
Pmat <- Vimat - Vimat %*% Xmat %*% solve(t(Xmat) %*% Vimat %*% Xmat) %*% t(Xmat) %*% Vimat

yPy.check <- as.numeric(t(ymat) %*% Pmat %*% ymat)
logDetR.check <- sum(log(eigen(Rmat)$values))
logDetG.check <- sum(log(eigen(Gmat)$values))
logDetC.check <- sum(log(eigen(Cmat)$values))

all.equal(yPy, yPy.check, check.attributes = FALSE)
all.equal(logDetR, logDetR.check, check.attributes = FALSE)
all.equal(logDetG, logDetG.check, check.attributes = FALSE)
all.equal(logDetC, logDetC.check, check.attributes = FALSE)

Cimat <- solve(Cmat)
bumat <- Cimat %*% t(Wmat) %*% Rimat %*% ymat
umat <- bumat[15:nrow(bumat), 1, drop = FALSE]

all.equal(bumat, as.matrix(bu), check.attributes = FALSE)
all.equal(Cimat[as.matrix(C) != 0], as.matrix(C)[as.matrix(C) != 0], check.attributes = FALSE)

wv15 <- as.numeric(Zmat %*% kronecker(as.matrix(derivs[[1]][[15]]), Kmat) %*% Gimat %*% umat)
all.equal(wv15, Wu[,15])

wv51 <- as.numeric(Zmat %*% kronecker(as.matrix(derivs[[1]][[51]]), Kmat) %*% Gimat %*% umat)
all.equal(wv51, Wu[,51])

wv12 <- as.numeric(Zmat %*% kronecker(as.matrix(derivs[[1]][[12]]), Kmat) %*% Gimat %*% umat)
all.equal(wv12, Wu[,12])

wv35 <- as.numeric(Zmat %*% kronecker(as.matrix(derivs[[1]][[35]]), Kmat) %*% Gimat %*% umat)
all.equal(wv35, Wu[,35])

# Manual derivative:
dLmat <- matrix(0, nrow(Lmat), ncol(Lmat))
dLmat[35 - 14] <- 1
deriv <- rbind(cbind(Lmat %*% t(dLmat) + dLmat %*% t(Lmat),  -dLmat),
               cbind(-t(dLmat),                              matrix(0, ncol(Lmat), ncol(Lmat))))
diag(deriv) <- 0
wv35 <- as.numeric(Zmat %*% kronecker(varsmat * deriv, Kmat) %*% Gimat %*% umat)
all.equal(wv35, Wu[,35])

# Residual working variate:
emat <- ymat - Wmat %*% bumat
all.equal(emat, as.matrix(e))
wv57 <- as.numeric(diag(nrow(Rimat)) %*% Rimat %*% emat)
all.equal(wv57, Wu[,57])

AImat <- t(as.matrix(Wu)) %*% Pmat %*% as.matrix(Wu) / 2
all.equal(as.matrix(AI), AImat, check.attributes = FALSE)
CimatZZ <- Cimat[15:nrow(Cimat), 15:ncol(Cimat)]

s1 <- as.numeric(-0.5 * (sum(diag(kronecker(derivs[[1]][[1]], Kmat) %*% Gimat)) - sum(diag(Gimat %*% kronecker(derivs[[1]][[1]], Kmat) %*% Gimat %*% CimatZZ)) - (t(ymat) %*% Pmat %*% matrix(Wu[, 1]))))
all.equal(s1, scores[1])

s35 <- as.numeric(-0.5 * (sum(diag(kronecker(derivs[[1]][[35]], Kmat) %*% Gimat)) - sum(diag(Gimat %*% kronecker(derivs[[1]][[35]], Kmat) %*% Gimat %*% CimatZZ)) - (t(ymat) %*% Pmat %*% matrix(Wu[, 35]))))
all.equal(s35, scores[35])

# Residual score:
s57 <- as.numeric(-0.5 * (sum(diag(diag(nrow(Rimat)) %*% Rimat)) - sum(diag(t(Wmat) %*% Rimat %*% diag(nrow(Rimat)) %*% Rimat %*% Wmat %*% Cimat)) - (t(ymat) %*% Pmat %*% matrix(Wu[, 57]))))
all.equal(s57, scores[57])





