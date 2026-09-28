# Old:
Cchol.old <- readRDS("dev/fix_updatesis/Cchol.old.rds")
Ci.old <- readRDS("dev/fix_updatesis/Ci.old.rds")
SC.old <- readRDS("dev/fix_updatesis/SC.old.rds")
Q.old <- readRDS("dev/fix_updatesis/Q.old.rds")

system.time(Ci.old <- SIS(Cchol.old))
# Ci.old is 1463 x 1463, Q.old is 1463 x 1, SC.old is 1 x 1, so W_mat is 1463 x 1
system.time(inplaceUpdateSIS(x_ptr = Ci.old@x, p_ptr = Ci.old@p, i_ptr = Ci.old@i, W_mat = as.matrix(Q.old %*% solve(SC.old)), Q = as.matrix(Q.old)))

# New:
Ccholss.new <- readRDS("dev/fix_updatesis/Ccholss.new.rds")
SCchol.new <- readRDS("dev/fix_updatesis/SCchol.new.rds")
Q.new <- readRDS("dev/fix_updatesis/Q.new.rds")

system.time(Ciss.new <- SIS(Ccholss.new))
system.time(Cidd.new <- chol2inv(as(SCchol.new, "dtrMatrix")))
# Ciss.new is 763 x 763, Q.new is 763 x 701, Cidd.new is 701 x 701, so W_mat is 763 x 701
system.time(W_mat.new <- as.matrix(Q.new %*% Cidd.new))
system.time(Q_mat.new <- as.matrix(Q.new))
Ciss.new2 <- Ciss.new * 1
Ciss.new3 <- Ciss.new * 1
system.time(inplaceUpdateSIS(x_ptr = Ciss.new@x, p_ptr = Ciss.new@p, i_ptr = Ciss.new@i, W_mat = W_mat.new, Q = Q_mat.new))
tW_mat.new <- t(W_mat.new)
tQ.new <- t(Q_mat.new)
system.time(inplaceUpdateSISnew(x_ptr = Ciss.new2@x, p_ptr = Ciss.new2@p, i_ptr = Ciss.new2@i, tW_mat = tW_mat.new, tQ = tQ.new))
system.time(Cisd.new <- -Q.new %*% Cidd.new)



system.time(inplaceUpdateSISnewnew(Ciss.new3, tcrossprod(Cidd.new, Q.new), t(Q.new)))
