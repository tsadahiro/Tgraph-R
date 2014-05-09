K <- matrix(c(1,0,1, 1,1,-1, 0,1,1),ncol=3, byrow=T)
Kd <- K %*% diag(as.vector(solve(K) %*% matrix(c(1,0,0),ncol=1)))
Kt <- diag(as.vector(matrix(c(1,1i,-1),nrow=1) %*% solve(Kd)) ) %*% Kd
