# Embedding and decomposition stages of univariate functional singular spectrum analysis
ufssa <- function(Y, L) {
  N <- ncol(Y@coefs[[1]])
  basis <- Y@B[[1]]
  d <- ncol(Y@B[[1]])
  grid <- Y@grid[[1]]
  K <- N - L + 1L
  if(ncol(Y@grid[[1]])==1){
    C_tilde <-t(onedG(A=basis%*%Y@coefs[[1]],B=basis,grid=Y@grid[[1]]))
    G <- onedG(A=basis,B=basis,grid=Y@grid[[1]])
  }else{
    C_tilde <-t(twodG(A=basis%*%Y@coefs[[1]],B=basis,grid=Y@grid[[1]]))
    G <- twodG(A=basis,B=basis,grid=Y@grid[[1]])
  }
  S0 <- SS(K, L, C_tilde, d)
  H <- solve(Gram(K, L, G, d))
  Q <- eigen(H %*% S0)
  Q$vectors <- Re(Q$vectors)
  out <- list(NA)
  d1 <- sum(Re(Q$values) > 0.00001)
  for (i in 1L:d1) out[[i]] <- Y@B[[1]]%*%Cofmat(d, L, Q$vectors[, i])
  out$values <- Re(Q$values[1L:d1])
  out$L <- L
  out$N <- N
  out$Y <- Y
  out$RVectrs <- uV(out,d1)
  return(out)
}

