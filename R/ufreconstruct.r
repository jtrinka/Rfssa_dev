# Reconstruction stage (including Hankelization) of univariate functional singular spectrum analysis

ufreconstruct <- function(U, groups = as.list(1L:10L)) {
  N <- ncol(U$Y@C[[1]])
  Y <- U$Y
  d <- ncol(U$Y@B[[1]])
  L <- U$L
  K <- N - L + 1L
  basis <- Y@B[[1]]
  m <- length(groups)
  out <- list()
  for (i in 1L:m) {
    Cx <- matrix(NA, nrow = d, ncol = N)
    g <- groups[[i]]
    S <- 0L
    for (j in 1L:length(g)) S <- S + ufproj(U, g[j], d)
    S <- fH(S, d)
    Cx[, 1L:L] <- S[, 1L, ]
    Cx[, L:N] <- S[, ,L]
    out[[i]] <- fts(list(basis%*%Cx),list(Y@B[[1]]),list(Y@grid[[1]]))
  }
  out$values <- sqrt(U$values)
  return(out)
}
