# Univariate and multivariate weighted correlation used to find weighted correlation matrix for the grouping stage
# of ufssa and mfssa.

ufwcor <- function(U, groups) {
  if(is.numeric(groups)) groups <- as.list(groups)
  d <- length(groups)
  Q <- freconstruct(U, groups = groups)
  N <- U$N
  L <- U$L
  K <- N - L + 1L
  w <- 1L:N
  basis <- U$Y@B[[1]]
  if(ncol(U$Y@grid[[1]])==1){
    G <- onedG(A=basis,B=basis,grid=U$Y@grid[[1]])
  }else{
    G <- twodG(A=basis,B=basis,grid=U$Y@grid[[1]])
  }
  L1 <- min(L, K)
  K1 <- max(K, L)
  w[L1:K1] <- L1
  w[(K1 + 1L):N] <- N + 1L - ( (K1 +1L):N)
  out <- matrix(1L, nrow = d, ncol = d)
  for (i in 1L:(d - 1)) {
    for (j in (i + 1L):d) {
      out[i, j] <- winprod(Q[[i]]@C[[1]], Q[[j]]@C[[1]], w, G)/sqrt(winprod(Q[[i]]@C[[1]],Q[[i]]@C[[1]], w, G) * winprod(Q[[j]]@C[[1]],Q[[j]]@C[[1]], w, G))
    }
  }
  for (i in 2:d) for (j in 1:(i - 1)) out[i, j] <- out[j, i]
  return(out)
}

mfwcor <- function(U, groups) {
  if(is.numeric(groups)) groups <- as.list(groups)
  d <- length(groups)
  Q <- mfreconstruct(U, groups = groups)
  N <- U$N
  L <- U$L
  K <- N - L + 1L
  w <- 1L:N
  p = length(U$Y@C)
  Y = U$Y
  G = list()
  for(i in 1:p){
    if(ncol(Y@grid[[i]])==1){
      G[[i]] <- t(onedG(A=Y@B[[i]],B=Y@B[[i]],grid = Y@grid[[i]]))
    }else{
      G[[i]] <- t(twodG(A=Y@B[[i]],B=Y@B[[i]],grid = Y@grid[[i]]))

    }
  }
  L1 <- min(L, K)
  K1 <- max(K, L)
  w[L1:K1] <- L1
  w[(K1 + 1L):N] <- N + 1L - ( (K1 +1L):N)
  wcor <- matrix(1L, nrow = d, ncol = d)
  for (i in 1L:(d - 1)){
    Q_i=Q[[i]]
    Q_i_l <- list()
    for(k in 1:p){
      Q_i_l[[k]]=Q_i@C[[k]]
    }
    for (j in (i + 1L):d){
      Q_j=Q[[j]]
      Q_j_l <- list()
      for(k in 1:p){
        Q_j_l[[k]]=Q_j@C[[k]]
      }
      wcor[i, j] <- mwinprod(Q_i_l, Q_j_l, w, G, p)/sqrt(mwinprod(Q_i_l,Q_i_l, w, G, p) * mwinprod(Q_j_l,Q_j_l, w, G, p))
    }
  }
  for (i in 2:d) for (j in 1:(i - 1)) wcor[i, j] <- wcor[j, i]
  return(wcor)
}

#' Weighted Correlation Matrix
#'
#' This function returns the weighted correlation (w-correlation) matrix for functional time series (\code{\link{fts}}) objects
#' that were reconstructed from functional singular spectrum analysis (\code{\link{fssa}}) objects.
#' @return a square matrix of w-correlation values for the reconstructed \code{\link{fts}} objects that were built from
#' \code{\link{fssa}} components
#' @param U an object of class \code{\link{fssa}}
#' @param groups a list or vector of indices which determines the grouping used for the reconstruction
#' in pairwise w-correlations matrix
#' @examples
#'
#' \dontrun{
#'
#' ## Univariate W-Correlation Example on Callcenter data
#' data("Callcenter")
#' require(fda)
#' require(Rfssa)
#' ## Define functional objects
#' D <- matrix(sqrt(Callcenter$calls),nrow = 240)
#' N <- ncol(D)
#' time <- 1:N
#' K <- nrow(D)
#' u <- seq(0,K,length.out =K)
#' d <- 22 #Optimal Number of basis elements
#' basis <- create.bspline.basis(c(min(u),max(u)),d)
#' Ysmooth <- smooth.basis(u,D,basis)
#' ## Define functional time series
#' Y <- fts(Ysmooth$fd)
#' ## Decomposition stage of univariate functional singular spectrum analysis
#' L <- 28
#' U <- fssa(Y,L)
#' ufwcor=fwcor(U = U,groups = list(1,2,3))
#' wplot(W=ufwcor)
#'
#' ## Multivariate W-Correlation Example on Bivariate Satelite Image Data
#' require(fda)
#' require(Rfssa)
#' ## Raw image data
#' NDVI=Jambi$NDVI
#' EVI=Jambi$EVI
#' ## Kernel density estimation of pixel intensity
#' D0_NDVI <- matrix(NA,nrow = 512, ncol = 448)
#' D0_EVI <- matrix(NA,nrow =512, ncol = 448)
#' for(i in 1:448){
#'   D0_NDVI[,i] <- density(NDVI[,,i],from=0,to=1)$y
#'   D0_EVI[,i] <- density(EVI[,,i],from=0,to=1)$y
#' }
#' ## Define functional objects
#' d <- 11
#' basis <- create.bspline.basis(c(0,1),d)
#' u <- seq(0,1,length.out = 512)
#' y_NDVI <- smooth.basis(u,as.matrix(D0_NDVI),basis)$fd
#' y_EVI <- smooth.basis(u,as.matrix(D0_EVI),basis)$fd
#' y=list(y_NDVI,y_EVI)
#' ## Define functional time series
#' Y=fts(y)
#' plot(Y)
#' L=45
#' ## Decomposition stage of multivariate functional singular spectrum analysis
#' U=fssa(Y,L)
#' mfwcor=fwcor(U = U,groups = list(1,2,3,4))
#' wplot(W=mfwcor)
#' }
#'
#'
#' @seealso \code{\link{fssa}}, \code{\link{freconstruct}}, \code{\link{fts}}, \code{\link{wplot}}
#' @export
fwcor <- function(U, groups) {
  if(length(U$Y@C)==1) out <- ufwcor(U, groups) else out <- mfwcor(U, groups)
  return(out)
}
