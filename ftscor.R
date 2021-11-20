#' Correlation for Functional Time Series Objects
#'
#' This function finds the correlation between univarite or multivariate functional time series (\code{\link{fts}}) objects.
#' @return a scalar that is the correlation between \code{\link{fts}} objects
#' @param Y1 an object of class \code{\link{fts}}
#' @param Y2 an object of class \code{\link{fts}}
#'
#'
#' @seealso \code{\link{fts}}
#'
#' @examples
#'
#' \dontrun{
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
#' Y_NDVI <- fts(y_NDVI)
#' Y_EVI <- fts(y_EVI)
#' cor.fts(Y_NDVI,Y_EVI)
#' }
#'
#' @export
cor.fts <- function(Y1, Y2) {
  if(ncol(Y1@C[[1]])!=ncol(Y2@C[[1]])){

    stop("Functional time series are of different length")

  }
  if(length(Y1)!=length(Y2)){

    stop("Functional time series have different number of covariates")

  }
  N <- ncol(Y1@C[[1]])
  w <- matrix(nrow=1,ncol=N,data=1)
  p = length(Y1)
  G = list()
  Y1_list=list()
  Y2_list=list()
  for(j in 1:p){
    G[[j]] = inprod(Y1@B[[j]],Y1@B[[j]])
    Y1_list[[j]]=Y1@C[[j]]
    Y2_list[[j]]=Y2@C[[j]]
  }

  wcor <- mwinprod(Y1_list, Y2_list, w, G, p)/sqrt(mwinprod(Y1_list,Y1_list, w, G, p) * mwinprod(Y2_list,Y2_list, w, G, p))

  return(wcor)
}
