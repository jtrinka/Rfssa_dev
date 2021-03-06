#' Functional Singular Spectrum Analysis
#'
#' This is a function which performs the decomposition (including embedding
#'  and  functional SVD steps) stage for univariate functional singular spectrum analysis (ufssa)
#'  or multivariate functional singular spectrum analysis (mfssa)
#'  depending on whether the supplied input is a univariate or
#'  multivariate functional time series (\code{\link{fts}}) object. Also note that the variables of the \code{\link{fts}} maybe observed over different dimensional domains where the maximum dimension currently supported is two.
#' @return An object of class \code{fssa}, which is a list of functional objects and the following components:
#' \item{values}{A numeric vector of eigenvalues.}
#' \item{L}{The specified window length.}
#' \item{N}{The length of the functional time series.}
#' \item{Y}{The original functional time series.}
#' @param Y An object of class \code{\link{fts}}.
#' @param L A positive integer giving the window length.
#' @param ntriples A positive integer specifying the number of eigentriples to calculate in the decomposition.
#' @param type A string indicating which type of fssa to perform. Use \code{type="ufssa"} to perform univariate fssa (default for univariate fts). Use \code{type="mfssa"} to perform multivariate fssa (default for multivariate fts).
#' @importFrom fda fd inprod eval.fd smooth.basis is.fd create.bspline.basis
#' @importFrom RSpectra eigs
#' @examples
#' \dontrun{
#' ## Univariate FSSA Example on Callcenter data
#' require(Rfssa)
#' load_github_data("https://github.com/jtrinka/Rfssa_dev/blob/main/data/Callcenter.RData")
#' ## Define functional objects
#' D <- matrix(sqrt(Callcenter$calls), nrow = 240)
#' N <- ncol(D)
#' time <- seq(ISOdate(1999, 1, 1), ISOdate(1999, 12, 31), by = "day")
#' K <- nrow(D)
#' u <- seq(0, K, length.out = K)
#' d <- 22 # Optimal Number of basis elements
#' ## Define functional time series
#' Y <- fts(list(D), list(list(d, "bspline")), list(u))
#' Y
#' plot(Y, mains = c("Sqrt of Call Center Data"))
#' ## Univariate functional singular spectrum analysis
#' L <- 28
#' U <- fssa(Y, L)
#' plot(U, d = 13)
#' plot(U, d = 9, type = "lheats")
#' plot(U, d = 9, type = "lcurves")
#' plot(U, d = 9, type = "vectors")
#' plot(U, d = 10, type = "periodogram")
#' plot(U, d = 10, type = "paired")
#' plot(U, d = 10, type = "wcor")
#' gr <- list(1, 2:3, 4:5, 6:7, 8:20)
#' Q <- freconstruct(U, gr)
#' plot(Y, mains = "Sqrt of Call Center Data")
#' plot(Q[[1]], mains = "1st Component")
#' plot(Q[[2]], mains = "2nd Component")
#' plot(Q[[3]], mains = "3rd Component")
#' plot(Q[[4]], mains = "4th Component")
#' plot(Q[[5]], mains = "5th Component (Noise)")
#'
#' ## Other visualisation types for object of class "fts":
#'
#' plot(Q[[1]], type = "3Dsurface", xlabels = "Intraday intervals", tlabels = "Day", zlabels = "Output")
#' # Visualizing the first 60 observations in the reconstructed fts.
#' plot(Q[[2]][1:60], type = "heatmap", xlabels = "Intraday intervals")
#' plot(Q[[3]][1:60], type = "3Dline", xlabels = "Intraday intervals", tlabels = "Day", zlabels = "Output")
#'
#' ## Multivariate FSSA Example on bivariate intraday
#' ## temperature curves and smoothed images of vegetation
#' require(Rfssa)
#' load_github_data("https://github.com/jtrinka/Rfssa_dev/blob/main/data/Montana.RData")
#' Temp <- Montana$Temp
#' NDVI <- Montana$NDVI
#' d_temp <- 11
#' d_NDVI <- 13
#' ## Define functional time series
#' Y <- fts(
#'   list(Temp / sd(Temp), NDVI), list(
#'     list(d_temp, "bspline"),
#'     list(d_NDVI, d_NDVI, "bspline", "bspline")
#'   ),
#'   list(c(0, 23), list(c(1, 33), c(1, 33)))
#' )
#' # Plot the first 100 observations
#' plot(Y[1:100],
#'   xlabels = c("Time", "Lon."), ylabels = c("Temperature (\u00B0C)", "Lat."),
#'   zlabels = c("", "NDVI"), mains = c("Temperature Curves", "NDVI Images")
#' )
#' plot(Y, types = c("3Dline", "heatmap"), vars = c(1, 1))
#' plot(Y, types = "heatmap", vars = 2)
#' plot(Y, vars = c(2, 1))
#' L <- 45
#' ## Multivariate functional singular spectrum analysis
#' U <- fssa(Y, L)
#' plot(U, type = "values", d = 10)
#' plot(U, type = "vectors", d = 4)
#' plot(U, type = "lheats", d = 4)
#' plot(U, type = "lcurves", d = 4, vars = c(1))
#' plot(U, type = "paired", d = 6)
#' plot(U, type = "wcor", d = 10)
#' plot(U, type = "periodogram", d = 4)
#' # Reconstruction of multivariate fts observed over different dimensional domains
#' Q <- freconstruct(U = U, groups = list(c(1), c(2:3), c(4)))
#' # Plotting reconstructions to show accuracy
#' plot(Q[[1]]) # mean
#' plot(Q[[2]]) # periodic
#' plot(Q[[3]]) # trend
#' }
#' @useDynLib Rfssa
#' @export
fssa <- function(Y, L = NA, ntriples = 20, type = "ufssa") {
  N <- ncol(Y@C[[1]])
  if (is.na(L)) L <- floor(N / 2L)
  if (length(Y@C) == 1 && type == "ufssa") {
    out <- ufssa(Y, L, ntriples)
  } else if (length(Y@C) > 1 || type == "mfssa") {
    out <- mfssa(Y, L, ntriples)
  } else {
    stop("Error in type or dimension.")
  }
  class(out) <- "fssa"
  return(out)
}
