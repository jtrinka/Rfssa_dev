#' Weighted-Correlations Plot
#'
#' This function displays a plot of the weighted-correlation (w-correlation) matrix of functional time series (\code{\link{fts}}) objects that were reconstructed
#' from functional singular spectrum analysis (\code{\link{fssa}}) objects.
#' @param W a w-correlation matrix
#' @param cuts an integer that is the number of levels the range of w-correlation values will be divided into
#' @note refer to \code{\link{fwcor}} for an example on how to run this function starting from a w-correlation matrix
#' @seealso \code{\link{fssa}}, \code{\link{freconstruct}}, \code{\link{fts}}, \code{\link{fwcor}}

#' @export
wplot <- function(W, cuts=20) {
  at <- pretty(c(0,1), n = cuts)
  d <- nrow(W)
  W0 <- abs(W)
  a <- min(W0)
  b <- max(W0 - diag(1, d))
  s <- stats::sd(W0 - diag(1, d))
  diag(W0) <- min(1, b + 3 * s)
  xylabels <- paste0("F", 1:d)
  p1 <- lattice::levelplot(1 - W0, xlab = "", cex.main=2.5,at=at,
                           ylab = "", colorkey = NULL,
                           main = list(paste("(B) W-Correlation Matrix"),cex=2),
                           scales = list(x = list(at = 1:d,
                                                  lab = xylabels,cex=1.6),
                                         y = list(at = 1:d,
                                                  lab = xylabels,cex=1.6)),
                           col.regions = grDevices::gray(seq(0,                                                                                                    1, length = 100)))
  graphics::plot(p1)
}
