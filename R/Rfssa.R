#' Rfssa: A Package for Functional Singular Spectrum Analysis and Related Methods.
#'
#' The Rfssa package provides the collection of necessary functions to
#' implement functional singular spectrum analysis (FSSA)-based methods for
#' analyzing univariate and multivariate functional time series (FTS).
#' Univariate and multivariate FSSA are novel, non-parametric methods used to perform decomposition and reconstruction of univariate and multivariate FTS
#' respectively. In addition, the FSSA-based routines may be performed on FTS
#' whose variables are observed over a one or two-dimensional domain. Finally,
#' one may perform FSSA recurrent or fssa vector forecasting of univariate or
#' multivariate FTS observed over one-dimensional domains. Forecasting of FTS
#' whose variables are observed over domains of dimension greater than one is
#' under development.
#'
#'
#' @details
#' The use of the package starts with the decomposition of functional time
#' series (\code{\link{fts}}) objects using the \code{\link{fssa}} routine. Then a suitable grouping of the
#' principal components is required for reconstruction (\code{\link{freconstruct}}) or
#' forecasting (\code{\link{fforecast}}) which can be done heuristically by looking at the
#' plots of the decomposition. Once a suitable grouping is chosen,
#' one may perform reconstruction where the sum of all the elements between the disjoint
#' groups approximates the original FTS. One may also choose to perform
#' forecasting after a grouping is chosen which returns future observations in
#' each FTS specified by the groups.
#'
#' This version of the package leverages a new S4 object for FTS objects (\code{\link{fts}}).
#' Along with providing the raw, sampled data, the new object may be specified using a provided basis and grid, a requested
#' basis and grid, or a mixture of provided and requested elements. We note that
#' the FTS object may be univariate or multivariate and variables may be observed
#' over one or two-dimensional domains. Validity checking of the S4 object
#' constructor inputs was also added to help guide the user. The plotting of FTS
#' objects was also updated to allow the user to plot FTS variables
#' observed over two-dimensional domains. Next, the FSSA routine (\code{\link{fssa}}) was
#' updated to perform faster by leveraging the RSpectra and RcppEigen R packages,
#' and the Eigen C++ package. We achieved a roughly 20 times speed up for
#' certain data examples. We updated the plotting of fssa objects to
#' allow for plotting of left singular functions that correspond with FTS
#' variables observed over a two-dimensional domain.
#' We updated FSSA reconstruction \code{\link{freconstruct}} to handle
#' FTS whose variables are observed over one or two-dimensional domains. We also
#' updated FTS arithmetic (such as FTS addition, FTS subtraction, etc.) to allow
#' the user to perform scalar-FTS arithmetic on different variables of a
#' multivariate FTS.
#'
#' The first piece of new functionality that has been added is that the user
#' may now specify univariate or multivariate FTS comprised of variables observed
#' over one or two-dimensional domains. In addition, forecasting of univariate
#' and multivariate FTS observed over one-dimensional domains by FSSA/MFSSA
#' recurrent forecasting and FSSA/MFSSA vector forecasting has also been added.
#' We have also added in a new data set (\code{\link{Montana}}) which provides the data for a
#' multivariate FTS observed over different dimensional domains.
#'
#' The package update also includes updates to the shiny app (\code{\link{launchApp}}) that can be used for demonstrations of univariate or multivariate FSSA
#' depending on the type that is specified.
#' The app allows the user to explore FSSA with simulated data, data that is provided on the server, or data that the user provides.
#' It allows the user to change parameters as they please, gives visual results of the methods, and also allows the user to compare FSSA results to other
#' spectrum analysis methods such as multivariate singular spectrum analysis. The tool is easy to use and can act as a nice starting point for a user that wishes to
#' perform FSSA as a part of their data analysis.
#'
#'
#' @seealso
#'  \code{\link{fssa}}, \code{\link{freconstruct}}, \code{\link{fforecast}}
#'  \code{\link{fwcor}}, \code{\link{wplot}}, \code{\link{fts}}, \code{\link{plot.fts}}, \code{\link{plot.fssa}},
#'  \code{\link{cor.fts}}, \code{\link{launchApp}}
#'
#'
#'
#' @references
#'   Haghbin, H., Morteza Najibi, S., Mahmoudvand, R., Trinka, J., and Maadooliat, M. (2021). Functional
#'   singular spectrum analysis. Stat. e330 STAT-20-0240.R1.
#'
#'   Trinka J., Haghbin H., Maadooliat M. (Accepted) Multivariate Functional Singular Spectrum Analysis: A Nonparametric Approach for Analyzing Multivariate Functional Time Series. In: Bekker A., Ferreira, J., Arashi M., Chen D. (eds) Innovations in Multivariate Statistical Modeling: Navigating Theoretical and Multidisciplinary Domains. Emerging Topics in Statistics and Biostatistics. Springer, Cham.
#'
#'   Trinka J. (2021) Functional Singular Spectrum Analysis: Nonparametric Decomposition and Forecasting Approaches for Functional Time Series [Doctoral dissertation, Marquette University]. ProQuest Dissertations Publishing.
#'
#'   Trinka, J., Haghbin, H., and Maadooliat, M. (2021). Functional time series forecasting: Functional
#'   singular spectrum analysis approaches. Version 4 retrieved from https://arxiv.org/abs/2011.
#'   13077.
#'
#'
#' @docType package
#' @name Rfssa
NULL
