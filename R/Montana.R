#' Montana Intraday Temperature Curves and NDVI Images Data Set
#'
#' This data set contains the intraday hourly temperature curves measured in degrees celcius and normalized difference vegetation index (NDVI) image data where both types of data are recorded near Saint Mary, Montana, USA. The NDVI images are taken of a region located between longitudes of 113.30 degrees West and 113.56 degrees west and latitudes of 48.71 degrees North and 48.78 degrees North. The intraday temperature curves are available for download from Diamond et al. (2013) and the NDVI images were attained leveraging resources provided by Tuck et al. (2014).
#' For each recorded intraday temperature curve, an NDVI image was recorded on the same day, every 16 days, starting January 1, 2008 and ending September 30, 2013.
#' With the the threat of global warming damaging various ecosystems, the goal of the study was to analyze trends in the temperature and to investigate how changes in temperature effects the amount of vegetation in the region. We discovered that leveraging both types of variables in a multivariate analysis revealed a stronger signal extraction result and more informative patterns. The data is hosted on GitHub and \code{\link{load_github_data}} may be used to load the data.
#' @name Montana
#' @format A list which contains a 24 by 133 matrix of discrete samplings of intraday hourly temperature curves and an array that is 33 by 33 by 133 where one 33 by 33 slice of the array is an NDVI image.
#' @references
#' \enumerate{
#' \item Diamond, H. J., Karl, T., Palecki, M. A., Baker, C. B., Bell, J. E., Leeper, R. D.,
#' Easterling, D. R., Lawrimore, J. H., Meyers, T. P., Helfert, M. R., Goodge, G.,
#' and Thorne, P.W. (2013). U.S. climate reference network after one decade of operations:
#' status and assessment. https://www.ncdc.noaa.gov/crn/qcdatasets.html.
#' Last accessed April 2020.
#'
#' \item Tuck, S. L., Phillips, H. R., Hintzen, R. E., Scharlemann, J. P., Purvis, A., and
#' Hudson, L. N. (2014). MODISTools – downloading and processing MODIS
#' remotely sensed data in R. Ecology and Evolution, 4(24):4658–4668.
#'
#'
#'
#' }
#' @seealso \code{\link{fssa}}
NULL
