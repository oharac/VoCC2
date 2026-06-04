#' Long-term local climatic trends
#'
#' Function to calculate temporal trend from a raster series
#' of a climatic variable. This trend is to be used for the calculation of the
#' gradient-based climate velocity using gVoCC.
#'
#' @usage tempTrend(r, th)
#'
#' @param r \code{SpatRaster} containing a time series of (annual, seasonal, monthly...) values of
#' the climatic variable for the period of interest.
#' @param th \code{Integer} minimum number of observations in the series needed to
#' calculate the trend at each cell.
#'
#' @return A \code{SpatRaster} containing the cell-specific temporal trends
#' extracted from simple linear regressions of the climatic variable against time
#' ("slpTrends" in degree Celsius per year), together with their standard
#' errors ("seTrends") and statistical significance ("sigTrends").
#'
#' @seealso{\code{\link{spatGrad}}, \code{\link{gVoCC}}}
#'
#' @export
#' @author Jorge Garcia Molinos and Christopher J. Brown
#' @examples
#'
#' yrSST <- sumSeries(HSST, p = "1969-01/2009-12", yr0 = "1955-01-01", l = nlayers(HSST),
#' fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years")
#'
#' # Mean annual SST trend (minimum threshold of 10 years of data), with SE and p-values.
#'
#' tr <- tempTrend(yrSST, th = 10)
#'
#' plot(tr)
#'
#' @rdname tempTrend

tempTrend <- function(r, th) {
  # Extract values as matrix
  y <- terra::values(r)

  # Identify ocean cells (cells with at least one non-NA value)
  ocean <- which(rowSums(is.na(y)) != nlyr(r))
  y <- t(y[ocean, ])

  # Count non-NA observations per cell
  N <- apply(y, 2, function(x) sum(!is.na(x)))

  # Keep only cells with at least th observations
  ind <- which(N >= th)
  y <- y[ , ind]
  N <- apply(y, 2, function(x) sum(!is.na(x)))

  # Create time index matrix
  x <- matrix(nrow = nlyr(r), ncol = ncol(y))
  x[] <- 1:nlyr(r)

  # Put NA values into x to correspond with y
  x1 <- y
  x1[!is.na(x1)] <- 1
  x <- x * x1

  # Calculate sum terms for linear regression
  sumx   <- apply(x, 2, sum, na.rm = TRUE)
  sumy   <- apply(y, 2, sum, na.rm = TRUE)
  sumxx  <- apply(x, 2, function(x) sum(x^2, na.rm = TRUE))
  sumyy  <- apply(y, 2, function(x) sum(x^2, na.rm = TRUE))
  prodxy <- x * y
  sumxy  <- apply(xy, 2, sum, na.rm = TRUE)

  # Estimate slope coefficients and associated statistics
  slope <- (N * sumxy - (sumx * sumy)) / (N * sumxx - sumx^2)
  sres  <- (N * sumyy - sumy^2 - slope^2 * (N * sumxx - sumx^2)) / (N * (N - 2))
  se    <- suppressWarnings(sqrt((N * sres) / (N * sumxx - sumx^2)))
  test  <- slope / se
  p     <- mapply(function(x, y) (2 * pt(abs(x), df = y - 2, lower.tail = FALSE)),
                  x = test, y = N)

  # Create output rasters
  slp_trends <- sig_trends <- se_trends <- rast(r[[1]])
  slp_trends[ocean[ind]]   <- slope
  se_trends[ocean[ind]]    <- se
  sig_trends[ocean[ind]]   <- p

  # Combine into SpatRaster stack
  output <- c(slp_trends, se_trends, sig_trends)
  names(output) <- c("slpTrends", "seTrends", "sigTrends")

  return(output)
}
