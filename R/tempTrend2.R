#' Long-term local climatic trends
#'
#' Function to calculate temporal trend from a raster series
#' of a climatic variable. This trend is to be used for the calculation of the
#' gradient-based climate velocity using gVoCC.
#'
#' @usage tempTrend(clim_r, th)
#'
#' @param clim_r \code{SpatRaster} containing a time series of values of
#'   the climatic variable for the period of interest.  Note that the calculated
#'   trend will be in units of (climate variable units) per (time interval of
#'   input raster layers).  For example, if the input raster \code{r} contains
#'   annual values of temperature in degrees Celsius, the resulting trend will
#'   be in degrees Celsius per year.
#' @param th \code{Integer} minimum number of observations in the series needed to
#'   calculate the trend at each cell; default 10.
#'
#' @return A \code{SpatRaster} containing the cell-specific temporal trends
#'   extracted from simple linear regressions of the climatic variable against time
#'   ("slp_trends" in climate variable units per time interval of input raster
#'   \code{r}), together with their standard
#'   errors ("se_trends") and statistical significance ("pval_trends").
#'
#' @seealso{\code{\link{spatGrad}}, \code{\link{gVoCC}}}
#'
#' @import terra
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

tempTrend2 <- function(clim_r, th = 10) {

  if(class(clim_r) != 'SpatRaster') stop('This function expects spatial data in terra::SpatRaster format')

  ### Extract values as matrix
  y <- terra::values(clim_r)

  ### Identify valid cells (cells with at least one non-NA value)
  valid_cells <- which(rowSums(is.na(y)) != nlyr(clim_r))
  y <- t(y[valid_cells, ])

  ### Count non-NA observations per cell
  non_na <- apply(y, 2, function(x) sum(!is.na(x)))

  ### Keep only cells with at least th observations
  ind <- which(non_na >= th)
  y <- y[ , ind]
  n_obs <- apply(y, 2, function(x) sum(!is.na(x)))

  ### Create time index matrix
  x <- matrix(nrow = nlyr(clim_r), ncol = ncol(y))
  x[] <- 1:nlyr(clim_r)

  ### Put NA values into x to correspond with y
  x1 <- y
  x1[!is.na(x1)] <- 1
  x <- x * x1

  ### Calculate sum terms for linear regression
  sumx   <- apply(x, 2, sum, na.rm = TRUE)
  sumy   <- apply(y, 2, sum, na.rm = TRUE)
  sumxx  <- apply(x, 2, function(x) sum(x^2, na.rm = TRUE))
  sumyy  <- apply(y, 2, function(x) sum(x^2, na.rm = TRUE))
  prodxy <- x * y
  sumxy  <- apply(prodxy, 2, sum, na.rm = TRUE)

  ### Estimate slope coefficients and associated statistics
  slope <- (n_obs * sumxy - (sumx * sumy)) / (n_obs * sumxx - sumx^2)
  sres  <- (n_obs * sumyy - sumy^2 - slope^2 * (n_obs * sumxx - sumx^2)) / (n_obs * (n_obs - 2))
  se    <- suppressWarnings(sqrt((n_obs * sres) / (n_obs * sumxx - sumx^2)))
  test  <- slope / se
  pval  <- mapply(function(x, y) (2 * pt(abs(x), df = y - 2, lower.tail = FALSE)),
                  x = test, y = n_obs)

  ### Create output rasters
  slp_trends <- sig_trends <- se_trends <- rast(clim_r[[1]])
  slp_trends[valid_cells[ind]] <- slope
  se_trends[valid_cells[ind]]  <- se
  sig_trends[valid_cells[ind]] <- pval

  # Combine into SpatRaster stack
  output <- c(slp_trends, se_trends, sig_trends)
  names(output) <- c("slp_trends", "se_trends", "pval_trends")

  return(output)
}
