#' Gradient-based climate velocity
#'
#' Function to calculate the velocity of climate change after Burrows et al. (2011)
#' based on local climatic temporal trends and spatial gradients.
#'
#' @usage gVoCC(clim_trend, clim_grad)
#'
#' @param clim_trend A \code{SpatRaster} output from the \code{tempTrend}
#'   function containing the long-term linear climatic trends, in units of
#'   climate variable (e.g., temp in deg Celsius) per time variable (e.g.,
#'   year)
#' @param clim_grad A \code{SpatRaster} output from the \code{spatGrad}
#'   function containing the magnitudes and angles for the spatial climatic
#'   gradient.
#'
#' @return A \code{SpatRaster} containing the climate velocity magnitude ("vocc_mag",
#'   km/temporal unit for unprojected rasters and spatial unit/temporal unit for projected rasters)
#'   and angle("vocc_ang" in degrees north: 0N, 90E, 180S and 270W).
#'
#' @references \href{http://science.sciencemag.org/content/334/6056/652}{Burrows et al. 2011}. The pace of shifting climate
#' in marine and terrestrial ecosystems. Science, 334, 652-655.
#'
#' @seealso{\code{\link{tempTrend}}, \code{\link{spatGrad}}}
#'
#' @export
#' @author Jorge Garcia Molinos
#' @examples
#'
#' ?HSST
#' yrSST <- sumSeries(HSST, p = "1960-01/2009-12", yr0 = "1955-01-01", l = nlayers(HSST),
#' fun = function(x) colMeans(x, na.rm = TRUE),
#' freqin = "months", freqout = "years")
#' tr <- tempTrend(yrSST, th = 10)
#' sg <- spatGrad(yrSST, th = 0.0001, projected = FALSE)
#'
#' # Magnitude and angle of the climate velocity (km/yr) 1960-2009
#'
#' v <- gVoCC(tr,sg)
#' plot(v)
#'
#' @rdname gVoCC

gVoCC2 <- function(clim_trend, clim_grad) {
  ### check for SpatRaster
  check_terra(clim_trend); check_terra(clim_grad)

  vocc <- clim_trend[[1]]/clim_grad[[1]]
  # velocity angles have opposite direction to the spatial climatic gradient if warming and same direction (cold to warm) if cooling
  ind <- which(values(vocc) > 0)
  vocc_ang <- clim_grad[[2]]
  vocc_ang[ind] <- clim_grad[[2]][ind] + 180
  vocc_ang[] <- ifelse(vocc_ang[] >= 360, vocc_ang[] - 360, vocc_ang[])
  output <- c(vocc, vocc_ang)
  names(output) <- c("vocc_mag", "vocc_ang")
  return(output)
}


