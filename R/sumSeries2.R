#' Summarize climatic series to higher temporal resolution
#'
#' Function to convert climatic series (provided as \code{RasterStack}) into a
#' coarser time frequency series for a period of interest. This function transforms the \code{RasterStack}
#' into an \code{xts} time series object to extract the values for the period of interest and
#' apply some summary function. It is mainly a wrapper from the \code{apply.} function family
#' in the package xts (Ryan and Ulrich 2017).
#'
#' @usage sumSeries(r, p, yr0, l, fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months",
#' freqout = "years")
#'
#' @param r \code{RasterStack} containing the time series of the climatic variable.
#' @param p \code{character string} defining the period to extract for the calculation
#' of the series (see examples).
#' @param yr0 \code{character string} specifying the first (yr0) year in the series (see examples).
#' @param l \code{integer} length of the input time series.
#' @param fun \code{logical} summary function to be computed. Summary functions need to be applied by cell (columns)
#' so should have the structure 'function(x) apply(x, 2, function(y){...})'. For convinience, sumSeries imports
#' colMaxs, and colMins from package ‘matrixStats’ (Bengtsson 2018) so they can be called in directly.
#' @param freqin \code{character string} specifying the original time frequency of the series.
#' @param freqout \code{character string} specifying the desired time frequency of the new series.
#' Must be one of the following: "weeks", "months", "quarters", "years", "other". Argument "other"
#' allows for user-defined functions to be applied on the 'xts' time series object over the period of interest (see examples).
#'
#' @return A \code{RasterStack} with the new series.
#'
#' @references \href{https://CRAN.R-project.org/package=xts}{Ray and Ulrich. 2017}. xts: eXtensible Time Series. R package version 0.10-1. \cr
#' \href{https://CRAN.R-project.org/package=matrixStats}{Bengtsson 2018}. matrixStats: Functions that Apply to Rows and Columns
#' of Matrices (and to Vectors). R package version 0.53.1.
#'
#' @importFrom matrixStats colMaxs colMins
#' @importFrom xts xts apply.weekly apply.monthly apply.quarterly apply.yearly
#' @export
#' @author Jorge Garcia Molinos
#' @examples
#' # Monthly mean SST (HadISST) data for Europe Jan-1950 to Dec-2010
#'
#' ?HSST
#'
#' # Calculate mean annual monthly SST
#'
#' yrSST <- sumSeries(HSST, p = "1969-01/2009-12", yr0 = "1955-01-01", l = nlayers(HSST),
#' fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years")
#'
#' # Extract Jul Aug mean SST each year (xts months are indexed from 0 to 11)
#'
#' myf = function(x, m = c(7,8)){
#' x[xts::.indexmon(x) %in% (m-1)]
#' }
#'
#' JlAugSST <- sumSeries(HSST, p = "1969-01/2009-12", yr0 = "1950-01-01", l = raster::nlayers(HSST),
#' fun = myf, freqin = "months", freqout = "other")
#'
#' # Same but calculating the annual variance of the two months
#'
#' myf = function(x, m = c(7,8)){
#' x1 <- x[xts::.indexmon(x) %in% (m-1)]
#' xts::apply.yearly(x1, function(y) apply(y, 2, function(y){var(y, na.rm = TRUE)}))
#' }
#'
#' meanJASST <- sumSeries(HSST, p = "1969-01/2009-12", yr0 = "1950-01-01", l = raster::nlayers(HSST),
#' fun = myf, freqin = "months", freqout = "other")
#'
#' @rdname sumSeries
#'

sumSeries2 <- function(r,           ### SpatRast object
                       p    = NULL, ### auto extract from yyyy-mm-dd format
                       yr0  = NULL, ### auto extract
                       l    = NULL, ### auto detect
                       fun = function(x) colMeans(x, na.rm = TRUE),
                       freqin = c('weeks', 'months', 'quarters', 'years', 'other')[2], ### default months
                       freqout = "years") {
  ### get layer names for possible date detection
  layer_names <- names(r)
  if(is.null(l)) l <- nlyr(r)
  if(is.null(p)) {
    message('No period provided; extracting from layer names expecting yyyy-mm-dd')
    name_dates <- stringr::str_extract(layer_names, '[0-9]{4}.[0-9]{2}.[0-9]{2}') %>%
      stringr::str_replace_all('[^0-9]', '-')
    if(length(name_dates) == 0) stop('No dates detected! Expecting format yyyy-mm-dd')
    if(length(name_dates) != l) stop('Extracted dates n = ', length(name_dates), ' but layers n = ', l)
    name_dates <- as.Date(name_dates)
    period <- paste(name_dates[1], name_dates[length(name_dates)], sep = '/')
  } else {
    period <- p
  }

  ### construct xts object
  r_mtx <- t(values(r))  ### pretty slow
  if(is.null(p)) {
    dates <- as.Date(name_dates)
  } else dates <- seq(as.Date(yr0), length = l, by = freqin)

  ts1 <- xts::xts(r_mtx, order.by = dates)

  ### subset for the period of interest using xts indexing
  ts1_sub <- ts1[period]

  ### calculate the annual series using apply.X functions from xts
  if(freqout == "weeks")    s <- apply.weekly(ts1_sub, fun)
  if(freqout == "months")   s <- apply.monthly(ts1_sub, fun)
  if(freqout == "quarters") s <- apply.quarterly(ts1_sub, fun)  ### Jan-Mar (Q1); Apr-Jn (Q2); Jl-Sep(Q3); Oct-Dec (Q4)
  if(freqout == "years")    s <- apply.yearly(ts1_sub, fun)
  if(freqout == "other")    s <- fun(ts1_sub)

  ### create raster stack
  r1 <- stack()
  for(i in 1:nrow(s)){
    r2 <- raster(r[[1]])
    r2[] <-  as.numeric(s[i,])
    r1 <- addLayer(r1, r2)
  }
  if(freqout != "other") {
    names(r1) <- seq(start(ts1_sub), length = nlayers(r1), by = freqout)
  }
  return(r1)
}

