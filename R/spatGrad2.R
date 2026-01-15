#' Local spatial climatic gradients
#'
#' Function to calculate the magnitude and direction of the spatial gradient
#' associated to a climatic variable after Burrows et al. (2011). This trend is
#' to be used for the calculation of the gradient-based climate velocity using gVoCC.
#'
#' @usage spatGrad(r, th = -Inf, projected = FALSE)
#'
#' @param r \code{SpatRaster} with the annual climatic values for the period of interest.
#'   Alternatively, a \code{raster} with the annual climatic values averaged
#'   over the period of interest.
#' @param th \code{Integer} indicating a lower threshold to truncate the spatial
#'   gradient with. Use -Inf (default) if no threshold required.
#' @param projected \code{Logical} is the source raster in a projected coordinate system?
#'   If FALSE (default) a correction will be made to account for latitudinal distortion.
#'
#' @return A \code{SpatRaster} with the magnitude of the spatial gradient
#'   (grad in C per km for unprojected rasters and C per spatial unit for projected rasters),
#'   and the associated angle (ang in degrees).
#'
#' @references \href{http://science.sciencemag.org/content/334/6056/652}{Burrows et al. 2011}. The pace of shifting climate in marine and terrestrial ecosystems. Science, 334, 652-655.
#'
#' @seealso{\code{\link{tempTrend}}, \code{\link{gVoCC}}}
#'
#' @import data.table
#' @import terra
#' @importFrom CircStats rad deg
#' @export
#' @author Jorge Garcia Molinos, David S. Schoeman, and Michael T. Burrows
#' @examples
#'
#' yrSST <- sumSeries(HSST, p = "1969-01/2009-12", yr0 = "1955-01-01", l = nlayers(HSST),
#' fun = function(x) colMeans(x, na.rm = TRUE), freqin = "months", freqout = "years")
#'
#' # Spatial gradient (magnitude and angle) for the average mean annual SST.
#'
#' sg <- spatGrad(yrSST, th = 0.0001, projected = FALSE)
#'
#' plot(sg)
#'
#' @rdname spatGrad

spatGrad2 <- function(clim_r, th = -Inf, projected = FALSE) {

  if(class(clim_r) != 'SpatRaster') stop('This function expects spatial data in SpatRaster format')
  if(nlyr(clim_r) > 1) {
    r <- mean(clim_r, na.rm=T)
  } else {
    r <- clim_r
  }
  ### get resolution of the raster
  re <- res(r)

  ### Create a columns for focal and each of its 8 adjacent cells
  y <- data.table(adjacent(r, 1:ncell(r), directions = 8, pairs = TRUE))

  ### Get value for focal cell, order the table by raster sequence and omit NAs (land cells)
  y <- y[ , clim_focal := values(r)[from]][order(from, to)] |>
    na.omit()

  ### Insert values for adjacent cells
  y[ , clim := values(r)[to]]

  ### UPDATE:  if there are any isolated pixels present (i.e., pixel with no
  ### adjacent pixels), the cell will be dropped from the angle raster, which
  ### may cause problems when comparing to the magnitude raster...
  ### Put a check in for that!  How to deal - drop them (set to NA) or
  ### throw an error?
  y_check <- y[, .(adj_count = sum(!is.na(clim))), by = from]
  if(any(y_check$adj_count == 0)) {
    orphan_cells <- y_check$from[y_check$adj_count == 0]
    # stop('Isolated non-NA pixel(s) with no adjacent non-NA cells: ',
    #      paste0(orphan_cells, collapse = ', '))
    warning('Isolated non-NA pixel(s) with no adjacent non-NA cells: ',
            paste0(orphan_cells, collapse = ', '),
            '... setting values to NA')
    ### set local climate raster to NA
    values(r)[orphan_cells] <- NA
    ### drop rows from data table
    y <- y[!from %in% orphan_cells]
  }


  ### Column to identify rows in the raster (N = 1, mid = 0, S = -1)
  y[ , sy := rowFromCell(r, from) - rowFromCell(r, to)]
  ### Same for columns (E = 1, mid = 0, W = -1)
  y[ , sx := colFromCell(r, to)   - colFromCell(r, from)]

  ### Sort out the W-E wrap at the dateline, part I and part II
  y[sx > 1, sx := -1]
  y[sx < -1, sx := 1]

  ### Make a unique code for each of the eight neighbouring cells
  y[ , code := paste0(sx, sy)]

  ### Code cells with positions
  y[.(code = c("10","-10","-11","-1-1","11","1-1","01","0-1"),
      to = c("clim_e","clim_w","clim_nw","clim_sw","clim_ne","clim_se","clim_n","clim_s")),
    on = "code", code := i.to]
  y <- dcast(y[ , c("from", "code", "clim")], from ~ code, value.var = "clim")
  y[ , clim_focal := values(r)[from]]   ### Put clim_focal back in
  y[ , lat := yFromCell(r, from)]      ### Add focal cell latitude

  ### Calculate individual spatial temperature gradients: grads (degC per km)
  ### WE gradients difference in temperatures for each western and eastern
  ### pairs divided by the distance between the cells in each pair (corrected
  ### for latitudinal distortion if unprojected)
  ### Positive values indicate an increase in clim from W to E (i.e., in line with the Cartesian x axis)

  ifelse(projected == TRUE, d <- 1, d <- 111.325)
  ifelse(projected == TRUE, co <- 0, co <- 1)

  y[ , grad_we1 := (clim_n - clim_nw) / (cos(co * rad(lat + re[2])) * (d * re[1]))]
  y[ , grad_we2 := (clim_focal - clim_w) / (cos(co * rad(lat)) * (d * re[1]))]
  y[ , grad_we3 := (clim_s - clim_sw) / (cos(co * rad(lat - re[2])) * (d * re[1]))]
  y[ , grad_we4 := (clim_ne - clim_n) / (cos(co * rad(lat + re[2])) * (d * re[1]))]
  y[ , grad_we5 := (clim_e - clim_focal) / (cos(co * rad(lat)) * (d * re[1]))]
  y[ , grad_we6 := (clim_se - clim_s) / (cos(co * rad(lat - re[2])) * (d * re[1]))]

  ### NS gradients difference in temperatures for each northern and southern pairs divided by the distance between them (111.325 km per degC *re[2] degC)
  ### Positive values indicate an increase in sst from S to N (i.e., in line with the Cartesian y axis)
  y[ , grad_ns1 := (clim_nw - clim_w) / (d * re[2])]
  y[ , grad_ns2 := (clim_n - clim_focal) / (d * re[2])]
  y[ , grad_ns3 := (clim_ne - clim_e) / (d * re[2])]
  y[ , grad_ns4 := (clim_w - clim_sw) / (d * re[2])]
  y[ , grad_ns5 := (clim_focal - clim_s) / (d * re[2])]
  y[ , grad_ns6 := (clim_e - clim_se) / (d * re[2])]

  ### Calulate NS and WE gradients. NOTE: for angles to work (at least using
  ### simple positive and negative values on Cartesian axes), S-N & W-E
  ### gradients need to be positive)
  y[ , we_grad := apply(.SD, 1,
                        function(x) stats::weighted.mean(x, c(1, 2, 1, 1, 2, 1), na.rm = T)),
     .SDcols = 12:17]
  y[ , ns_grad := apply(.SD, 1,
                        function(x) stats::weighted.mean(x, c(1, 2, 1, 1, 2, 1), na.rm = T)),
     .SDcols = 18:23]
  y[is.na(we_grad) & !is.na(ns_grad), we_grad := 0L] ### Where ns_grad does not exist, but we_grad does, make ns_grad 0
  y[!is.na(we_grad) & is.na(ns_grad), ns_grad := 0L] ### same the other way around

  ### Calculate angles of gradients (degrees) - adjusted for quadrant (0 deg is North)
  y[ , angle := angulo(.SD$we_grad, .SD$ns_grad), .SDcols = c("we_grad", "ns_grad")]

  ### Calculate the vector sum of gradients (C/km)
  y[ , grad := sqrt(apply(cbind((y$we_grad^2), (y$ns_grad^2)), 1, sum, na.rm = TRUE))]

  ### Merge the reduced file back into the main file to undo the initial na.omit
  from <- data.table(1:ncell(r)) ### Make ordered from cells
  y <- y[from]   ### merge both

  ang_r <- grad_r <- rast(r)
  ang_r[y$from] <- y$angle
  grad_r[y$from] <- y$grad
  grad_r[grad_r[] < th] <- th
  output <- c(grad_r, ang_r)
  names(output) <- c("grad", "ang")
  return(output)

}
