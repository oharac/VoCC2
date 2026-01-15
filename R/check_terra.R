#' Internal. Angle associated to the spatial gradient
#' @param r A raster type object to be checked for compatibility with terra::SpatRaster
#' @author Casey O'Hara
#' check_terra(r)

check_terra <- function(r) {
  if(class(r) != 'SpatRaster') {
    pkg <- attr(class(r), 'package')
    if(pkg == 'raster') {
      warning('This function prefers SpatRaster objects; coercing raster to SpatRaster')
      r <- terra::rast(r)
    } else {
      stop('This function expect SpatRaster objects or raster objects, not objects from ', pkg, ' package')
    }
  }
  return(r)
}
