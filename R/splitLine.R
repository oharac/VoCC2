#' Internal. Split a line segment defined by points A-B into n parts
#' @param A \code{numeric} 2-d matrix or data.frame giving coordinates of first point
#' @param B \code{numeric} 2-d matrix or data.frame giving coordinates of second point
#' @param n \code{numeric} vector: number of segments to divide the distance between points with
#' @author Jorge Garcia Molinos
#' @rdname splitLine

splitLine <- function(A, B, n) {
  if(any(is.na(n))) stop('splitLine: NAs detected in `n`, not allowed!')
  if(nrow(A) != nrow(B) | nrow(A) != length(n)) stop('splitLine: Dimension mismatch!')
  # Remove sign change for dateline, if needed
  B[ , 1] <- with(B, ifelse(abs(B[ , 1] - A[ , 1]) > 180 &  B[ , 1] < 0,
                            B[ , 1] + 360,
                            B[ , 1]))
  A[ , 1] <- with(A, ifelse(abs(B[ , 1] - A[ , 1]) > 180 &  A[ , 1] < 0,
                            A[ , 1] + 360,
                            A[ , 1]))

  xstp <- (B[ , 1] - A[ , 1]) / n
  ystp <- (B[ , 2] - A[ , 2]) / n

  if(length(xstp) > 1) {
    points <- mapply(function(X, Y, xstp, ystp, n) {
      x <- as.numeric(X) + (as.numeric(xstp) * (0:n))
      # Adjust for dateline jumps / return to -180o to 180o format
      x <- x - (360 * floor((x + 180) / 360))
      y <- as.numeric(Y) + (as.numeric(ystp) * (0:n))
      do.call(cbind, list(x, y))
    }, A[ , 1], A[ , 2], xstp, ystp, n, SIMPLIFY = FALSE)
  }
  # if there is only 1 trajectory active no need for mapply
  if(length(xstp) == 1){
    x <- A[ , 1] + (xstp * (0:n))
    # Adjust for dateline jumps / return to -180o to 180o format
    x <- x - (360 * floor((x + 180) / 360))
    y <- A[ , 2] + (ystp * (0:n))
    points <- cbind(x, y)
  }
  return(points)
}


