.densKNearestNeighbour.xy <- function(x, y, k, hx, hy)
{
  output <- .C("RdensKNearestNeighbourXY",
    n = as.integer(length(x)),
    x = as.double(x),
    y = as.double(y),
    z = double(length(x)),
    k = as.integer(k),
    hx = as.double(hx),
    hy = as.double(hy),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densKNearestNeighbour.xy!", call. = FALSE); return(NA)
  }

  i <- order(output$z)

  output$x <- output$x[i]
  output$y <- output$y[i]
  output$z <- output$z[i]

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densKNearestNeighbour.xy
