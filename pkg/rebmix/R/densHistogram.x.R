.densHistogram.x <- function(k, x, x0, xmin, xmax, hx, cx, px)
{
  output <- .C(C_RdensHistogramX,
    k = as.integer(k),
    n = as.double(length(x)),
    x = as.double(x),
    y = double(length(x)),
    x0 = as.double(x0),
    xmin = as.double(xmin),    
    xmax = as.double(xmax),
    hx = as.double(hx),
    px = as.character(px),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RdensHistogramX!", call. = FALSE); return(NA)
  }

  length(output$x) <- output$k
  length(output$y) <- output$k

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densHistogram.x
