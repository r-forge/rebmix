.densHistogram.xy <- function(k, x, y, x0, xmin, xmax, y0, ymin, ymax, hx, hy, cx, cy, px, py)
{
  output <- .C(C_RdensHistogramXY,
    k = as.integer(k),
    n = as.double(length(x)),
    x = as.double(x),
    y = as.double(y),
    z = double(length(x)),
    x0 = as.double(x0),
    xmin = as.double(xmin),    
    xmax = as.double(xmax),    
    y0 = as.double(y0),
    ymin = as.double(ymin),    
    ymax = as.double(ymax),    
    hx = as.double(hx),
    hy = as.double(hy),
    px = as.character(px),
    py = as.character(py),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RdensHistogramXY!", call. = FALSE); return(NA)
  }

  length(output$x) <- output$k
  length(output$y) <- output$k
  length(output$z) <- output$k

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densHistogram.xy
