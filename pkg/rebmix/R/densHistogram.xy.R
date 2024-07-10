.densHistogram.xy <- function(k, x, y, x0, xmin, xmax, y0, ymin, ymax, hx, hy, cx, cy, px, py)
{
  output <- .C(C_RdensHistogramXY,
    k = as.integer(k),
    n = as.integer(length(x)),
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
    error = integer(9),
    PACKAGE = "rebmix")

  error <- error.to.string(output$error);
      
  if (error[1] != "") {
    stop(error[1], call. = FALSE); return(NA)
  }
    
  if (error[2] != "") {
    warning(error[2], call. = FALSE, immediate. = TRUE)
  }  
    
  if (error[3] != "") {
    warning(error[3], call. = FALSE, immediate. = TRUE)
  }

  length(output$x) <- output$k
  length(output$y) <- output$k
  length(output$z) <- output$k

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densHistogram.xy
