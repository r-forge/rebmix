.densHistogram.x <- function(k, x, x0, xmin, xmax, hx, cx, px)
{
  output <- .C(C_RdensHistogramX,
    k = as.integer(k),
    n = as.integer(length(x)),
    x = as.double(x),
    y = double(length(x)),
    x0 = as.double(x0),
    xmin = as.double(xmin),    
    xmax = as.double(xmax),
    hx = as.double(hx),
    px = as.character(px),
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

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densHistogram.x
