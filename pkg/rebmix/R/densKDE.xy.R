.densKDE.xy <- function(x, y, hx, hy, npts)
{
  output <- .C(C_RdensKDEXY,
    n = as.integer(length(x)),
    x = as.double(x),
    y = as.double(y),
    z = double(length(x)),
    hx = as.double(hx),
    hy = as.double(hy),
    error = character(1),
    PACKAGE = "rebmix")

  error <- unlist(strsplit(output$error, "\n"));
      
  if (error[1] != "") {
    stop(error[1], call. = FALSE); return(NA)
  }
    
  if (error[2] != "") {
    warning(error[2], call. = FALSE, immediate. = TRUE)
  }  
    
  if (error[3] != "") {
    warning(error[3], call. = FALSE, immediate. = TRUE)
  }

  i <- !duplicated(data.frame(output$x, output$y))

  output$x <- output$x[i]
  output$y <- output$y[i]
  output$z <- output$z[i]

  n <- length(output$z)

  if (n > npts) {
    i <- sample.int(n, npts, replace = FALSE, prob = NULL)

    output$x <- output$x[i]
    output$y <- output$y[i]
    output$z <- output$z[i]
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densKDE.xy
