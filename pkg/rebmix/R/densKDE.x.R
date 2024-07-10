.densKDE.x <- function(x, hx, npts)
{
  output <- .C(C_RdensKDEX,
    n = as.integer(length(x)),
    x = as.double(x),
    y = double(length(x)),
    hx = as.double(hx),
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

  i <- !duplicated(output$x)

  output$x <- output$x[i]
  output$y <- output$y[i]

  n <- length(output$y)

  if (n > npts) {
    i <- sample.int(n, npts, replace = FALSE, prob = NULL)

    output$x <- output$x[i]
    output$y <- output$y[i]
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densKDE.x
