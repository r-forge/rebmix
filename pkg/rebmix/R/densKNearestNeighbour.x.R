.densKNearestNeighbour.x <- function(x, k, hx, npts)
{
  output <- .C(C_RdensKNearestNeighbourX,
    n = as.integer(length(x)),
    x = as.double(x),
    y = double(length(x)),
    k = as.integer(k),
    hx = as.double(hx),
    error = character(3075),
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
} ## .densKNearestNeighbour.x
