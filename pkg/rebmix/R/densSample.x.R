.densSample.x <- function(x, ymin, npts)
{
  i <- !duplicated(x)
  
  output <- list()

  output$x <- x[i]
  
  n <- length(output$x)
  
  output$y <- rep(ymin, n)

  if (n > npts) {
    i <- sample.int(n, npts, replace = FALSE, prob = NULL)

    output$x <- output$x[i]
    output$y <- output$y[i]
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densSample.x
