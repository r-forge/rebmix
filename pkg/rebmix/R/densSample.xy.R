.densSample.xy <- function(x, y, zmin, npts)
{
  i <- !duplicated(data.frame(x, y))
  
  output <- list()

  output$x <- x[i]
  
  n <- length(output$x)
  
  output$y <- y[i]
  
  output$z <- rep(zmin, n)

  if (n > npts) {
    i <- sample.int(n, npts, replace = FALSE, prob = NULL)

    output$x <- output$x[i]
    output$y <- output$y[i]
    output$z <- output$z[i]
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densSample.xy
