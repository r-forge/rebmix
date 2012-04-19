.dist.xy <- function(x, y)
{
  n <- length(x)
  
  z <- array(data = 0.0, dim = n, dimnames = NULL)

  for (i in 1:n) {
    z[i] <- sum((x <= x[i]) & (y <= y[i]))
  }
  
  z <- z / n

  i <- order(z)
  
  output <- list()

  output$x <- x[i]
  output$y <- y[i]
  output$z <- z[i]

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .dist.xy
