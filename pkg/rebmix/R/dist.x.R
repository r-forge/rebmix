.dist.x <- function(x)
{
  n <- length(x)
  
  y <- array(data = 0.0, dim = n, dimnames = NULL)

  for (i in 1:n) {
    y[i] <- sum(x <= x[i])
  }
  
  y <- y / n

  i <- order(y)
  
  output <- list()

  output$x <- x[i]
  output$y <- y[i]

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .dist.x
