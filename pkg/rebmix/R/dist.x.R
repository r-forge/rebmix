.dist.x <- function(x, k, npts)
{
  set <- order(x)
  
  x <- x[set]
  k <- k[set]

  n <- length(x); y <- array()
  
  if (is.null(k)) {
    y <- seq(from = 1 / n, to = 1.0, length = n)
  }
  else {
    y[1] <- k[1]
    
    for (i in 2:n) {
      y[i] <- y[i - 1] + k[i]
    }
    
    y <- y / y[n]
  }

  if (n > npts) {
    i <- sample.int(n, npts, replace = FALSE, prob = NULL)

    x <- x[i]
    y <- y[i]
  }

  output <- list()

  output$x <- x
  output$y <- y

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .dist.x
