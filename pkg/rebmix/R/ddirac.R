ddirac <- function(x, location)
{
  f <- array(data = 0.0, dim = length(x), dimnames = NULL)

  f[x == location] <- 1.0

  rm(list = ls()[!(ls() %in% c("f"))])
 
  return(f)
} ## ddirac
