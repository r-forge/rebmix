pdirac <- function(q, location, lower.tail = TRUE, log.p = FALSE)
{
  f <- array(data = 0.0, dim = length(q), dimnames = NULL)

  if (lower.tail == TRUE) {
    f[q >= location] <- 1.0
  }
  else {
    f[q < location] <- 1.0
  }
  
  if (log.p == TRUE) {
    f <- log(f)
  }

  rm(list = ls()[!(ls() %in% c("f"))])
 
  return(f)
} ## pdirac
