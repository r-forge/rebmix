.densParzenWindow.xy <- function(x, y, hx, hy)
{
  output <- .C("RdensParzenWindowXY",
    n = as.integer(length(x)),
    x = as.double(x),
    y = as.double(y),
    z = double(length(x)),
    hx = as.double(hx),
    hy = as.double(hy),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densParzenWindow.xy!", call. = FALSE); return(NA)
  }

  i <- order(output$z)

  output$x <- output$x[i]
  output$y <- output$y[i]
  output$z <- output$z[i]

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densParzenWindow.xy
