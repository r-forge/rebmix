.densK.xy <- function(v, x, y, k, hx, hy)
{
  output <- .C(C_RdensKXY,
    v = as.integer(v),
    x = as.double(x),
    y = as.double(y),
    k = as.integer(k),
    z = double(length(x)),
    hx = as.double(hx),
    hy = as.double(hy),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densK.xy!", call. = FALSE); return(NA)
  }

  length(output$x) <- output$v
  length(output$y) <- output$v
  length(output$z) <- output$v

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densK.xy
