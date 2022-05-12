.densK.x <- function(v, x, k, hx)
{
  output <- .C(C_RdensKX,
    v = as.integer(v),
    x = as.double(x),
    k = as.double(k),
    y = double(v),
    hx = as.double(hx),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RdensKX!", call. = FALSE); return(NA)
  }

  length(output$x) <- output$v
  length(output$y) <- output$v

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densK.x
