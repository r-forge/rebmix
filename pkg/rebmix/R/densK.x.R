.densK.x <- function(v, x, k, hx)
{
  output <- .C(C_RdensKX,
    v = as.integer(v),
    x = as.double(x),
    k = as.integer(k),
    y = double(v),
    hx = as.double(hx),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densK.x!", call. = FALSE); return(NA)
  }

  length(output$x) <- output$v
  length(output$y) <- output$v

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densK.x
