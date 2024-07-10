.densK.xy <- function(v, x, y, k, hx, hy)
{
  output <- .C(C_RdensKXY,
    v = as.integer(v),
    x = as.double(x),
    y = as.double(y),
    k = as.double(k),
    z = double(length(x)),
    hx = as.double(hx),
    hy = as.double(hy),
    error = integer(9),
    PACKAGE = "rebmix")

  error <- error.to.string(output$error);
      
  if (error[1] != "") {
    stop(error[1], call. = FALSE); return(NA)
  }
    
  if (error[2] != "") {
    warning(error[2], call. = FALSE, immediate. = TRUE)
  }  
    
  if (error[3] != "") {
    warning(error[3], call. = FALSE, immediate. = TRUE)
  }

  length(output$x) <- output$v
  length(output$y) <- output$v
  length(output$z) <- output$v

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densK.xy
