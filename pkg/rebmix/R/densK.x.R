.densK.x <- function(v, x, k, hx)
{
  output <- .C(C_RdensKX,
    v = as.integer(v),
    x = as.double(x),
    k = as.double(k),
    y = double(v),
    hx = as.double(hx),
    error = character(1),
    PACKAGE = "rebmix")

  error <- unlist(strsplit(output$error, "\n"));
      
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

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densK.x
