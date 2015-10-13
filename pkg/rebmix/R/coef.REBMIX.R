coef.REBMIX <- function(x, pos = 1, ...)
{
  if (missing(x) || (class(x) != "REBMIX")) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  if ((pos < 1) || (pos > nrow(x$summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x$summary), "!", call. = FALSE)
  }
  
  cat(paste("$w", "\n", sep = ""))

  print(as.number(x$w[[pos]]), quote = FALSE, ...)  

  cat(paste("\n", sep = ""))
  
  cat(paste("$Theta", "\n", sep = ""))

  print(x$Theta[[pos]], quote = FALSE, ...) 
  
  rm(list = ls())   
} ## coef.REBMIX
