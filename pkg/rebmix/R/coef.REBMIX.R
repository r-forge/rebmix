coef.REBMIX <- function(object, pos = 1, ...)
{
  if (missing(object) || (class(object) != "REBMIX")) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  if ((pos < 1) || (pos > nrow(object$summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(object$summary), "!", call. = FALSE)
  }
  
  w <- as.data.frame(rbind(object$w[[pos]]), stringsAsFactors = FALSE)
  
  Theta <- as.data.frame(rbind(object$Theta[[pos]]), stringsAsFactors = FALSE)

  print(rbind(apply(w, c(1, 2), as.number), apply(Theta, c(1, 2), as.number)), quote = FALSE, ...)
  
  rm(list = ls())   
} ## coef.REBMIX
