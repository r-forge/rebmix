summary.REBMIX <- function(object, ...)
{
  if (missing(object) || (class(object) != "REBMIX")) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object$summary), nomatch = 0)

  print(object$summary[p], quote = FALSE, ...)

  cat(paste("Maximum logL = ", object$summary[object$pos, "logL"], " at pos = ", object$pos, ".\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
} ## summary.REBMIX
