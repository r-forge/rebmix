summary.REBMIX <- function(x, ...)
{
  if (missing(x) || (class(x) != "REBMIX")) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(x$summary), nomatch = 0)

  print(x$summary[p], quote = FALSE, ...)

  cat(paste("Maximum logL = ", x$summary[x$pos, "logL"], " at pos = ", x$pos, ".\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
} ## summary.REBMIX
