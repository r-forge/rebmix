print.boot.REBMIX <- function(x, ...)
{
  if (missing(x) || (class(x) != "boot.REBMIX")) {
    stop(sQuote("x"), " object of class boot.REBMIX is requested!", call. = FALSE)
  }
  
  names <- names(x)
  
  for (i in 1:length(names)) {
    cat("$", names[i], "\n", sep = "")

    print(x[[names[i]]], quote = FALSE, ...)

    cat("\n", sep = "")  
  }

  cat("attr(\"class\")", "\n", sep = "")

  print(attr(x, "class"), ...)
  
  rm(list = ls())
} ## print.boot.REBMIX
