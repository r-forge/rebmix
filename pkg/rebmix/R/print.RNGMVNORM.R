print.RNGMVNORM <- function(x, ...)
{
  if (missing(x) || (class(x) != "RNGMVNORM")) {
    stop(sQuote("x"), " object of class RNGMVNORM is requested!", call. = FALSE)
  }

  cat("$w", "\n", sep = "")

  print(x$w, quote = FALSE, ...)

  cat("\n", sep = "")

  cat("$ymin", "\n", sep = "")

  print(x$ymin, quote = FALSE, ...)

  cat("\n", sep = "")

  cat("$ymax", "\n", sep = "")

  print(x$ymax, quote = FALSE, ...)

  cat("\n", sep = "")

  cat("attr(\"class\")", "\n", sep = "")

  print(attr(x, "class"), ...)
  
  rm(list = ls())
} ## print.RNGMVNORM
