print.RNGMIX <- function(x, ...)
{
  if (missing(x) || (class(x) != "RNGMIX")) {
    stop(sQuote("x"), " object of class RNGMIX is requested!", call. = FALSE)
  }

  cat(paste("$w", "\n", sep = ""))

  print(apply(x$w, c(1, 2), as.number), quote = FALSE, ...)

  cat(paste("\n", sep = ""))

  cat(paste("$Theta", "\n", sep = ""))

  print(apply(x$Theta, c(1, 2), as.number), quote = FALSE, ...)

  cat(paste("\n", sep = ""))

  cat(paste("$ymin", "\n", sep = ""))

  print(as.number(x$ymin), quote = FALSE, ...)

  cat(paste("\n", sep = ""))

  cat(paste("$ymax", "\n", sep = ""))

  print(as.number(x$ymax), quote = FALSE, ...)

  cat(paste("\n", sep = ""))

  cat(paste("attr(,\"class\")", "\n", sep = ""))

  print(attr(x, "class"), ...)

  invisible(x)
} ## print.RNGMIX
