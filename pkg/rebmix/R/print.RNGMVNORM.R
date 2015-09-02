print.RNGMVNORM <- function(x, ...)
{
  if (missing(x) || (class(x) != "RNGMVNORM")) {
    stop(sQuote("x"), " object of class RNGMVNORM is requested!", call. = FALSE)
  }

  cat(paste("$w", "\n", sep = ""))

  print(as.number(x$w), quote = FALSE, ...)

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
} ## print.RNGMVNORM
