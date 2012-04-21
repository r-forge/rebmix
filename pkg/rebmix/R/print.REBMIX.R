print.REBMIX <- function(x, pos = 1, ...)
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

  print(apply(x$w[[pos]], c(1, 2), as.number), quote = FALSE, ...)

  cat(paste("\n", sep = ""))

  cat(paste("$Theta", "\n", sep = ""))

  print(apply(x$Theta[[pos]], c(1, 2), as.number), quote = FALSE, ...)

  cat(paste("\n", sep = ""))

  cat(paste("$summary", "\n", sep = ""))

  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(x$summary), nomatch = 0)

  print(apply(x$summary[pos, p], c(1, 2), as.number), quote = FALSE, ...)

  cat(paste("\n", sep = ""))

  cat(paste("attr(,\"class\")", "\n", sep = ""))

  print(attr(x, "class"), ...)

  invisible(x)
} ## print.REBMIX
