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

  print(x$w[[pos]], quote = FALSE, ...)  

  cat("\n", sep = "")
  
  cat("$Theta", "\n", sep = "")

  print(x$Theta[[pos]], quote = FALSE, ...)  

  cat("$summary", "\n", sep = "")

  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(x$summary), nomatch = 0)

  print(x$summary[pos, p], quote = FALSE, ...)

  cat("\n", sep = "")

  cat("attr(\"class\")", "\n", sep = "")

  print(attr(x, "class"), ...)
  
  rm(list = ls())
} ## print.REBMIX
