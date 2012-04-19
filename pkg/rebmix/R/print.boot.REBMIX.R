print.boot.REBMIX <- function(x, ...)
{
  if (missing(x) || (class(x) != "boot.REBMIX")) {
    stop(sQuote("x"), " object of class boot.REBMIX is requested!", call. = FALSE)
  }
  
  names <- names(x); names <- names[order(names)]
  
  cat(paste("$c", "\n", sep = ""))

  print(as.number(x$c), quote = FALSE, ...)

  cat(paste("\n", sep = ""))  
  
  cat(paste("$c.mode", "\n", sep = ""))

  print(as.number(x$c.mode), quote = FALSE, ...)

  cat(paste("\n", sep = ""))
  
  cat(paste("$c.prob", "\n", sep = ""))

  print(as.number(x$c.prob), quote = FALSE, ...)

  cat(paste("\n", sep = ""))  
  
  for (i in grep(".se", names, fixed = TRUE)) {
    cat(paste("$", names[i], "\n", sep = ""))

    print(as.number(x[[names[i]]]), quote = FALSE, ...)

    cat(paste("\n", sep = ""))  
  }
  
  for (i in grep(".cv", names, fixed = TRUE)) {
    cat(paste("$", names[i], "\n", sep = ""))

    print(as.number(x[[names[i]]]), quote = FALSE, ...)

    cat(paste("\n", sep = ""))  
  }  

  cat(paste("attr(,\"class\")", "\n", sep = ""))

  print(attr(x, "class"), ...)

  invisible(x)
} ## print.boot.REBMIX
