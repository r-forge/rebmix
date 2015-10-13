summary.boot.REBMIX <- function(x, ...)
{
  if (missing(x) || (class(x) != "boot.REBMIX")) {
    stop(sQuote("x"), " object of class boot.REBMIX is requested!", call. = FALSE)
  }
  
  w.cv <- matrix(x$w.cv, nrow = 1) 
  
  rownames(w.cv) <- "w.cv"
  colnames(w.cv) <- paste("comp", if (x$c.mode > 1) 1:x$c.mode else "", sep = "") 
  
  print(w.cv, quote = FALSE, ...)

  cat("\n", sep = "")  
  
  names <- names(x)
  
  names <- names[grep("theta", names, fixed = TRUE)]
  names <- names[grep(".cv", names, fixed = TRUE)]  
  
  theta.cv <- NULL
  
  for (i in 1:length(names)) {
    theta.cv <- c(theta.cv, x[[names[i]]])
  }
  
  theta.cv <- matrix(theta.cv, ncol = length(x[[names[1]]]), byrow = TRUE)
  
  rownames(theta.cv) <- names
  colnames(theta.cv) <- names(x[[names[1]]])
  
  print(theta.cv, quote = FALSE, ...)

  cat(paste("Mode probability = ", x$c.prob, " at c = ", x$c.mode, " components.\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
} ## summary.boot.REBMIX
