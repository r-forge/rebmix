summary.boot.REBMIX <- function(object, ...)
{
  if (missing(object) || (class(object) != "boot.REBMIX")) {
    stop(sQuote("object"), " object of class boot.REBMIX is requested!", call. = FALSE)
  }
  
  w.cv <- matrix(object$w.cv, nrow = 1) 
  
  rownames(w.cv) <- "w.cv"
  colnames(w.cv) <- paste("comp", if (object$c.mode > 1) 1:object$c.mode else "", sep = "") 
  
  print(w.cv, quote = FALSE, ...)

  cat("\n", sep = "")  
  
  names <- names(object)
  
  names <- names[grep("theta", names, fixed = TRUE)]
  names <- names[grep(".cv", names, fixed = TRUE)]  
  
  theta.cv <- NULL
  
  for (i in 1:length(names)) {
    theta.cv <- c(theta.cv, object[[names[i]]])
  }
  
  theta.cv <- matrix(theta.cv, ncol = length(object[[names[1]]]), byrow = TRUE)
  
  rownames(theta.cv) <- names
  colnames(theta.cv) <- names(object[[names[1]]])
  
  print(theta.cv, quote = FALSE, ...)

  cat(paste("Mode probability = ", object$c.prob, " at c = ", object$c.mode, " components.\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
} ## summary.boot.REBMIX
