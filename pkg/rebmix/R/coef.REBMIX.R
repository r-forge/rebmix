coef.REBMIX <- function(x, pos = 1, ...)
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
  
  w <- matrix(x$w[[pos]], nrow = 1) 
  
  rownames(w) <- "w"
  colnames(w) <- paste("comp", if (x$summary[pos, "c"] > 1) 1:x$summary[pos, "c"] else "", sep = "") 
  
  print(w, quote = FALSE, ...)

  cat("\n", sep = "")
  
  names <- names(x$Theta[[pos]])
  
  names <- names[grep("theta", names, fixed = TRUE)]
  
  theta <- NULL
  
  for (i in 1:length(names)) {
    theta <- c(theta, x$Theta[[pos]][[names[i]]])
  }  
  
  theta <- matrix(theta, ncol = length(x$Theta[[pos]][[names[1]]]), byrow = TRUE)
  
  rownames(theta) <- names
  colnames(theta) <- x$Theta[[pos]]$pdf1
  
  print(theta, quote = FALSE, ...)  
  
  rm(list = ls())   
} ## coef.REBMIX
