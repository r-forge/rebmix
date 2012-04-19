summary.boot.REBMIX <- function(object, ...)
{
  if (missing(object) || (class(object) != "boot.REBMIX")) {
    stop(sQuote("object"), " object of class boot.REBMIX is requested!", call. = FALSE)
  }
  
  names <- names(object); 
  
  names <- names[grep(".cv", names, fixed = TRUE)]

  names <- names[-which("c.cv" == names)]
  
  theta.cv <- NULL
  
  for (i in 1:length(names)) {
    theta.cv <- rbind(theta.cv, as.numeric(object[[names[i]]]))
  }
  
  theta.cv <- as.data.frame(rbind(theta.cv), stringsAsFactors = FALSE)
  
  rownames(theta.cv) <- names[grep(".cv", names, fixed = TRUE)];
  colnames(theta.cv) <- paste("comp", if (object$c.mode > 1) 1:object$c.mode else "", sep = "")
  
  print(apply(theta.cv, c(1, 2), as.number), quote = FALSE, ...)

  cat(paste("Mode probability = ", as.number(object$c.prob), " at c = ", as.number(object$c.mode), " components.\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
} ## summary.boot.REBMIX
