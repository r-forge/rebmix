setMethod("coef",
          signature(x = "REBMIX"),
function(x, pos, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  w <- matrix(x@w[[pos]], nrow = 1)
  
  c <- ncol(w)  
  
  rownames(w) <- "w"
  colnames(w) <- paste("comp", if (x@summary[pos, "c"] > 1) 1:x@summary[pos, "c"] else "", sep = "") 
  
  print(w, quote = FALSE, ...)
  
  d <- length(x@Variables)   
  
  names <- names(x@Theta[[pos]]) 

  Names <- names[grep("theta1", names, fixed = TRUE)]
  
  theta1 <- matrix(unlist(x@Theta[[pos]][Names]), ncol = d, byrow = TRUE)
  
  rownames(theta1) <- Names
  colnames(theta1) <- paste(1:d, sep = "")
  
  print(theta1, quote = FALSE, ...) 

  Names <- names[grep("theta2", names, fixed = TRUE)]
  
  theta2 <- matrix(unlist(x@Theta[[pos]][Names]), ncol = d, byrow = TRUE)
  
  rownames(theta2) <- Names
  colnames(theta2) <- paste(1:d, sep = "")
  
  print(theta2, quote = FALSE, ...)
  
  rm(list = ls())    
}) ## coef

setMethod("coef",
          signature(x = "REBMVNORM"),
function(x, pos, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMVNORM is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  w <- matrix(x@w[[pos]], nrow = 1)
  
  c <- ncol(w)  
  
  rownames(w) <- "w"
  colnames(w) <- paste("comp", if (x@summary[pos, "c"] > 1) 1:x@summary[pos, "c"] else "", sep = "") 
  
  print(w, quote = FALSE, ...)
  
  d <- length(x@Variables)   
  
  names <- names(x@Theta[[pos]]) 

  Names <- names[grep("theta1", names, fixed = TRUE)]
  
  theta1 <- matrix(unlist(x@Theta[[pos]][Names]), ncol = d, byrow = TRUE)
  
  rownames(theta1) <- Names
  colnames(theta1) <- paste(1:d, sep = "")
  
  print(theta1, quote = FALSE, ...) 

  Names <- names[grep("theta2", names, fixed = TRUE)]
  
  theta2 <- matrix(unlist(x@Theta[[pos]][Names]), ncol = d * d, byrow = TRUE)
  
  rownames(theta2) <- Names
  colnames(theta2) <- paste(rep(1:d, each = d), rep(1:d, d), sep = "-") 
  
  print(theta2, quote = FALSE, ...)
  
  rm(list = ls())   
}) ## coef
