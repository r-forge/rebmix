setMethod("coef",
          signature(object = "REBMIX"),
function(object, pos, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(object@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(object@summary), "!", call. = FALSE)
  }
  
  w <- matrix(object@w[[pos]], nrow = 1) 
  
  rownames(w) <- "w"
  colnames(w) <- paste("comp", if (object@summary[pos, "c"] > 1) 1:object@summary[pos, "c"] else "", sep = "") 
  
  print(w, quote = FALSE, ...)

  cat("\n", sep = "")
  
  names <- names(object@Theta[[pos]])
  
  names <- names[grep("theta", names, fixed = TRUE)]
  
  theta <- NULL
  
  for (i in 1:length(names)) {
    theta <- c(theta, object@Theta[[pos]][[names[i]]])
  }  
  
  theta <- matrix(theta, ncol = length(object@Theta[[pos]][[names[1]]]), byrow = TRUE)
  
  rownames(theta) <- names
  colnames(theta) <- object@Theta[[pos]]$pdf1
  
  if (nrow(theta) > ncol(theta)) {
    theta <- t(theta)
  }  
  
  print(theta, quote = FALSE, ...)  
  
  rm(list = ls())   
}) ## coef

setMethod("coef",
          signature(object = "REBMVNORM"),
function(object, pos, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMVNORM is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(object@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(object@summary), "!", call. = FALSE)
  }
  
  w <- matrix(object@w[[pos]], nrow = 1)
  
  c <- ncol(w)  
  
  rownames(w) <- "w"
  colnames(w) <- paste("comp", if (object@summary[pos, "c"] > 1) 1:object@summary[pos, "c"] else "", sep = "") 
  
  print(w, quote = FALSE, ...)

  names <- names(object@Theta[[pos]])
  
  names <- names[grep("theta", names, fixed = TRUE)]
  
  d <- length(object@Theta[[pos]][[names[1]]])
  
  theta <- NULL
  
  for (i in 1:length(names)) {
    theta <- c(theta, object@Theta[[pos]][[names[i]]])
  }  
  
  theta <- matrix(theta, ncol = d, byrow = TRUE)
  
  Names <- array(data = "", dim = nrow(theta), dimnames = NULL)
  
  i <- 1
  
  for (j in 1:c) {
    Names[d * (j - 1) + j] <- names[i]
    
    i <- i + 1
    
    Names[d * (j - 1) + j + 1] <- names[i]    
    
    i <- i + 1
  }
  
  rownames(theta) <- Names
  colnames(theta) <- rep("", d)
  
  print(theta, quote = FALSE, ...)  
  
  rm(list = ls())   
}) ## coef
