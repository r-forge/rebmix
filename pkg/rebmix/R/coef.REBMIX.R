setMethod("coef",
          signature(object = "REBMIX"),
function(object, pos, ...)
{
  if (missing(object) || (class(object) != "REBMIX")) {
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
  
  print(theta, quote = FALSE, ...)  
  
  rm(list = ls())   
}) ## coef
