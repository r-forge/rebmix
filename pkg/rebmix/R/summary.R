setMethod("summary", 
          signature(object = "REBMIX"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object@summary), nomatch = 0)

  print(object@summary[p], quote = FALSE, ...)

  cat(paste("Maximum logL = ", object@summary[object@pos, "logL"], " at pos = ", object@pos, ".\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
}) ## summary

setMethod("summary", 
          signature(object = "REBMVNORM"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMVNORM is requested!", call. = FALSE)
  }
  
  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object@summary), nomatch = 0)

  print(object@summary[p], quote = FALSE, ...)

  cat(paste("Maximum logL = ", object@summary[object@pos, "logL"], " at pos = ", object@pos, ".\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
}) ## summary

setMethod("summary", 
          signature(object = "REBMIX.boot"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX.boot is requested!", call. = FALSE)
  }
  
  w.cv <- matrix(object@w.cv, nrow = 1) 
  
  rownames(w.cv) <- "w.cv"
  colnames(w.cv) <- paste("comp", if (object@c.mode > 1) 1:object@c.mode else "", sep = "") 
  
  print(w.cv, quote = FALSE, ...)

  cat("\n", sep = "")  
  
  names <- names(object@Theta.cv) 

  names <- names[grep("theta", names, fixed = TRUE)]
  
  d <- length(object@x@Variables)

  theta.cv <- matrix(unlist(object@Theta.cv[names]), ncol = d, byrow = TRUE)  
  
  rownames(theta.cv) <- names
  colnames(theta.cv) <- names(object@Theta.cv[[names[1]]])
  
  print(theta.cv, quote = FALSE, ...)

  cat(paste("Mode probability = ", object@c.prob, " at c = ", object@c.mode, " components.\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
}) ## summary

setMethod("summary", 
          signature(object = "REBMVNORM.boot"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMVNORM.boot is requested!", call. = FALSE)
  }
  
  w.cv <- matrix(object@w.cv, nrow = 1) 
  
  rownames(w.cv) <- "w.cv"
  colnames(w.cv) <- paste("comp", if (object@c.mode > 1) 1:object@c.mode else "", sep = "") 
  
  print(w.cv, quote = FALSE, ...)

  cat("\n", sep = "")  
  
  names <- names(object@Theta.cv) 

  names <- names[grep("theta", names, fixed = TRUE)]
  
  d <- length(object@x@Variables)

  theta.cv <- matrix(unlist(object@Theta.cv[names]), ncol = d, byrow = TRUE)  
  
  rownames(theta.cv) <- names
  colnames(theta.cv) <- names(object@Theta.cv[[names[1]]])
  
  print(theta.cv, quote = FALSE, ...)

  cat(paste("Mode probability = ", object@c.prob, " at c = ", object@c.mode, " components.\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
}) ## summary
