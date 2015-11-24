# Class RNGMIX

validity.RNGMIX <- function(object)
{ 
  # Dataset.name.
  
  if (length(object@Dataset.name) == 0) {
    stop(sQuote("Dataset.name"), " must not be empty!", call. = FALSE)
  }
  
  # rseed.    
  
  if (!is.wholenumber(object@rseed)) {
    stop(sQuote("rseed"), " integer is requested!", call. = FALSE)
  }
  
  if (object@rseed > -1) {
    stop(sQuote("rseed"), " must be less than 0!", call. = FALSE)
  }
  
  return(TRUE)
} ## validity.RNGMIX

setClass(Class = "RNGMIX",
slots = c(Dataset.name = "character",
  rseed = "numeric",
  n = "numeric",
  Theta = "list",
  Dataset = "list",
  w = "numeric",
  Variables = "character",
  ymin = "numeric",
  ymax = "numeric"),
validity = validity.RNGMIX)

setMethod("initialize", "RNGMIX", 
function(.Object, ...,
  n,
  Theta)
{
  # n.
  
  if (length(n) == 0) {
    stop(sQuote("n"), " must not be empty!", call. = FALSE)
  }

  if (!is.wholenumber(n)) {
    stop(sQuote("n"), " integer is requested!", call. = FALSE)
  }
  
  if (!all(n > 0)) {
    stop("all ", sQuote("n"), " must be greater than 0!", call. = FALSE)
  }
  
  c <- length(n)
  
  # Theta.  
  
  if (length(Theta) == 0) {
    stop(sQuote("Theta"), " must not be empty!", call. = FALSE)
  }
  
  Names <- names(Theta)
  
  j <- 0; length.pdf <- length(Theta[[1]])
  
  for (i in grep("pdf", Names)) {  
    pdf <- as.character(Theta[[i]])
  
    pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)
    
    if (length(pdf) != length.pdf) {
      stop("lengths of ", sQuote("pdfi"), " must be equal!", call. = FALSE)
    }    
  
    Theta[[i]] <- pdf; j <- j + 1
  }

  if ((length.pdf > 1) && (j != c)) {
    stop(sQuote("pdfi"), " and ", sQuote("n"), " must match!", call. = FALSE)
  } 
  
  j <- 0; length.theta1 <- length(Theta[[2]])
  
  for (i in grep("theta1", Names)) {  
    theta1 <- as.numeric(Theta[[i]])
   
    if (length(theta1) != length.theta1) {
      stop("lengths of ", sQuote("theta1.i"), " must be equal!", call. = FALSE)
    }    

    j <- j + 1
  }

  if ((length.pdf > 1) && (j != c)) {
    stop(sQuote("theta1.i"), " and ", sQuote("n"), " must match!", call. = FALSE)
  } 
  
  j <- 0; length.theta2 <- length(Theta[[3]])
  
  for (i in grep("theta2", Names)) {  
    theta2 <- as.numeric(Theta[[i]])
   
    if (length(theta2) != length.theta2) {
      stop("lengths of ", sQuote("theta2.i"), " must be equal!", call. = FALSE)
    }

    j <- j + 1    
  }

  if ((length.pdf > 1) && (j != c)) {
    stop(sQuote("theta2.i"), " and ", sQuote("n"), " must match!", call. = FALSE)
  } 
 
  # Variables.

  for (i in 1:length(.rebmix$pdf)) {
    .Object@Variables[which(pdf == .rebmix$pdf[i])] <- .rebmix$pdf.Variables[i]
  }
  
  callNextMethod(.Object, ...,
    n = n,
    Theta = Theta)
})
                          
setClass("RNGMVNORM", contains = "RNGMIX")

setMethod("show",
          signature(object = "RNGMIX"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RNGMIX is requested!", call. = FALSE)
  }
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")
  
  cat("Slot \"w\":", "\n", sep = "")

  print(object@w, quote = FALSE)

  cat("Slot \"ymin\":", "\n", sep = "")

  print(object@ymin, quote = FALSE)

  cat("Slot \"ymax\":", "\n", sep = "")

  print(object@ymax, quote = FALSE)

  rm(list = ls())
}) ## show

setMethod("show",
          signature(object = "RNGMVNORM"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RNGMVNORM is requested!", call. = FALSE)
  }
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")
  
  cat("Slot \"w\":", "\n", sep = "")

  print(object@w, quote = FALSE)

  cat("Slot \"ymin\":", "\n", sep = "")

  print(object@ymin, quote = FALSE)

  cat("Slot \"ymax\":", "\n", sep = "")

  print(object@ymax, quote = FALSE)

  rm(list = ls())
}) ## show

# Class REBMIX

validity.REBMIX <- function(object)
{ 
  # Dataset.
  
  if (length(object@Dataset) == 0) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }

  if (!all(unlist(lapply(object@Dataset, is.data.frame)))) {
    stop(sQuote("Dataset"), " list of data frames is requested!", call. = FALSE)
  }

  d <- ncol(object@Dataset[[1]])

  if (d < 1) {
    stop(sQuote("d"), " must be greater than 0!", call. = FALSE)
  }

  if (!all(unlist(lapply(object@Dataset, ncol)) == d)) {
    stop(sQuote("Dataset"), " numbers of columns in data frames must be equal!", call. = FALSE)
  }

  if (!all(unlist(lapply(object@Dataset, nrow)) > 1)) {
    stop(sQuote("Dataset"), " numbers of rows in data frames must be greater than 1!", call. = FALSE)
  }

  # cmax.

  if (!is.wholenumber(object@cmax)) {
    stop(sQuote("cmax"), " integer is requested!", call. = FALSE)
  }

  if (object@cmax < 1) {
    stop(sQuote("cmax"), " must be greater than 0!", call. = FALSE)
  }

  # pdf. 

  if (length(object@pdf) != d) {
    stop("lengths of ", sQuote("pdf"), " and ", sQuote("d"), " must match!", call. = FALSE)
  } 

  # theta1.  

  if (length(object@theta1) > 0) {
    object@theta1[is.na(object@theta1)] <- 0
      
    if (length(object@theta1) != d) {
      stop("lengths of ", sQuote("theta1"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
  }

  # theta2.  

  if (length(object@theta2) > 0) {
    object@theta2[is.na(object@theta2)] <- 0
      
    if (length(object@theta2) != d) {
      stop("lengths of ", sQuote("theta2"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
  }

  # K. 

  if (length(object@K) == 0) {
    stop(sQuote("K"), " must not be empty!", call. = FALSE)
  }

  if (is.list(object@K)) {
    for (i in 1:length(object@K)) {
      if (!is.wholenumber(object@K[[i]])) {
        stop(sQuote("K"), " integer vector is requested!", call. = FALSE)
      }

      if (!all(object@K[[i]] > 0)) {
        stop("all ", sQuote("K"), " must be greater than 0!", call. = FALSE)
      }
    }

    if (length(object@K) != length(object@Preprocessing)) {
      stop("lengths of ", sQuote("Preprocessing"), " and ", sQuote("K"), " must match!", call. = FALSE)
    }
  }
  else 
  if (is.numeric(object@K)) {
    if (!is.wholenumber(object@K)) {
      stop(sQuote("K"), " integer vector is requested!", call. = FALSE)
    }

    if (!all(object@K > 0)) {
      stop("all ", sQuote("K"), " must be greater than 0!", call. = FALSE)
    }
  }
  else {
    stop(sQuote("K"), " list of integer vectors or integer vector is requested!", call. = FALSE)
  }

  # y0. 

  if (length(object@y0) > 0) {
    object@y0[is.na(object@y0)] <- 0
      
    if (length(object@y0) != d) {
      stop("lengths of ", sQuote("y0"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }      
  } 

  # ymin.    

  if (length(object@ymin) > 0) {
    object@ymin[is.na(object@ymin)] <- 0
      
    if (length(object@ymin) != d) {
      stop("lengths of ", sQuote("ymin"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }      
  }

  # ymax. 

  if (length(object@ymax) > 0) {
    object@ymax[is.na(object@ymax)] <- 0
      
    if (length(object@ymax) != d) {
      stop("lengths of ", sQuote("ymax"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }      
  }

  if ((object@ar <= 0.0) || (object@ar > 1.0)) {
    stop(sQuote("ar"), " must be greater than 0.0 and less or equal than 1.0!", call. = FALSE)
  }

  return(TRUE)
} ## validity.REBMIX

setClass(Class = "REBMIX",
slots = c(Dataset = "list",
  Preprocessing = "character",
  cmax = "numeric",
  Criterion = "character",
  Variables = "character",
  pdf = "character",
  theta1 = "numeric",
  theta2 = "numeric",
  K = "ANY",
  y0 = "numeric",
  ymin = "numeric",
  ymax = "numeric",                                      
  ar = "numeric",   
  Restraints = "character",
  w = "list",
  Theta = "list",
  summary = "ANY",
  pos = "numeric",
  opt.c = "list",
  opt.IC = "list",
  opt.logL = "list",
  opt.D = "list",
  all.K = "list",
  all.IC = "list"),
prototype = list(pos = 1),
validity = validity.REBMIX) 

setMethod("initialize", "REBMIX", 
function(.Object, ...,
  Dataset,
  Preprocessing,
  Criterion,
  pdf,
  Restraints)
{
  # Dataset.
  
  if (is.null(names(Dataset))) {
    names(Dataset) <- paste("dataset", 1:length(Dataset), sep = "")
  }

  # Preprocessing.

  Preprocessing <- match.arg(Preprocessing, .rebmix$Preprocessing, several.ok = TRUE)

  # Criterion.

  Criterion <- match.arg(Criterion, .rebmix$Criterion, several.ok = TRUE)

  # pdf.

  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)

  # Restraints.

  Restraints <- match.arg(Restraints, .rebmix$Restraints, several.ok = FALSE)
  
  # Variables.

  for (i in 1:length(.rebmix$pdf)) {
    .Object@Variables[which(pdf == .rebmix$pdf[i])] <- .rebmix$pdf.Variables[i]
  }
  
  callNextMethod(.Object, ...,
    Dataset = Dataset,
    Preprocessing = Preprocessing,
    Criterion = Criterion,
    pdf = pdf,
    Restraints = Restraints)
})
                    
setClass("REBMVNORM", contains = "REBMIX")

setMethod("show",
          signature(object = "REBMIX"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  if (!is.wholenumber(object@pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(object@pos) <- 1

  if ((object@pos < 1) || (object@pos > nrow(object@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(object@summary), "!", call. = FALSE)
  }  
  
  cat("An object of class ", "\"", class(object), "\"", "\n", sep = "")
  
  cat("Slot \"w\":", "\n", sep = "")

  print(object@w[[object@pos]], quote = FALSE)

  cat("Slot \"Theta\":", "\n", sep = "")

  print(object@Theta[[object@pos]], quote = FALSE)

  cat("Slot \"summary\":", "\n", sep = "")
  
  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object@summary), nomatch = 0)

  print(object@summary[object@pos, p], quote = FALSE)  

  rm(list = ls())
}) ## show

## Class REBMIX.boot

setClass("REBMIX.boot",
slots = c(x = "REBMIX",
  pos = "numeric",
  Bootstrap = "character",
  B = "numeric", 
  n = "numeric",
  replace = "logical", 
  prob = "numeric",
  c = "numeric",
  c.se = "numeric",
  c.cv = "numeric",
  c.mode = "numeric",
  c.prob = "numeric",
  w = "matrix",
  w.se = "numeric",
  w.cv = "numeric",
  Theta = "list",
  Theta.se = "list",
  Theta.cv = "list"),
prototype = list(pos = 1))

setMethod("initialize", "REBMIX.boot", 
function(.Object, ...,
  x,
  pos,
  Bootstrap,
  B,
  n)
{
  # x

  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  # pos.

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }  
  
  # Bootstrap.
  
  Bootstrap <- match.arg(Bootstrap, .rebmix.boot$Bootstrap, several.ok = FALSE) 

  # B.
  
  if (!is.wholenumber(B)) {
    stop(sQuote("B"), " integer is requested!", call. = FALSE)
  }
  
  length(B) <- 1

  if (B < 1) {
    stop(sQuote("B"), " must be greater than 0!", call. = FALSE)
  }
  
  # n.
  
  nmax <- nrow(as.matrix(x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]]))
  
  if (length(n) == 0) {
    n <- nmax
  }
  else {
    if (!is.wholenumber(n)) {
      stop(sQuote("n"), " integer is requested!", call. = FALSE)
    }
  
    if ((n < 1) || (n > nmax)) {
      stop(sQuote("n"), " must be greater than 0 and less or equal than ", nmax, "!", call. = FALSE)
    }
  }    
  
  callNextMethod(.Object, ...,
    x = x,
    pos = pos,
    Bootstrap = Bootstrap,
    B = B,
    n = n)
})



