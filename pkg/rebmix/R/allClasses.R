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
prototype = list(Dataset.name = NULL,
  rseed = -1,
  n = NULL,
  Theta = NULL))                          
                          
setClass("RNGMVNORM", contains = "RNGMIX")

setMethod("show",
          signature(object = "RNGMIX"),
function(object)
{
  if (missing(object) || (class(object) != "RNGMIX")) {
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
  if (missing(object) || (class(object) != "RNGMVNORM")) {
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
  all.length = "list",
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
  if (missing(object) || (class(object) != "REBMIX")) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  if (!is.wholenumber(object@pos)) {
    stop(sQuote("object@pos"), " integer is requested!", call. = FALSE)
  }
  
  length(object@pos) <- 1

  if ((object@pos < 1) || (object@pos > nrow(object@summary))) {
    stop(sQuote("object@pos"), " must be greater than 0 and less or equal than ", nrow(object@summary), "!", call. = FALSE)
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
