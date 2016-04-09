setMethod("pemix",
          signature(x = "REBMIX"),
function(x, 
  pos = 1, 
  variables = numeric(), 
  lower.tail = TRUE, 
  log.p = FALSE, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
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
  
  Dataset <- x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]]
  
  d <- ncol(Dataset); dini <- d
  n <- nrow(Dataset)
  
  if (length(variables) != 0) {
    if (!is.wholenumber(variables)) {
      stop(sQuote("variables"), " integer is requested!", call. = FALSE)
    }

    if ((min(variables) < 1) || (max(variables) > d)) {
      stop(sQuote("variables"), " must be greater than 0 and less or equal than ", d, "!", call. = FALSE)
    }
    
    variables <- unique(variables); d <- length(variables) 
  }
  else {
    variables <- 1:d
  }
  
  Dataset <- as.matrix(Dataset[, variables])
   
  F <- array(data = 0.0, dim = n, dimnames = NULL)

  if (lower.tail == TRUE) {
    for (i in 1:n) {
      F[i] <- sum(apply(Dataset <= Dataset[i,], 1, all))
    }
  }
  else {
    for (i in 1:n) {
      F[i] <- sum(apply(Dataset > Dataset[i,], 1, all))
    }
  }

  F <- F / n

  if (log.p == TRUE) {
    F <- log(F)
  }
  
  output <- as.data.frame(cbind(Dataset, F), stringsAsFactors = FALSE)

  colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "F")  

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## pemix

setMethod("pemix",
          signature(x = "REBMVNORM"),
function(x, 
  pos = 1, 
  variables = numeric(), 
  lower.tail = TRUE, 
  log.p = FALSE, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
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
  
  Dataset <- x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]]
  
  d <- ncol(Dataset); dini <- d
  n <- nrow(Dataset)
  
  if (length(variables) != 0) {
    if (!is.wholenumber(variables)) {
      stop(sQuote("variables"), " integer is requested!", call. = FALSE)
    }

    if ((min(variables) < 1) || (max(variables) > d)) {
      stop(sQuote("variables"), " must be greater than 0 and less or equal than ", d, "!", call. = FALSE)
    }
    
    variables <- unique(variables); d <- length(variables)
  }
  else {
    variables <- 1:d
  }
  
  Dataset <- as.matrix(Dataset[, variables])
   
  F <- array(data = 0.0, dim = n, dimnames = NULL)

  if (lower.tail == TRUE) {
    for (i in 1:n) {
      F[i] <- sum(apply(Dataset <= Dataset[i,], 1, all))
    }
  }
  else {
    for (i in 1:n) {
      F[i] <- sum(apply(Dataset > Dataset[i,], 1, all))
    }
  }

  F <- F / n

  if (log.p == TRUE) {
    F <- log(F)
  }
  
  output <- as.data.frame(cbind(Dataset, F), stringsAsFactors = FALSE)

  colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "F")  

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## pemix
