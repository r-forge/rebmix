setMethod("demix",
          signature(x = "REBMIX"),
function(x, 
  pos = 1, 
  variables = numeric(), ...)
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
  
  Dataset <- as.matrix(x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]])
  
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
  
  Dataset <- Dataset[, variables]
   
  Preprocessing <- x@summary[pos, "Preprocessing"]
  
  Names <- names(x@Theta[[pos]])  
  
  pdf <- unlist(x@Theta[[pos]][grep("pdf", Names)])
  
  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)
  
  pdf <- pdf[variables]
  
  k <- as.numeric(x@summary[pos, "v/k"])
  
  Names <- names(x@summary)
  
  if (Preprocessing == .rebmix$Preprocessing[1]) {
    h <- x@summary[pos, grep("h", Names)]; h <- h[variables] 
    y0 <- x@summary[pos, grep("y0", Names)]; y0 <- y0[variables]

    output <- .C("RPreprocessingHMIX",
      h = as.double(h),
      y0 = as.double(y0),
      length.pdf = as.integer(d),
      pdf = as.character(pdf),
      k = as.integer(k),
      n = as.integer(n),
      d = as.integer(d),
      x = as.double(unlist(Dataset)),
      y = double(n * (d + 1)),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in preprocessing!", call. = FALSE); return(NA)
    }

    length(output$y) <- output$k * (output$d + 1); dim(output$y) <- c(output$k, output$d + 1)

    output$y[, d + 1] <- output$y[, d + 1] / prod(output$h) / n

    output <- as.data.frame(output$y, stringsAsFactors = FALSE)

    colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")
  } 
  else 
  if (Preprocessing == .rebmix$Preprocessing[2]) {
    h <- x@summary[pos, grep("h", Names)]; h <- h[variables]   
      
    output <- .C("RPreprocessingPWMIX",
      h = as.double(h),
      n = as.integer(n),
      d = as.integer(d),
      x = as.double(unlist(Dataset)),
      y = double(n * (d + 2)),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in preprocessing!", call. = FALSE); return(NA)
    }

    dim(output$y) <- c(n, d + 2)

    output$y[, d + 2] <- output$y[, d + 2] / prod(output$h) / n

    output <- as.data.frame(output$y[, -(d + 1)], stringsAsFactors = FALSE)

    colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")
  } 
  else
  if (Preprocessing == .rebmix$Preprocessing[3]) {
    h <- x@summary[pos, grep("h", Names)]; h <- h[variables] 

    output <- .C("RPreprocessingKNNMIX",
      k = as.integer(k),
      h = as.double(h),
      n = as.integer(n),
      d = as.integer(d),
      x = as.double(unlist(Dataset)),
      y = double(n * (d + 3)),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in preprocessing!", call. = FALSE); return(NA)
    }

    dim(output$y) <- c(n, d + 3)

    output$y[, d + 2] <- k / output$y[, d + 2] / n

    output <- as.data.frame(output$y[, c(-(d + 1), -(d + 3))], stringsAsFactors = FALSE)

    colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")
  }
  
  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## demix

setMethod("demix",
          signature(x = "REBMVNORM"),
function(x, 
  pos = 1, 
  variables = numeric(), ...)
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
  
  Dataset <- as.matrix(x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]])
  
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
  
  Dataset <- Dataset[, variables]
   
  Preprocessing <- x@summary[pos, "Preprocessing"]
  
  Names <- names(x@Theta[[pos]])  
  
  pdf <- unlist(x@Theta[[pos]][grep("pdf", Names)])
  
  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)
  
  pdf <- pdf[variables]
  
  k <- as.numeric(x@summary[pos, "v/k"])
  
  Names <- names(x@summary)
  
  if (Preprocessing == .rebmix$Preprocessing[1]) {
    h <- x@summary[pos, grep("h", Names)]; h <- h[variables] 
    y0 <- x@summary[pos, grep("y0", Names)]; y0 <- y0[variables]

    output <- .C("RPreprocessingHMVNORM",
      h = as.double(h),
      y0 = as.double(y0),
      length.pdf = as.integer(d),
      pdf = as.character(pdf),
      k = as.integer(k),
      n = as.integer(n),
      d = as.integer(d),
      x = as.double(unlist(Dataset)),
      y = double(n * (d + 1)),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in preprocessing!", call. = FALSE); return(NA)
    }

    length(output$y) <- output$k * (output$d + 1); dim(output$y) <- c(output$k, output$d + 1)

    output$y[, d + 1] <- output$y[, d + 1] / prod(output$h) / n

    output <- as.data.frame(output$y, stringsAsFactors = FALSE)

    colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")
  } 
  else 
  if (Preprocessing == .rebmix$Preprocessing[2]) {
    h <- x@summary[pos, grep("h", Names)]; h <- h[variables]   
      
    output <- .C("RPreprocessingPWMVNORM",
      h = as.double(h),
      n = as.integer(n),
      d = as.integer(d),
      x = as.double(unlist(Dataset)),
      y = double(n * (d + 2)),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in preprocessing!", call. = FALSE); return(NA)
    }

    dim(output$y) <- c(n, d + 2)

    output$y[, d + 2] <- output$y[, d + 2] / prod(output$h) / n

    output <- as.data.frame(output$y[, -(d + 1)], stringsAsFactors = FALSE)

    colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")
  } 
  else
  if (Preprocessing == .rebmix$Preprocessing[3]) {
    h <- x@summary[pos, grep("h", Names)]; h <- h[variables] 

    output <- .C("RPreprocessingKNNMVNORM",
      k = as.integer(k),
      h = as.double(h),
      n = as.integer(n),
      d = as.integer(d),
      x = as.double(unlist(Dataset)),
      y = double(n * (d + 3)),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in preprocessing!", call. = FALSE); return(NA)
    }

    dim(output$y) <- c(n, d + 3)

    output$y[, d + 2] <- k / output$y[, d + 2] / n

    output <- as.data.frame(output$y[, c(-(d + 1), -(d + 3))], stringsAsFactors = FALSE)

    colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")
  }
  
  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## demix
