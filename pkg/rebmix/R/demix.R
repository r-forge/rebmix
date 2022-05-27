setMethod("demix",
          signature(x = "REBMIX"),
function(x, pos, variables, ...)
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

  if (as.character(class(Dataset)) == "data.frame") {
    Y.type <- 0
    
    Dataset <- as.matrix(Dataset)
    
    d <- ncol(Dataset)
  }  
  else
  if (as.character(class(Dataset)) == "Histogram") {
    Y.type <- 1
    
    Dataset <- as.matrix(Dataset@Y)
    
    d <- ncol(Dataset) - 1
  }
  
  dini <- d; variables <- eval(variables)

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

  Names <- names(x@Theta[[pos]])

  k <- as.numeric(x@summary[pos, "v/k"])

  Names <- names(x@summary)

  if (Y.type == 0) {
    Dataset <- Dataset[, variables]
    
    Preprocessing <- x@summary[pos, "Preprocessing"] 
      
    if (Preprocessing == .rebmix$Preprocessing[1]) {
      h <- x@summary[pos, grep("h", Names)]; h <- h[variables]
      y0 <- x@summary[pos, grep("y0", Names)]; y0 <- y0[variables]
      ymin <- x@summary[pos, grep("ymin", Names)]; ymin <- ymin[variables]
      ymax <- x@summary[pos, grep("ymax", Names)]; ymax <- ymax[variables]    

      output <- .C(C_RPreprocessingHMIX,
        h = as.double(h),
        y0 = as.double(y0),
        ymin = as.double(ymin),
        ymax = as.double(ymax),      
        k = as.integer(k),
        n = as.integer(n),
        d = as.integer(d),
        x = as.double(Dataset),
        y = double(n * (d + 1)),
        error = integer(1),
        PACKAGE = "rebmix")

      if (output$error == 1) {
        stop("in RPreprocessingHMIX!", call. = FALSE); return(NA)
      }

      length(output$y) <- output$k * (output$d + 1); dim(output$y) <- c(output$k, output$d + 1)

      output$y[, d + 1] <- output$y[, d + 1] / prod(output$h) / n

      output <- as.data.frame(output$y, stringsAsFactors = FALSE)

      colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")
    }
    else
    if (Preprocessing == .rebmix$Preprocessing[2]) {
      h <- x@summary[pos, grep("h", Names)]; h <- h[variables]

      output <- .C(C_RPreprocessingKDEMIX,
        h = as.double(h),
        n = as.integer(n),
        d = as.integer(d),
        x = as.double(Dataset),
        y = double(n * (d + 2)),
        error = integer(1),
        PACKAGE = "rebmix")

      if (output$error == 1) {
        stop("in RPreprocessingKDEMIX!", call. = FALSE); return(NA)
      }

      dim(output$y) <- c(n, d + 2)

      output$y[, d + 2] <- output$y[, d + 2] / prod(output$h) / n
 
      output <- as.data.frame(output$y[, -(d + 1)], stringsAsFactors = FALSE)

      colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")
    }
    else
    if (Preprocessing == .rebmix$Preprocessing[3]) {
      h <- x@summary[pos, grep("h", Names)]; h <- h[variables]

      output <- .C(C_RPreprocessingKNNMIX,
        k = as.integer(k),
        h = as.double(h),
        n = as.integer(n),
        d = as.integer(d),
        x = as.double(Dataset),
        y = double(n * (d + 3)),
        error = integer(1),
        PACKAGE = "rebmix")

      if (output$error == 1) {
        stop("in RPreprocessingKNNMIX!", call. = FALSE); return(NA)
      }

      dim(output$y) <- c(n, d + 3)

      output$y[, d + 2] <- k / output$y[, d + 2] / n

      output <- as.data.frame(output$y[, c(-(d + 1), -(d + 3))], stringsAsFactors = FALSE)

      colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")
    }
  }
  else
  if (Y.type == 1) {
    Dataset <- Dataset[, c(variables, dini + 1)]
  
    h <- x@summary[pos, grep("h", Names)]; h <- h[variables]
    
    output <- .C(C_RPreprocessingKMIX,
      h = as.double(h),
      d = as.integer(d),
      n = as.integer(n),
      x = as.double(Dataset),
      error = integer(1),
      PACKAGE = "rebmix")    
    
    if (output$error == 1) {
      stop("in RPreprocessingKMIX!", call. = FALSE); return(NA)
    }
    
    dim(output$x) <- c(n, d + 1); 
    
    output$x <- output$x[1:output$n, ]
    
    dim(output$x) <- c(output$n, d + 1)    

    output$x[, d + 1] <- output$x[, d + 1] / prod(output$h) / sum(output$x[, d + 1])

    output <- as.data.frame(output$x, stringsAsFactors = FALSE)

    colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f") 
  }

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## demix

setMethod("demix",
          signature(x = "REBMVNORM"),
function(x, pos, variables, ...)
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

  if (as.character(class(Dataset)) == "data.frame") {
    Y.type <- 0
    
    Dataset <- as.matrix(Dataset)
    
    d <- ncol(Dataset)
  }  
  else
  if (as.character(class(Dataset)) == "Histogram") {
    Y.type <- 1
    
    Dataset <- as.matrix(Dataset@Y)
    
    d <- ncol(Dataset) - 1
  }

  dini <- d; variables <- eval(variables)
  
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

  k <- as.numeric(x@summary[pos, "v/k"])

  Names <- names(x@summary)

  if (Y.type == 0) {
    Dataset <- Dataset[, variables]    
    
    Preprocessing <- x@summary[pos, "Preprocessing"]  
  
    if (Preprocessing == .rebmix$Preprocessing[1]) {
      h <- x@summary[pos, grep("h", Names)]; h <- h[variables]
      y0 <- x@summary[pos, grep("y0", Names)]; y0 <- y0[variables]
      ymin <- x@summary[pos, grep("ymin", Names)]; ymin <- ymin[variables]
      ymax <- x@summary[pos, grep("ymax", Names)]; ymax <- ymax[variables]

      output <- .C(C_RPreprocessingHMVNORM,
        h = as.double(h),
        y0 = as.double(y0),
        ymin = as.double(ymin),
        ymax = as.double(ymax),
        k = as.integer(k),
        n = as.integer(n),
        d = as.integer(d),
        x = as.double(Dataset),
        y = double(n * (d + 1)),
        error = integer(1),
        PACKAGE = "rebmix")

      if (output$error == 1) {
        stop("in RPreprocessingHMVNORM!", call. = FALSE); return(NA)
      }

      length(output$y) <- output$k * (output$d + 1); dim(output$y) <- c(output$k, output$d + 1)

      output$y[, d + 1] <- output$y[, d + 1] / prod(output$h) / n

      output <- as.data.frame(output$y, stringsAsFactors = FALSE)

      colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")
    }
    else
    if (Preprocessing == .rebmix$Preprocessing[2]) {
      h <- x@summary[pos, grep("h", Names)]; h <- h[variables]

      output <- .C(C_RPreprocessingKDEMVNORM,
        h = as.double(h),
        n = as.integer(n),
        d = as.integer(d),
        x = as.double(Dataset),
        y = double(n * (d + 2)),
        error = integer(1),
        PACKAGE = "rebmix")

      if (output$error == 1) {
        stop("in RPreprocessingKDEMVNORM!", call. = FALSE); return(NA)
      }

      dim(output$y) <- c(n, d + 2)

      output$y[, d + 2] <- output$y[, d + 2] / prod(output$h) / n

      output <- as.data.frame(output$y[, -(d + 1)], stringsAsFactors = FALSE)

      colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")
    }
    else
    if (Preprocessing == .rebmix$Preprocessing[3]) {
      h <- x@summary[pos, grep("h", Names)]; h <- h[variables]

      output <- .C(C_RPreprocessingKNNMVNORM,
        k = as.integer(k),
        h = as.double(h),
        n = as.integer(n),
        d = as.integer(d),
        x = as.double(Dataset),
        y = double(n * (d + 3)),
        error = integer(1),
        PACKAGE = "rebmix")

      if (output$error == 1) {
        stop("in RPreprocessingKNNMVNORM!", call. = FALSE); return(NA)
      }

      dim(output$y) <- c(n, d + 3)

      output$y[, d + 2] <- k / output$y[, d + 2] / n

      output <- as.data.frame(output$y[, c(-(d + 1), -(d + 3))], stringsAsFactors = FALSE)

      colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")
    }
  }
  else
  if (Y.type == 1) {
    Dataset <- Dataset[, c(variables, dini + 1)]
  
    h <- x@summary[pos, grep("h", Names)]; h <- h[variables]
    
    output <- .C(C_RPreprocessingKMIX,
      h = as.double(h),
      d = as.integer(d),
      n = as.integer(n),
      x = as.double(Dataset),
      error = integer(1),
      PACKAGE = "rebmix")    
    
    if (output$error == 1) {
      stop("in RPreprocessingKMIX!", call. = FALSE); return(NA)
    }
    
    dim(output$x) <- c(n, d + 1); 
    
    output$x <- output$x[1:output$n, ]
    
    dim(output$x) <- c(output$n, d + 1)    

    output$x[, d + 1] <- output$x[, d + 1] / prod(output$h) / sum(output$x[, d + 1])

    output <- as.data.frame(output$x, stringsAsFactors = FALSE)

    colnames(output) <- c(paste("x", if (dini > 1) variables else "", sep = ""), "f")     
  }

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## demix
