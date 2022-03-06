setMethod("histogram",
          signature(Dataset = "data.frame"),
function(x, Dataset, K, ymin, ymax, shrink, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  # Dataset.

  if (missing(Dataset) || (length(Dataset) == 0)) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }

  if (!is.data.frame(Dataset)) {
    stop(sQuote("Dataset"), " data frame is requested!", call. = FALSE)
  }

  d <- ncol(Dataset)

  if (d < 1) {
    stop(sQuote("Dataset"), " number of columns in data frame must be greater than 0!", call. = FALSE)
  }

  ny <- nrow(Dataset)

  if (ny < 1) {
    stop(sQuote("Dataset"), " number of rows in data frame must be greater than 0!", call. = FALSE)
  }
  
  # x
  
  if (missing(x) || is.null(x)) {
    x <- new("Histogram", 
      Y = as.data.frame(matrix(0.0, nrow = 1, ncol = d + 1)), 
      K = K, 
      ymin = ymin, 
      ymax = ymax)
      
    nz <- prod(x@K)
    
    names <- names(x@Y)
    
    x@Y <- as.data.frame(matrix(0.0, nrow = nz, ncol = d + 1))

    colnames(x@Y) <- names    
  }
  else {
    if (class(x) != "Histogram") {
      stop(sQuote("x"), " object of class Histogram is requested!", call. = FALSE)
    }
    
    if (ncol(x@Y) != d + 1) {
      stop(sQuote("Dataset"), " number of columns in data frame must equal ", ncol(x@Y), "!", call. = FALSE)
    }
    
    nz <- prod(x@K)
    
    if (nrow(x@Y) != nz) {
      stop(sQuote("x"), " has allready bin shrunk!", call. = FALSE)
    }    
  }     
  
  y <- as.matrix(Dataset)
  
  output <- .C(C_Rhistogram,
    K = as.integer(x@K),
    y0 = as.double(x@y0),    
    h = as.double(x@h),
    d = as.integer(d),
    ny = as.integer(ny),
    y = as.double(y),
    nz = as.integer(nz),
    z = as.double(unlist(x@Y)),
    shrink = as.integer(shrink),    
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in histogram!", call. = FALSE); return(NA)
  }
  
  dim(output$z) <- c(nz, d + 1)
  
  if (shrink) output$z <- output$z[1:output$nz, ]
  
  colnames(output$z) <- colnames(x@Y)
  
  x@Y <- as.data.frame(output$z)
  
  output <- x
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)  
}) ## histogram
