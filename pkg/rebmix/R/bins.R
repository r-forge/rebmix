### Panic Branislav & Marko Nagode.  
setMethod("bins",
          signature(Dataset = "list"),
function(Dataset, K, ymin, ymax, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  # Dataset.

  if (missing(Dataset) || (length(Dataset) == 0)) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }

  if (!is.list(Dataset)) {
    stop(sQuote("Dataset"), " list is requested!", call. = FALSE)
  }

  if (is.null(names(Dataset))) {
    names(Dataset) <- paste("dataset", 1:length(Dataset), sep = "")
  }

  if (!all(unlist(lapply(Dataset, is.data.frame)))) {
    stop(sQuote("Dataset"), " list of data frames is requested!", call. = FALSE)
  }

  d <- unique(unlist(lapply(Dataset, ncol)))

  if (length(d) != 1) {
    stop(sQuote("Dataset"), " numbers of columns in data frames must be equal!", call. = FALSE)
  }

  if (!all(unlist(lapply(Dataset, ncol)) > 0)) {
    stop(sQuote("Dataset"), " numbers of columns in data frames must be greater than 0!", call. = FALSE)
  }

  for (j in 1:length(Dataset)) {
    Dataset[[j]] <- as.data.frame(Dataset[[j]][complete.cases(Dataset[[j]]), ])
  }

  if (!all(unlist(lapply(Dataset, nrow)) > 1)) {
    stop(sQuote("Dataset"), " numbers of rows in data frames must be greater than 1!", call. = FALSE)
  }
  
  # K.
  
  if (missing(K) || (length(K) == 0)) {
    stop(sQuote("K"), " must not be empty!", call. = FALSE)
  }

  if (!is.matrix(K)) {
    stop(sQuote("K"), " integer matrix is requested!", call. = FALSE)
  }
  
  if (!all(unlist(lapply(K, is.wholenumber) == TRUE))) {
    stop(sQuote("K"), " matrix of integer values is requested!", call. = FALSE)
  }

  if (!all(unlist(lapply(K, function(x) all(x > 0)))) == TRUE) {
    stop("all ", sQuote("K"), " must be greater than 0!", call. = FALSE)
  }
    
  if(ncol(K) != d) {
    stop(sQuote("K"), " number of columns in matrix must equal ", d, "!", call. = FALSE)    
  }
    
  if(nrow(K) != length(Dataset)) {
    stop(sQuote("K"), " number of rows in matrix must equal ", length(Dataset), "!", call. = FALSE)    
  }
  
  # ymin.

  if (missing(ymin) || (length(ymin) == 0)) {
    ymin <- numeric()
  }
  else {
    if (!is.numeric(ymin)) {
      stop(sQuote("ymin"), " numeric is requested!", call. = FALSE)
    }

    if (length(ymin) != d) {
      stop("lengths of ", sQuote("ymin"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
  }
    
  # ymax.

  if (missing(ymax) || (length(ymax) == 0)) {
     ymax <- numeric()
  }
  else {
    if (!is.numeric(ymax)) {
       stop(sQuote("ymax"), " numeric is requested!", call. = FALSE)
    }

    if (length(ymax) != d) {
      stop("lengths of ", sQuote("ymax"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
      
    if (any(ymax <= ymin)) {
      stop(sQuote("ymax"), " must be greater than ", sQuote("ymin"), "!", call. = FALSE)
    }       
  }
  
  output <- list()

  for (i in 1:length(Dataset)) {
    x <- as.matrix(Dataset[[i]])
    
    n <- nrow(x)
    
    temp <- .C(C_Rbins,
      d = as.integer(d),
      n = as.integer(n),
      x = as.double(x),
      length.ymin = as.integer(length(ymin)),
      ymin = as.double(ymin),
      length.ymax = as.integer(length(ymax)),
      ymax = as.double(ymax),            
      k = as.integer(K[i, ]),
      length.y = integer(1),
      y = double(n * (d + 1)),
      error = character(1),
      PACKAGE = "rebmix")

    error <- unlist(strsplit(temp$error, "\n"));
      
    if (error[1] != "") {
      stop(error[1], call. = FALSE); return(NA)
    }
   
    if (error[2] != "") {
      warning(error[2], call. = FALSE, immediate. = TRUE)
    }  
    
    if (error[3] != "") {
      warning(error[3], call. = FALSE, immediate. = TRUE)
    }
    
    length(temp$y) <- temp$length.y * (d + 1); dim(temp$y) <- c(temp$length.y, temp$d + 1)
    
    output[[i]] <- as.data.frame(temp$y, stringsAsFactors = FALSE)
    
    colnames(output[[i]]) <- c(paste("y", if (d > 1) 1:d else "", sep = ""), "k")
  }
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)  
}) ## bins
### End
