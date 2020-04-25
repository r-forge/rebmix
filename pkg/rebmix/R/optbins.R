### Panic Branislav & Marko Nagode.  

setMethod("optbins",
          signature(Dataset = "list"),
function(Dataset, Rule = "Knuth equal", y0, ymin, ymax, kmin, kmax, ...)
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
  
  # Rule.

  if (!is.character(Rule)) {
    stop(sQuote("Rule"), " character is requested!", call. = FALSE)
  }

  Rule <- match.arg(Rule, .optbins$Rule, several.ok = FALSE)
  
  # kmin.   
  
  if (missing(kmin) || length(kmin) == 0) {
    if (Rule %in% c(.optbins$Rule[4], .optbins$Rule[5])) {
      stop(sQuote("kmin"), " must not be empty!", call. = FALSE)
    }
  }
  else {
    if (!is.wholenumber(kmin)) {
      stop(sQuote("kmin"), " integer is requested!", call. = FALSE)
    }

    length(kmin) <- 1

    if (kmin < 1) {
      stop(sQuote("kmin"), " must be greater than 0!", call. = FALSE)
    }  
  }
    
  # kmax.
    
  if (missing(kmax) || length(kmax) == 0) {
    if (Rule %in% c(.optbins$Rule[4], .optbins$Rule[5])) {  
      stop(sQuote("kmax"), " must not be empty!", call. = FALSE)
    }
  }
  else {
    if (!is.wholenumber(kmax)) {
      stop(sQuote("kmax"), " integer is requested!", call. = FALSE)
    }

    length(kmax) <- 1

    if (kmax < kmin) {
      stop(sQuote("kmax"), " must be greater or equal than ", kmin, "!", call. = FALSE)
    }
  }    
    
  # y0.

  if (missing(y0) || (length(y0) == 0)) {
    y0 <- numeric()
  }
  else {
    if (!is.numeric(y0)) {
      stop(sQuote("y0"), " numeric is requested!", call. = FALSE)
    }

    if (length(y0) != d) {
      stop("lengths of ", sQuote("y0"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
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
  
  output <- matrix(0, nrow = length(Dataset), ncol = d)

  for (i in 1:length(Dataset)) {
    X <- as.matrix(Dataset[[i]])

    n <- nrow(X)
    d <- ncol(X)

    temp <- .C(C_Roptbins,
      d = as.integer(d),
      n = as.integer(n),
      X = as.double(X),
      Rule = as.character(Rule),
      length.y0 = as.integer(length(y0)),
      y0 = as.double(y0),
      length.ymin = as.integer(length(ymin)),
      ymin = as.double(ymin),
      length.ymax = as.integer(length(ymax)),
      ymax = as.double(ymax),            
      kmin = as.integer(kmin),
      kmax = as.integer(kmax),
      opt.k = integer(d),
      error = integer(1),
      PACKAGE = "rebmix")

    if (temp$error == 1) {
      stop("in optbins!", call. = FALSE); return(NA)
    }
    
    output[i, ] <- temp$opt.k
  }
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)  
}) ## optbins

### End
