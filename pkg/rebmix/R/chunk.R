setMethod("chunk", 
          signature(x = "RCLS"),
function(x, variables)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class RCLS is requested!", call. = FALSE)
  }

  d <- ncol(x@test); variables <- eval(variables)
  
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
  
  for (i in 1:x@s) {
    x@train[[i]] <- as.data.frame(x@train[[i]][, variables])
  }
  
  x@test <- as.data.frame(x@test[, variables])

  rm(list = ls()[!(ls() %in% c("x"))]) 

  invisible(x)
}) ## chunk
