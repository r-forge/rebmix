setMethod("labelmoments",
          signature(Zp = "array"),
function(Zp, Sigma, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  # Zp.
  
  if (missing(Zp) || (length(Zp) == 0)) {
    stop(sQuote("Zp"), " must not be empty!", call. = FALSE)
  }

  if (!is.array(Zp)) {
    stop(sQuote("Zp"), " integer array is requested!", call. = FALSE)
  }
  
  nd <- dim(Zp); d <- length(nd)
  
  if ((d < 2) || (d > 3)) {
    stop(sQuote("Zp"), " two or three dimensional array is requested!", call. = FALSE)
  }
  
  # cmax.
  
  if (missing(cmax) || (length(cmax) == 0)) {
    stop(sQuote("cmax"), " must not be empty!", call. = FALSE)
  }  
  
  if (!is.wholenumber(cmax) || (length(cmax) > 1)) {
    stop(sQuote("cmax"), " integer is requested!", call. = FALSE)
  }

  if (cmax < 2) {
    stop(sQuote("cmax"), " must be greater than 1!", call. = FALSE)
  }  
  
  c <- max(cmax, Zp)
  
  # Sigma.
  
  if (!is.numeric(Sigma)) {
    stop(sQuote("Sigma"), " numeric is requested!", call. = FALSE)
  }

  length(Sigma) <- 1

  if (Sigma <= 0.0) {
    stop(sQuote("Sigma"), " must be greater than 0.0!", call. = FALSE)
  }
  
  if (d == 2) {
    output <- .C(C_RLabelMomentsXY,
      nx = as.integer(nd[2]),
      ny = as.integer(nd[1]),
      Zp = as.double(Zp),    
      c = as.integer(c),
      N = double(c),
      Mx = double(c),
      My = double(c),
      Mxy = double(c),
      A = double(c * c),
      Sigma = as.double(Sigma),
      error = character(3075),
      PACKAGE = "rebmix")
      
    error <- unlist(strsplit(output$error, "\n"));

    if (error[1] != "") {
      stop(error[1], call. = FALSE); return(NA)
    }
   
    if (error[2] != "") {
      warning(error[2], call. = FALSE, immediate. = TRUE)
    }  
    
    if (error[3] != "") {
      warning(error[3], call. = FALSE, immediate. = TRUE)
    }
    
    output$A <- matrix(output$A, nrow = c, ncol = c)
    
    rownames(output$A) <- paste(1:c, sep = "")
    colnames(output$A) <- paste(1:c, sep = "")
    
    clusters <- which(output$N != 0.0)
    
    set <- 1:c; set <- set[-clusters]; output$A[set, ] <- output$A[, set] <- NA
    
    output <- list(Mij = data.frame(l = clusters, 
      n = output$N[clusters], 
      M10 = output$Mx[clusters], 
      M01 = output$My[clusters], 
      M11 = output$Mxy[clusters]),
      A = output$A)
  }
  else {
    output <- .C(C_RLabelMomentsXYZ,
      nx = as.integer(nd[2]),
      ny = as.integer(nd[1]),
      nz = as.integer(nd[3]),
      Zp = as.double(Zp),    
      c = as.integer(c),
      N = double(c),
      Mx = double(c),
      My = double(c),
      Mz = double(c),
      Mxyz = double(c),
      A = double(c * c),
      Sigma = as.double(Sigma),
      error = character(3075),
      PACKAGE = "rebmix")
      
    error <- unlist(strsplit(output$error, "\n"));      

    if (error[1] != "") {
      stop(error[1], call. = FALSE); return(NA)
    }
   
    if (error[2] != "") {
      warning(error[2], call. = FALSE, immediate. = TRUE)
    }  
    
    if (error[3] != "") {
      warning(error[3], call. = FALSE, immediate. = TRUE)
    }
    
    output$A <- matrix(output$A, nrow = c, ncol = c)
    
    rownames(output$A) <- paste(1:c, sep = "")
    colnames(output$A) <- paste(1:c, sep = "")
    
    clusters <- which(output$N != 0.0)
    
    set <- 1:c; set <- set[-clusters]; output$A[set, ] <- output$A[, set] <- NA

    output <- list(Mijk = data.frame(l = clusters, 
      n = output$N[clusters], 
      M100 = output$Mx[clusters], 
      M010 = output$My[clusters], 
      M001 = output$Mz[clusters], 
      M111 = output$Mxyz[clusters]),
      A = output$A)
  }
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)  
}) ## labelmoments
