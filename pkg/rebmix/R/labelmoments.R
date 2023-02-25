setMethod("labelmoments",
          signature(Zp = "array"),
function(Zp, ...)
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
  
  c <- max(Zp)
  
  if (d == 2) {
    output <- .C(C_RLabelMomentsXY,
      nx = as.integer(nd[2]),
      ny = as.integer(nd[1]),
      Zp = as.integer(Zp),    
      c = as.integer(c),
      N = double(c),
      Mx = double(c),
      My = double(c),
      Mxy = double(c),
      Eud = double(c * c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RLabelMomentsXY!", call. = FALSE); return(NA)
    }
    
    output$Eud <- matrix(output$Eud, nrow = c, ncol = c)
    
    rownames(output$Eud) <- paste(1:c, sep = "")
    colnames(output$Eud) <- paste(1:c, sep = "")
    
    clusters <- which(output$N != 0)
    
    output <- list(Mij = data.frame(l = clusters, 
      n = output$N[clusters], 
      M10 = output$Mx[clusters], 
      M01 = output$My[clusters], 
      M11 = output$Mxy[clusters]),
      Eud = output$Eud[clusters, clusters])
  }
  else {
    output <- .C(C_RLabelMomentsXYZ,
      nx = as.integer(nd[2]),
      ny = as.integer(nd[1]),
      nz = as.integer(nd[3]),
      Zp = as.integer(Zp),    
      c = as.integer(c),
      N = double(c),
      Mx = double(c),
      My = double(c),
      Mz = double(c),
      Mxyz = double(c),
      Eud = double(c * c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RLabelMomentsXYZ!", call. = FALSE); return(NA)
    }
    
    output$Eud <- matrix(output$Eud, nrow = c, ncol = c)
    
    rownames(output$Eud) <- paste(1:c, sep = "")
    colnames(output$Eud) <- paste(1:c, sep = "")
    
    clusters <- which(output$N != 0)   
    
    output <- list(Mij = data.frame(l = clusters, 
      n = output$N[clusters], 
      M100 = output$Mx[clusters], 
      M010 = output$My[clusters], 
      M001 = output$Mz[clusters], 
      M111 = output$Mxyz[clusters]),
      Eud = output$Eud[clusters, clusters])  
  }
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)  
}) ## labelmoments
