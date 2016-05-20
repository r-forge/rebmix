setMethod("combine",
          signature(x = "REBMIX"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  model <- new("REBMIX.combine",
    x = x, 
    pos = pos)

  Dataset <- as.character(model@x@summary[pos, "Dataset"])

  X <- as.matrix(model@x@Dataset[[which(names(model@x@Dataset) == model@x@summary[pos, "Dataset"])]])

  n <- nrow(X)
  d <- ncol(X)

  h <- as.double(model@x@summary[pos, paste("h", if (d > 1) 1:d else "", sep = "")])
  
  Names <- names(model@x@Theta[[pos]])
  
  c <- length(model@x@w[[pos]])

  pdf <- unlist(model@x@Theta[[pos]][grep("pdf", Names)])
  
  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)
  
  theta1 <- unlist(model@x@Theta[[pos]][grep("theta1", Names)])
  
  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(model@x@Theta[[pos]][grep("theta2", Names)])
  
  theta2[is.na(theta2)] <- 0

  length(pdf) <- d

  C <- model@x@summary[pos, "Preprocessing"]

  if (C == .rebmix$Preprocessing[1]) {
    y0 <- as.double(model@x@summary[pos, paste("y0", if (d > 1) 1:d else "", sep = "")])

    output <- .C("RCombineComponentsHMIX",
      h = as.double(h),
      y0 = as.double(y0),
      k = as.integer(model@x@summary[pos, "v/k"]),
      c = as.integer(c),
      w = as.double(model@x@w[[pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(2),
      length.theta = as.integer(c(d, d)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(X),
      F = integer(c),
      T = integer(c),
      EN = double(c),
      ED = double(c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in combine!", call. = FALSE); return(NA)
    }
  } 
  else 
  if (C == .rebmix$Preprocessing[2]) {
    output <- .C("RCombineComponentsPWMIX",
      h = as.double(h),
      c = as.integer(c),
      w = as.double(model@x@w[[pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(2),
      length.theta = as.integer(c(d, d)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(X),
      F = integer(c),
      T = integer(c),
      EN = double(c),
      ED = double(c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in combine!", call. = FALSE); return(NA)
    }
  } 
  else
  if (C == .rebmix$Preprocessing[3]) {
    k <- as.integer(model@x@summary[pos, "v/k"]) 

    output <- .C("RCombineComponentsKNNMIX",
      h = as.double(h),
      k = as.integer(model@x@summary[pos, "v/k"]),
      c = as.integer(c),
      w = as.double(model@x@w[[pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(2),
      length.theta = as.integer(c(d, d)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(X),
      F = integer(c),
      T = integer(c),
      EN = double(c),
      ED = double(c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in combine!", call. = FALSE); return(NA)
    }
  }
  
  model@from <- output$F
  model@to <- output$T
  model@EN <- output$EN
  model@ED <- output$ED
  
  options(digits = digits)    
  
  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## combine

setMethod("combine",
          signature(x = "REBMVNORM"),
function(x, pos, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  model <- new("REBMVNORM.combine",
    x = x, 
    pos = pos)

  Dataset <- as.character(model@x@summary[pos, "Dataset"])

  X <- as.matrix(model@x@Dataset[[which(names(model@x@Dataset) == model@x@summary[pos, "Dataset"])]])

  n <- nrow(X)
  d <- ncol(X)

  h <- as.double(model@x@summary[pos, paste("h", if (d > 1) 1:d else "", sep = "")])
  
  Names <- names(model@x@Theta[[pos]])
  
  c <- length(model@x@w[[pos]])

  pdf <- unlist(model@x@Theta[[pos]][grep("pdf", Names)])
  
  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)
  
  theta1 <- unlist(model@x@Theta[[pos]][grep("theta1", Names)])
  
  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(model@x@Theta[[pos]][grep("theta2", Names)])
  
  theta2[is.na(theta2)] <- 0

  length(pdf) <- d

  C <- model@x@summary[pos, "Preprocessing"]

  if (C == .rebmix$Preprocessing[1]) {
    y0 <- as.double(model@x@summary[pos, paste("y0", if (d > 1) 1:d else "", sep = "")])

    output <- .C("RCombineComponentsHMVNORM",
      h = as.double(h),
      y0 = as.double(y0),
      k = as.integer(model@x@summary[pos, "v/k"]),
      c = as.integer(c),
      w = as.double(model@x@w[[pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(4),
      length.theta = as.integer(c(d, d * d, d * d, 1)),      
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(X),
      F = integer(c),
      T = integer(c),
      EN = double(c),
      ED = double(c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in combine!", call. = FALSE); return(NA)
    }
  } 
  else 
  if (C == .rebmix$Preprocessing[2]) {
    output <- .C("RCombineComponentsPWMVNORM",
      h = as.double(h),
      c = as.integer(c),
      w = as.double(model@x@w[[pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(4),
      length.theta = as.integer(c(d, d * d, d * d, 1)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(X),
      F = integer(c),
      T = integer(c),
      EN = double(c),
      ED = double(c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in combine!", call. = FALSE); return(NA)
    }
  } 
  else
  if (C == .rebmix$Preprocessing[3]) {
    k <- as.integer(model@x@summary[pos, "v/k"]) 

    output <- .C("RCombineComponentsKNNMVNORM",
      h = as.double(h),
      k = as.integer(model@x@summary[pos, "v/k"]),
      c = as.integer(c),
      w = as.double(model@x@w[[pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(4),
      length.theta = as.integer(c(d, d * d, d * d, 1)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(X),
      F = integer(c),
      T = integer(c),
      EN = double(c),
      ED = double(c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in combine!", call. = FALSE); return(NA)
    }
  }
  
  model@from <- output$F
  model@to <- output$T
  model@EN <- output$EN
  model@ED <- output$ED
  
  options(digits = digits)    
  
  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## combine

