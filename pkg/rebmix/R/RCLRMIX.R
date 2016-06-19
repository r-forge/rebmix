setMethod("RCLRMIX",
          signature(model = "RCLRMIX"),
function(model, ...)
{
  Names <- names(model@x@Theta[[model@pos]])
  
  pdf <- unlist(model@x@Theta[[model@pos]][grep("pdf", Names)])
    
  theta1 <- unlist(model@x@Theta[[model@pos]][grep("theta1", Names)])
      
  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(model@x@Theta[[model@pos]][grep("theta2", Names)])
      
  theta2[is.na(theta2)] <- 0
  
  c <- length(model@x@w[[model@pos]])

  w <- model@x@w[[model@pos]]
      
  d <- length(pdf) / c
  
  dataset <- as.matrix(model@x@Dataset[[model@pos]])
  
  n <- nrow(dataset)  

  h <- as.double(model@x@summary[model@pos, paste("h", if (d > 1) 1:d else "", sep = "")])
  
  C <- model@x@summary[model@pos, "Preprocessing"]

  if (C == .rebmix$Preprocessing[1]) {
    y0 <- as.double(model@x@summary[model@pos, paste("y0", if (d > 1) 1:d else "", sep = "")])

    output <- .C("RCombineComponentsHMIX",
      h = as.double(h),
      y0 = as.double(y0),
      k = as.integer(model@x@summary[model@pos, "v/k"]),
      c = as.integer(c),
      w = as.double(model@x@w[[model@pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(2),
      length.theta = as.integer(c(d, d)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(dataset),
      F = integer(c),
      T = integer(c),
      EN = double(c),
      ED = double(c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RCLRMIX!", call. = FALSE); return(NA)
    }
  } 
  else 
  if (C == .rebmix$Preprocessing[2]) {
    output <- .C("RCombineComponentsPWMIX",
      h = as.double(h),
      c = as.integer(c),
      w = as.double(model@x@w[[model@pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(2),
      length.theta = as.integer(c(d, d)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(dataset),
      F = integer(c),
      T = integer(c),
      EN = double(c),
      ED = double(c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RCLRMIX!", call. = FALSE); return(NA)
    }
  } 
  else
  if (C == .rebmix$Preprocessing[3]) {
    k <- as.integer(model@x@summary[model@pos, "v/k"]) 

    output <- .C("RCombineComponentsKNNMIX",
      h = as.double(h),
      k = as.integer(model@x@summary[model@pos, "v/k"]),
      c = as.integer(c),
      w = as.double(model@x@w[[model@pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(2),
      length.theta = as.integer(c(d, d)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(dataset),
      F = integer(c),
      T = integer(c),
      EN = double(c),
      ED = double(c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RCLRMIX!", call. = FALSE); return(NA)
    }
  }
  
  model@from <- output$F
  model@to <- output$T
  model@EN <- output$EN
  model@ED <- output$ED

  output <- .C("RCLRMIX",
    n = n,
    X = as.double(dataset),
    d = as.integer(d),
    c = as.integer(unlist(c)),
    w = as.double(unlist(w)),
    pdf = as.character(unlist(pdf)),
    theta1 = as.double(unlist(theta1)),
    theta2 = as.double(unlist(theta2)),
    Z = integer(n),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RCLRMIX!", call. = FALSE); return(NA)
  }
  
  model@Zp <- as.factor(output$Z) 

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RCLRMIX

setMethod("RCLRMIX",
          signature(model = "RCLRMVNORM"),
function(model, ...)
{
  Names <- names(model@x@Theta[[model@pos]])
    
  pdf <- unlist(model@x@Theta[[model@pos]][grep("pdf", Names)])
    
  theta1 <- unlist(model@x@Theta[[model@pos]][grep("theta1", Names)])
      
  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(model@x@Theta[[model@pos]][grep("theta2", Names)])
      
  theta2[is.na(theta2)] <- 0

  c <- length(model@x@w[[model@pos]])

  w <- model@x@w[[model@pos]]
      
  d <- length(pdf) / c
  
  dataset <- as.matrix(model@x@Dataset[[model@pos]])
  
  n <- nrow(dataset)
  
  h <- as.double(model@x@summary[model@pos, paste("h", if (d > 1) 1:d else "", sep = "")])
  
  C <- model@x@summary[model@pos, "Preprocessing"]

  if (C == .rebmix$Preprocessing[1]) {
    y0 <- as.double(model@x@summary[model@pos, paste("y0", if (d > 1) 1:d else "", sep = "")])

    output <- .C("RCombineComponentsHMVNORM",
      h = as.double(h),
      y0 = as.double(y0),
      k = as.integer(model@x@summary[model@pos, "v/k"]),
      c = as.integer(c),
      w = as.double(model@x@w[[model@pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(4),
      length.theta = as.integer(c(d, d * d, -d * d, -1)),      
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(dataset),
      F = integer(c),
      T = integer(c),
      EN = double(c),
      ED = double(c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RCLRMIX!", call. = FALSE); return(NA)
    }
  } 
  else 
  if (C == .rebmix$Preprocessing[2]) {
    output <- .C("RCombineComponentsPWMVNORM",
      h = as.double(h),
      c = as.integer(c),
      w = as.double(model@x@w[[model@pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(4),
      length.theta = as.integer(c(d, d * d, -d * d, -1)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(dataset),
      F = integer(c),
      T = integer(c),
      EN = double(c),
      ED = double(c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RCLRMIX!", call. = FALSE); return(NA)
    }
  } 
  else
  if (C == .rebmix$Preprocessing[3]) {
    k <- as.integer(model@x@summary[model@pos, "v/k"]) 

    output <- .C("RCombineComponentsKNNMVNORM",
      h = as.double(h),
      k = as.integer(model@x@summary[model@pos, "v/k"]),
      c = as.integer(c),
      w = as.double(model@x@w[[model@pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(4),
      length.theta = as.integer(c(d, d * d, -d * d, -1)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(dataset),
      F = integer(c),
      T = integer(c),
      EN = double(c),
      ED = double(c),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RCLRMIX!", call. = FALSE); return(NA)
    }
  }
  
  model@from <- output$F
  model@to <- output$T
  model@EN <- output$EN
  model@ED <- output$ED

  output <- .C("RCLRMVNORM",
    n = n,
    X = as.double(dataset),
    d = as.integer(d),
    c = as.integer(unlist(c)),
    w = as.double(unlist(w)),
    pdf = as.character(unlist(pdf)),
    theta1 = as.double(unlist(theta1)),
    theta2 = as.double(unlist(theta2)),
    Z = integer(n),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RCLRMIX!", call. = FALSE); return(NA)
  }
  
  model@Zp <- as.factor(output$Z)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RCLRMIX

setMethod("RCLRMIX",
          signature(model = "ANY"),
function(model,
  x,
  pos, 
  Zt, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  message("RCLRMIX Version 2.8.2")
 
  flush.console()
  
  model <- new(model,
    x = x,
    pos = pos,
    Zt = Zt)
     
  model <- RCLRMIX(model = model, ...)
  
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RCLRMIX
