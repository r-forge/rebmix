setMethod("mapclusters",
          signature(x = "RCLRMIX"),
function(x,
  Dataset,
  s, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  Names <- names(x@x@Theta[[x@pos]])

  pdf <- unlist(x@x@Theta[[x@pos]][grep("pdf", Names)])

  theta1 <- unlist(x@x@Theta[[x@pos]][grep("theta1", Names)])

  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(x@x@Theta[[x@pos]][grep("theta2", Names)])

  theta2[is.na(theta2)] <- 0
  
  theta3 <- unlist(x@x@Theta[[x@pos]][grep("theta3", Names)])

  theta3[is.na(theta3)] <- 0  

  c <- length(x@x@w[[x@pos]])

  w <- x@x@w[[x@pos]]

  d <- length(pdf) / c  
  
  # Dataset.
  
  if (missing(Dataset) || (length(Dataset) == 0)) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }  
  
  if (!is.data.frame(Dataset)) {
    stop(sQuote("Dataset"), " data frame is requested!", call. = FALSE)
  }

  if (ncol(Dataset) != d) {
    stop(sQuote("Dataset"), " number of columns in data frame must equal ", d, "!", call. = FALSE)
  }
      
  n <- nrow(Dataset)

  if (n < 1) {
    stop(sQuote("Dataset"), " number of rows in data frame must be greater than 0!", call. = FALSE)
  }
 
  # s. 
  
  s <- eval(s)

  if (!is.wholenumber(s)) {
    stop(sQuote("s"), " integer is requested!", call. = FALSE)
  }

  length(s) <- 1

  if ((s < 1) || (s > c)) {
    stop(sQuote("s"), " must be greater than 0 and less or equal than ", c, "!", call. = FALSE)
  }  

  Dataset <- as.matrix(Dataset)
    
  n <- nrow(Dataset)

  output <- .C(C_RCLRMIX,
    n = n,
    X = as.double(Dataset),
    d = as.integer(d),
    c = as.integer(unlist(c)),
    w = as.double(unlist(w)),
    pdf = as.character(unlist(pdf)),
    theta1 = as.double(unlist(theta1)),
    theta2 = as.double(unlist(theta2)),
    theta3 = as.double(unlist(theta3)),
    Z = integer(n),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in mapclusters!", call. = FALSE); return(NA)
  }
  
  output <- output$Z
  
  unique.output <- unique(output)
  
  set <- which(x@from %in% unique.output)

  from <- x@from[set]; to <- x@to[set]
  
  l <- length(unique.output)

  while (l > s) {
    l <- l - 1
    
    output[output == from[l]] <- to[l]
  }
  
  output <- as.factor(output)

  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
}) ## mapclusters

setMethod("mapclusters",
          signature(x = "RCLRMVNORM"),
function(x,
  Dataset,
  s, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  Names <- names(x@x@Theta[[x@pos]])

  pdf <- unlist(x@x@Theta[[x@pos]][grep("pdf", Names)])

  theta1 <- unlist(x@x@Theta[[x@pos]][grep("theta1", Names)])

  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(x@x@Theta[[x@pos]][grep("theta2", Names)])

  theta2[is.na(theta2)] <- 0

  c <- length(x@x@w[[x@pos]])

  w <- x@x@w[[x@pos]]

  d <- length(pdf) / c  
  
  # Dataset.
  
  if (missing(Dataset) || (length(Dataset) == 0)) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }  
  
  if (!is.data.frame(Dataset)) {
    stop(sQuote("Dataset"), " data frame is requested!", call. = FALSE)
  }

  if (ncol(Dataset) != d) {
    stop(sQuote("Dataset"), " number of columns in data frame must equal ", d, "!", call. = FALSE)
  }
      
  n <- nrow(Dataset)

  if (n < 1) {
    stop(sQuote("Dataset"), " number of rows in data frame must be greater than 0!", call. = FALSE)
  }
 
  # s. 
  
  s <- eval(s)

  if (!is.wholenumber(s)) {
    stop(sQuote("s"), " integer is requested!", call. = FALSE)
  }

  length(s) <- 1

  if ((s < 1) || (s > c)) {
    stop(sQuote("s"), " must be greater than 0 and less or equal than ", c, "!", call. = FALSE)
  }  

  Dataset <- as.matrix(Dataset)
    
  n <- nrow(Dataset)

  output <- .C(C_RCLRMVNORM,
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
    stop("in mapclusters!", call. = FALSE); return(NA)
  }
  
  output <- output$Z
  
  unique.output <- unique(output)
  
  set <- which(x@from %in% unique.output)

  from <- x@from[set]; to <- x@to[set]
  
  l <- length(unique.output)

  while (l > s) {
    l <- l - 1
    
    output[output == from[l]] <- to[l]
  }
  
  output <- as.factor(output)

  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
}) ## mapclusters
