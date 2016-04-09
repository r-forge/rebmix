setMethod("dfmix",
          signature(x = "REBMIX"),
function(x = NULL,
  Dataset = data.frame(),
  pos = 1, 
  variables = numeric(), ...)          
{
  digits <- getOption("digits"); options(digits = 15)
  
  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  if (missing(Dataset)) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }
  
  if (!is.data.frame(Dataset)) {
    stop(sQuote("Dataset"), " data frame is requested!", call. = FALSE)
  }

  d <- length(x@Variables)

  if (ncol(Dataset) != d) {
    stop(sQuote("Dataset"), " number of columns in data frame must equal ", d, "!", call. = FALSE)
  }

  n <- nrow(Dataset)

  if (n < 1) {
    stop(sQuote("Dataset"), " number of rows in data frame must be greater than 1!", call. = FALSE)
  }  

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  if (length(variables) != 0) {
    if (!is.wholenumber(variables)) {
      stop(sQuote("variables"), " integer is requested!", call. = FALSE)
    }

    if ((min(variables) < 1) || (max(variables) > d)) {
      stop(sQuote("variables"), " must be greater than 0 and less or equal than ", d, "!", call. = FALSE)
    }
    
    variables <- unique(variables)   
  }
  else {
    variables <- 1:d
  }
  
  w <- x@w[[pos]]
  
  c <- length(w)
  
  Theta <- x@Theta[[pos]]
  
  Names <- names(Theta)
  
  pdf <- Theta[grep("pdf", Names)]

  theta1 <- Theta[grep("theta1", Names)]
  
  theta2 <- Theta[grep("theta2", Names)]

  f <- array(data = 0.0, dim = n, dimnames = NULL)

  for (i in 1:c) {
    fi <- rep(1.0, n)
    
    for (j in variables) {
      if (pdf[[i]][j] == .rebmix$pdf[1]) {
        fi <- fi * dnorm(as.numeric(Dataset[, j]), mean = as.numeric(theta1[[i]][j]), sd = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[2]) {
        fi <- fi * dlnorm(as.numeric(Dataset[, j]), meanlog = as.numeric(theta1[[i]][j]), sdlog = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[3]) {
        fi <- fi * dweibull(as.numeric(Dataset[, j]), scale = as.numeric(theta1[[i]][j]), shape = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[4]) {
        fi <- fi * dbinom(as.integer(Dataset[, j]), size = as.integer(theta1[[i]][j]), prob = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[5]) {
        fi <- fi * dpois(as.integer(Dataset[, j]), lambda = as.numeric(theta1[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[6]) {
        fi <- fi * ddirac(as.numeric(Dataset[, j]), location = as.numeric(theta1[[i]][j]))
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[7]) {
        fi <- fi * dgamma(as.numeric(Dataset[, j]), scale = as.numeric(theta1[[i]][j]), shape = as.numeric(theta2[[i]][j]), ...)
      }      
    }
    
    f <- f + as.numeric(w[i]) * fi
  }
  
  output <- as.data.frame(cbind(Dataset[, variables], f), stringsAsFactors = FALSE)

  colnames(output) <- c(paste("x", if (d > 1) variables else "", sep = ""), "f")    
  
  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
}) ## dfmix

setMethod("dfmix",
          signature(x = "REBMVNORM"),
function(x = NULL,
  Dataset = data.frame(),
  pos = 1, 
  variables = numeric(), ...)          
{
  digits <- getOption("digits"); options(digits = 15)
  
  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMVNORM is requested!", call. = FALSE)
  }
  
  if (missing(Dataset)) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }
  
  if (!is.data.frame(Dataset)) {
    stop(sQuote("Dataset"), " data frame is requested!", call. = FALSE)
  }

  d <- length(x@Variables)

  if (ncol(Dataset) != d) {
    stop(sQuote("Dataset"), " number of columns in data frame must equal ", d, "!", call. = FALSE)
  }

  n <- nrow(Dataset)

  if (n < 1) {
    stop(sQuote("Dataset"), " number of rows in data frame must be greater than 1!", call. = FALSE)
  }  

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  if (length(variables) != 0) {
    if (!is.wholenumber(variables)) {
      stop(sQuote("variables"), " integer is requested!", call. = FALSE)
    }

    if ((min(variables) < 1) || (max(variables) > d)) {
      stop(sQuote("variables"), " must be greater than 0 and less or equal than ", d, "!", call. = FALSE)
    }
    
    variables <- unique(variables)  
  }
  else {
    variables <- 1:d
  }
  
  w <- x@w[[pos]]
  
  c <- length(w)
  
  Theta <- x@Theta[[pos]]
  
  Names <- names(Theta)
  
  pdf <- Theta[grep("pdf", Names)]

  theta1 <- Theta[grep("theta1", Names)]
  
  theta2 <- Theta[grep("theta2", Names)]

  f <- array(data = 0.0, dim = n, dimnames = NULL)

  for (i in 1:c) {
    fi <- 0.0
    
    if (all(pdf[[i]] == .rebmix$pdf[1])) {
      mean <- as.numeric(theta1[[i]][variables])
    
      sigma <- matrix(theta2[[i]], ncol = d, byrow = TRUE)
      
      sigma <- sigma[variables, variables]
          
      fi <- dmvnorm(as.matrix(Dataset[, variables]), mean = mean, sigma = sigma, ...)
    }
    
    f <- f + as.numeric(w[i]) * fi
  }
  
  output <- as.data.frame(cbind(Dataset[, variables], f), stringsAsFactors = FALSE)

  colnames(output) <- c(paste("x", if (d > 1) variables else "", sep = ""), "f")    
  
  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
}) ## dfmix
