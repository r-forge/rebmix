setMethod("pfmix",
          signature(x = "REBMIX"),
function(x,
  Dataset,
  pos,
  variables,
  lower.tail,
  log.p, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }  

  if (missing(Dataset)) {
    Dataset <- x@Dataset[[pos]]  
  
    if (as.character(class(Dataset)) == "Histogram") {
      d <- ncol(Dataset@Y) - 1
   
      Dataset <- as.data.frame(Dataset@Y[, 1:d])
    }
  }

  if (!is.data.frame(Dataset)) {
    stop(sQuote("Dataset"), " data frame is requested!", call. = FALSE)
  }

  d <- length(x@Variables); variables <- eval(variables)

  if (ncol(Dataset) != d) {
    stop(sQuote("Dataset"), " number of columns in data frame must equal ", d, "!", call. = FALSE)
  }

  n <- nrow(Dataset)

  if (n < 1) {
    stop(sQuote("Dataset"), " number of rows in data frame must be greater than 0!", call. = FALSE)
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

  if (!is.logical(lower.tail)) {
    stop(sQuote("lower.tail"), " logical is requested!", call. = FALSE)
  }

  if (!is.logical(log.p)) {
    stop(sQuote("log.p"), " logical is requested!", call. = FALSE)
  }

  w <- x@w[[pos]]

  c <- length(w)

  Theta <- x@Theta[[pos]]

  Names <- names(Theta)

  pdf <- Theta[grep("pdf", Names)]

  theta1 <- Theta[grep("theta1", Names)]

  theta2 <- Theta[grep("theta2", Names)]
  
  theta3 <- Theta[grep("theta3", Names)]

  f <- array(data = 0.0, dim = n, dimnames = NULL)

  for (i in 1:c) {
    fi <- rep(1.0, n)

    for (j in variables) {
      if (pdf[[i]][j] == .rebmix$pdf[1]) {
        fi <- fi * pnorm(as.numeric(Dataset[, j]), mean = as.numeric(theta1[[i]][j]), sd = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[2]) {
        fi <- fi * plnorm(as.numeric(Dataset[, j]), meanlog = as.numeric(theta1[[i]][j]), sdlog = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[3]) {
        fi <- fi * pweibull(as.numeric(Dataset[, j]), scale = as.numeric(theta1[[i]][j]), shape = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[4]) {
        fi <- fi * pbinom(as.integer(Dataset[, j]), size = as.integer(theta1[[i]][j]), prob = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[5]) {
        fi <- fi * ppois(as.integer(Dataset[, j]), lambda = as.numeric(theta1[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[6]) {
        fi <- fi * pdirac(as.numeric(Dataset[, j]), location = as.numeric(theta1[[i]][j]))
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[7]) {
        fi <- fi * pgamma(as.numeric(Dataset[, j]), scale = as.numeric(theta1[[i]][j]), shape = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[8]) {
        fi <- fi * punif(as.numeric(Dataset[, j]), min = as.numeric(theta1[[i]][j]), max = as.numeric(theta2[[i]][j]), ...)
      }
      else      
      if (pdf[[i]][j] == .rebmix$pdf[9]) {
        output <- .C(C_RvonMisesCdf,
          n = as.integer(n),
          y = as.double(Dataset[, j]),
          Mean = as.double(theta1[[i]][j]),
          Kappa = as.double(theta2[[i]][j]),
          F = double(n),
          error = integer(1),
          PACKAGE = "rebmix")
          
        if (output$error == 1) {
          stop("in pfmix!", call. = FALSE); return(NA)
        }          

        fi <- fi * output$F
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[10]) {
        output <- .C(C_RGumbelCdf,
          n = as.integer(n),
          y = as.double(Dataset[, j]),
          Mean = as.double(theta1[[i]][j]),
          Sigma = as.double(theta2[[i]][j]),
          Xi = as.double(theta3[[i]][j]),
          F = double(n),
          PACKAGE = "rebmix")

        fi <- fi * output$F      
      }        
    }

    f <- f + as.numeric(w[i]) * fi
  }

  output <- as.data.frame(cbind(Dataset[, variables], f), stringsAsFactors = FALSE)

  colnames(output) <- c(paste("x", if (d > 1) variables else "", sep = ""), "F")

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## pfmix

setMethod("pfmix",
          signature(x = "REBMVNORM"),
function(x,
  Dataset,
  pos,
  variables,
  lower.tail,
  log.p, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMVNORM is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }  

  if (missing(Dataset)) {
    Dataset <- x@Dataset[[pos]]  
  
    if (as.character(class(Dataset)) == "Histogram") {
      d <- ncol(Dataset@Y) - 1
   
      Dataset <- as.data.frame(Dataset@Y[, 1:d])
    }
  }
  
  if (!is.data.frame(Dataset)) {
    stop(sQuote("Dataset"), " data frame is requested!", call. = FALSE)
  }

  d <- length(x@Variables); variables <- eval(variables)

  if (ncol(Dataset) != d) {
    stop(sQuote("Dataset"), " number of columns in data frame must equal ", d, "!", call. = FALSE)
  }

  n <- nrow(Dataset)

  if (n < 1) {
    stop(sQuote("Dataset"), " number of rows in data frame must be greater than 0!", call. = FALSE)
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

  if (!is.logical(lower.tail)) {
    stop(sQuote("lower.tail"), " logical is requested!", call. = FALSE)
  }

  if (!is.logical(log.p)) {
    stop(sQuote("log.p"), " logical is requested!", call. = FALSE)
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
    fi <- rep(0.0, n)

    if (all(pdf[[i]] == .rebmix$pdf[1])) {
      mean <- as.numeric(theta1[[i]][variables])

      sigma <- matrix(theta2[[i]], ncol = d, byrow = TRUE)

      sigma <- sigma[variables, variables]

      for (j in 1:n) {
#       fi[j] <- pmvnorm(upper = as.numeric(Dataset[j, variables]), mean = mean, sigma = sigma, ...)
        fi[j] <- 0.0
      }
    }

    f <- f + as.numeric(w[i]) * fi
  }

  output <- as.data.frame(cbind(Dataset[, variables], f), stringsAsFactors = FALSE)

  colnames(output) <- c(paste("x", if (d > 1) variables else "", sep = ""), "F")

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## pfmix

