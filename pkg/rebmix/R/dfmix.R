setMethod("dfmix",
          signature(x = "REBMIX"),
function(x,
  Dataset,
  pos,
  variables, ...)
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
      else
      if (pdf[[i]][j] == .rebmix$pdf[8]) {
        fi <- fi * dunif(as.numeric(Dataset[, j]), min = as.numeric(theta1[[i]][j]), max = as.numeric(theta2[[i]][j]), ...)
      }
      else      
      if (pdf[[i]][j] == .rebmix$pdf[9]) {
        output <- .C(C_RvonMisesPdf,
          n = as.integer(n),
          y = as.double(Dataset[, j]),
          Mean = as.double(theta1[[i]][j]),
          Kappa = as.double(theta2[[i]][j]),
          f = double(n),
          PACKAGE = "rebmix")

        fi <- fi * output$f
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[10]) {
        output <- .C(C_RGumbelPdf,
          n = as.integer(n),
          y = as.double(Dataset[, j]),
          Mean = as.double(theta1[[i]][j]),
          Sigma = as.double(theta2[[i]][j]),
          Xi = as.double(theta3[[i]][j]),
          f = double(n),
          PACKAGE = "rebmix")

        fi <- fi * output$f      
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
function(x,
  Dataset,
  pos,
  variables, ...)
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

      output <- .C(C_RMvtNormalPdf,
        n = as.integer(n),
        X = as.double(unlist(Dataset[, variables])),
        d = as.integer(length(mean)),
        Mean = as.double(mean),
        Sigma = as.double(sigma),
        f = double(n),
        error = integer(9),
        PACKAGE = "rebmix")
        
      error <- error.to.string(output$error);
      
      if (error[1] != "") {
        stop(error[1], call. = FALSE); return(NA)
      }
   
      if (error[2] != "") {
        warning(error[2], call. = FALSE, immediate. = TRUE)
      }  
    
      if (error[3] != "") {
        warning(error[3], call. = FALSE, immediate. = TRUE)
      }        

      fi <- output$f        
    }

    f <- f + as.numeric(w[i]) * fi
  }

  output <- as.data.frame(cbind(Dataset[, variables], f), stringsAsFactors = FALSE)

  colnames(output) <- c(paste("x", if (d > 1) variables else "", sep = ""), "f")

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## dfmix
