dfmix <- function(x = NULL, w = NULL, Theta = NULL, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  if (is.null(x)) {
    stop(sQuote("x"), " must not be NULL!", call. = FALSE)
  }

  if ((!is.numeric(x)) && (!is.data.frame(x))) {
    stop(sQuote("x"), " numeric or data frame is requested!", call. = FALSE)
  }
  
  x <- rbind(x)
  
  d <- ncol(x)
  n <- nrow(x)  
  
  if (is.null(w)) {
    stop(sQuote("w"), " must not be NULL!", call. = FALSE)
  }  
  
  if ((!is.numeric(w)) && (!is.data.frame(w))) {
    stop(sQuote("w"), " numeric or data frame is requested!", call. = FALSE)
  }
  
  if (!all(w > 0.0)) {
    stop("all ", sQuote("w"), " must be greater than 0.0!", call. = FALSE)
  }
  
  if (is.null(Theta)) {
    stop(sQuote("Theta"), " must not be NULL!", call. = FALSE)
  }

  if ((!is.character(Theta)) && (!is.data.frame(Theta))) {
    stop(sQuote("Theta"), " character matrix or data frame is requested!", call. = FALSE)
  }
  
  if (length(w) != ncol(Theta)) {
    stop("number of columns in ", sQuote("w"), " and ", sQuote("Theta"), " must match!", call. = FALSE)
  }

  nrow <- nrow(Theta)
  ncol <- ncol(Theta)

  c <- ncol

  w <- rbind(w)

  Theta <- rbind(Theta)

  pdf <- array(data = NA, dim = c(nrow, ncol), dimnames = NULL)
  theta1 <- array(data = 0.0, dim = c(nrow, ncol), dimnames = NULL)
  theta2 <- array(data = 0.0, dim = c(nrow, ncol), dimnames = NULL)

  for (i in 1:ncol) {
    j <- 1; k <- 1
    
    while (j < nrow) {
      pdf[k, i] <- match.arg(Theta[j, i], .rebmix$pdf)

      if (pdf[k, i] %in% .rebmix$pdf[.rebmix$pdf.nargs == 2]) {
        theta1[k, i] <- as.numeric(Theta[j + 1, i])
        theta2[k, i] <- as.numeric(Theta[j + 2, i])

        j <- j + 3; k <- k + 1
      }
      else
      if (pdf[k, i] %in% .rebmix$pdf[.rebmix$pdf.nargs == 1]) {
        theta1[k, i] <- as.numeric(Theta[j + 1, i])
        theta2[k, i] <- as.numeric(0.0)

        j <- j + 2; k <- k + 1
      }
    }
  }
  
  d <- k - 1; c <- ncol
  
  pdf <- pdf[1:d, ]; dim(pdf) <- c(d, c)
  theta1 <- theta1[1:d, ]; dim(theta1) <- c(d, c)
  theta2 <- theta2[1:d, ]; dim(theta2) <- c(d, c)
      
  f <- array(data = 0.0, dim = n, dimnames = NULL)

  for (i in 1:c) {
    fi <- rep(1.0, n)
    
    for (j in 1:d) {
      if (pdf[j, i] == .rebmix$pdf[1]) {
        fi <- fi * dnorm(as.numeric(x[, j]), mean = as.numeric(theta1[j, i]), sd = as.numeric(theta2[j, i]), ...)
      }
      else
      if (pdf[j, i] == .rebmix$pdf[2]) {
        fi <- fi * dlnorm(as.numeric(x[, j]), meanlog = as.numeric(theta1[j, i]), sdlog = as.numeric(theta2[j, i]), ...)
      }
      else
      if (pdf[j, i] == .rebmix$pdf[3]) {
        fi <- fi * dweibull(as.numeric(x[, j]), scale = as.numeric(theta1[j, i]), shape = as.numeric(theta2[j, i]), ...)
      }
      else
      if (pdf[j, i] == .rebmix$pdf[4]) {
        fi <- fi * dbinom(as.integer(x[, j]), size = as.integer(theta1[j, i]), prob = as.numeric(theta2[j, i]), ...)
      }
      else
      if (pdf[j, i] == .rebmix$pdf[5]) {
        fi <- fi * dpois(as.integer(x[, j]), lambda = as.numeric(theta1[j, i]), ...)
      }
      else
      if (pdf[j, i] == .rebmix$pdf[6]) {
        fi <- fi * ddirac(as.numeric(x[, j]), location = as.numeric(theta1[j, i]))
      }
      else
      if (pdf[j, i] == .rebmix$pdf[7]) {
        fi <- fi * dgamma(as.numeric(x[, j]), scale = as.numeric(theta1[j, i]), shape = as.numeric(theta2[j, i]), ...)
      }      
    }
    
    f <- f + as.numeric(w[i]) * fi
  }
  
  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("f"))])
  
  invisible(f)
} ## dfmix
