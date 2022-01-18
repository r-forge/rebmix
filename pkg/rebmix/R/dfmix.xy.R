.dfmix.xy <- function(x, y, w, xTheta, yTheta, ...)
{
  n <- length(x)

  f <- array(data = 0.0, dim = n, dimnames = NULL)

  for (i in 1:length(w)) {
    if (xTheta[[i]]$pdf == .rebmix$pdf[1]) {
      fix <- dnorm(as.numeric(x), mean = as.numeric(xTheta[[i]]$theta1), sd = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[2]) {
      fix <- dlnorm(as.numeric(x), meanlog = as.numeric(xTheta[[i]]$theta1), sdlog = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[3]) {
      fix <- dweibull(as.numeric(x), scale = as.numeric(xTheta[[i]]$theta1), shape = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[4]) {
      fix <- dbinom(as.integer(x), size = as.integer(xTheta[[i]]$theta1), prob = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[5]) {
      fix <- dpois(as.integer(x), lambda = as.numeric(xTheta[[i]]$theta1), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[6]) {
      fix <- ddirac(as.numeric(x), location = as.numeric(xTheta[[i]]$theta1))
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[7]) {
      fix <- dgamma(as.numeric(x), scale = as.numeric(xTheta[[i]]$theta1), shape = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[8]) {
      fix <- dunif(as.numeric(x), min = as.numeric(xTheta[[i]]$theta1), max = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else    
    if (xTheta[[i]]$pdf == .rebmix$pdf[9]) {
      output <- .C(C_RvonMisesPdf,
        n = as.integer(n),
        y = as.double(x),
        Mean = as.double(xTheta[[i]]$theta1),
        Kappa = as.double(xTheta[[i]]$theta2),
        f = double(n),
        PACKAGE = "rebmix")

      fix <- output$f
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[10]) {
      output <- .C(C_RGumbelPdf,
        n = as.integer(n),
        y = as.double(x),
        Mean = as.double(xTheta[[i]]$theta1),
        Sigma = as.double(xTheta[[i]]$theta2),
        Xi = as.double(xTheta[[i]]$theta3),
        f = double(n),
        PACKAGE = "rebmix")

      fix <- output$f      
    }      

    if (yTheta[[i]]$pdf == .rebmix$pdf[1]) {
      fiy <- dnorm(as.numeric(y), mean = as.numeric(yTheta[[i]]$theta1), sd = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[2]) {
      fiy <- dlnorm(as.numeric(y), meanlog = as.numeric(yTheta[[i]]$theta1), sdlog = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[3]) {
      fiy <- dweibull(as.numeric(y), scale = as.numeric(yTheta[[i]]$theta1), shape = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[4]) {
      fiy <- dbinom(as.integer(y), size = as.integer(yTheta[[i]]$theta1), prob = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[5]) {
      fiy <- dpois(as.integer(y), lambda = as.numeric(yTheta[[i]]$theta1), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[6]) {
      fiy <- ddirac(as.numeric(y), location = as.numeric(yTheta[[i]]$theta1))
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[7]) {
      fiy <- dgamma(as.numeric(y), scale = as.numeric(yTheta[[i]]$theta1), shape = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[8]) {
      fiy <- dunif(as.numeric(y), min = as.numeric(yTheta[[i]]$theta1), max = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else    
    if (yTheta[[i]]$pdf == .rebmix$pdf[9]) {
      output <- .C(C_RvonMisesPdf,
        n = as.integer(n),
        y = as.double(y),
        Mean = as.double(yTheta[[i]]$theta1),
        Kappa = as.double(yTheta[[i]]$theta2),
        f = double(n),
        PACKAGE = "rebmix")

      fiy <- output$f
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[10]) {
      output <- .C(C_RGumbelPdf,
        n = as.integer(n),
        y = as.double(y),
        Mean = as.double(yTheta[[i]]$theta1),
        Sigma = as.double(yTheta[[i]]$theta2),
        Xi = as.double(yTheta[[i]]$theta3),
        f = double(n),
        PACKAGE = "rebmix")

      fiy <- output$f      
    }      

    f <- f + w[i] * fix * fiy
  }

  rm(list = ls()[!(ls() %in% c("f"))])

  return(f)
} ## .dfmix.xy
