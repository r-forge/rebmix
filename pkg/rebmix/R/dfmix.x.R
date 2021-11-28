.dfmix.x <- function(x, w, xTheta, ...)
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

    f <- f + w[i] * fix
  }

  rm(list = ls()[!(ls() %in% c("f"))])

  return(f)
} ## .dfmix.x
