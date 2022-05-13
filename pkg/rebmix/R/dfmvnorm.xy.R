.dfmvnorm.xy <- function(x, y, w, Theta, xi, yi, ...)
{
  n <- length(x)
  
  f <- array(data = 0.0, dim = n, dimnames = NULL)
  
  for (i in 1:length(w)) {
    output <- .C(C_RTvtNormalPdf,
      n = as.integer(n),
      x = as.double(x),
      y = as.double(y),
      Mean = as.double(Theta[[i]]$theta1[c(xi, yi)]),
      Sigma = as.double(Theta[[i]]$theta2[c(xi, yi), c(xi, yi)]),
      f = double(n),
      PACKAGE = "rebmix")

    fi <- output$f   

    f <- f + w[i] * fi
  }

  rm(list = ls()[!(ls() %in% c("f"))])
  
  return(f)
} ## .dfmvnorm.xy
