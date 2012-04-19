.extractTheta <- function(Theta)
{
  nrow <- nrow(Theta)
  ncol <- ncol(Theta)

  output <- array(data = list(NULL), dim = c(nrow, ncol), dimnames = NULL)

  for (i in 1:ncol) {
    M <- match(Theta[, i], .rebmix$pdf)

    j <- 1;

    for (k in 1:length(M)) {
      if (M[k] %in% which(.rebmix$pdf.nargs == 2)) {
        output[[j, i]]$pdf <- Theta[k, i]
        output[[j, i]]$theta1 <- as.numeric(Theta[k + 1, i])
        output[[j, i]]$theta2 <- as.numeric(Theta[k + 2, i])

        j <- j + 1
      }
      else
      if (M[k] %in% which(.rebmix$pdf.nargs == 1)) {
        output[[j, i]]$pdf <- Theta[k, i]
        output[[j, i]]$theta1 <- as.numeric(Theta[k + 1, i])

        j <- j + 1
      }
    }
  }

  j <- j - 1

  output <- output[1:j, ]; dim(output) <- c(j, ncol)

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .extractTheta
