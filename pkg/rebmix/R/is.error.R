is.error <- function(Zt, # A vector of true cluster membership.
                     Zp) # A vector of predictive cluster membership.
{
  if (length(Zt) == 0) {
    error <- array(data = 0, dim = length(Zp))
  }
  else {
    error <- array(data = 1, dim = length(Zp))
  
    zt <- as.numeric(names(table(Zt))); zp <- unique(Zp)

    iset <- 1:length(zt)

    while (1) {
      imax <- 0; jmax <- 0; which.not.error <- NULL; Pmax <- 0.0

      for (i in iset) {
        Zti <- Zt[which(Zt == zt[i])]
        Zpi <- Zp[which(Zt == zt[i])]

        j <- as.numeric(names(sort(table(Zpi), decreasing = TRUE)))
      
        k <- match(j, zp); k <- k[!is.na(k)]

        if (length(k) > 0) {
          j <- zp[k[1]]; which.not.error <- which(Zt == zt[i] & Zp == j)

          P <- length(which.not.error) / length(Zti)

          if (P > Pmax) {
            imax <- i; jmax <- j; which.not.errormax <- which.not.error; Pmax <- P 
          }
        }
      }

      if (imax == 0) {
        break
      }
      else { 
        error[which.not.errormax] <- 0

        zp <- zp[which(zp != jmax)]; iset <- iset[which(iset != imax)]
      }
    }
  }
  
  error
} ## is.error

is.error.old <- function(Zt, # A vector of true cluster membership.
                         Zp) # A vector of predictive cluster membership.
{
  error <- array(data = 0, dim = length(Zp))

  if (length(Zt) != 0) {
    zt <- as.numeric(names(sort(table(Zt), decreasing = TRUE)))
    zp <- sort(unique(Zp))
    
    for (i in 1:length(zt)) {
      Zti <- Zt[which(Zt == zt[i])]
      Zpi <- Zp[which(Zt == zt[i])]

      j <- as.numeric(names(sort(table(Zpi), decreasing = TRUE)))
      
      k <- match(j, zp); k <- k[!is.na(k)]

      if (length(k) == 0) {
        error[which(Zt == zt[i])] <- 1
      }
      else {
        j <- zp[k[1]] 
      
        error[which(Zt == zt[i] & Zp != j)] <- 1

        zp <- zp[which(zp != j)] 
      }
    }
  }
  
  error
} ## is.error.old
