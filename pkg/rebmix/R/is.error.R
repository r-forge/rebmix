is.error <- function(Zt, # A vector of true cluster membership.
                     Zp) # A vector of predictive cluster membership.
{
  error <- array(data = 0, dim = length(Zp))

  if (length(Zt) != 0) {
    zt <- sort(table(Zt), decreasing = TRUE)
    zp <- sort(unique(Zp))  

    for (i in 1:length(zt)) {
      Zti <- Zt[which(Zt == zt[i])]
      Zpi <- Zp[which(Zt == zt[i])]

      j <- as.numeric(names(sort(table(Zpi), decreasing = TRUE)))

      k <- which(j %in% zp)

      if (length(k) == 0) {
        error[which(Zt == zt[i])] = 1
      }
      else {
        j <- j[k[1]]
    
        error[which(Zt == zt[i] & Zp != j)] = 1

        zp <- zp[which(zp != j)] 
      }
    }
  }
  
  error
} ## is.error
