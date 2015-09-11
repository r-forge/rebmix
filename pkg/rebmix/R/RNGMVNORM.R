RNGMVNORM <- function(Dataset = NULL,
  rseed = -1,
  n = NULL,
  Theta = NULL)
{
  digits <- getOption("digits"); options(digits = 15)

  message("RNGMVNORM Version 2.7.3");
  
  flush.console()
  
  if (is.null(Dataset)) {
    stop(sQuote("Dataset"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(Dataset)) {
    stop(sQuote("Dataset"), " character vector is requested!", call. = FALSE)
  }

  if (!is.wholenumber(rseed)) {
    stop(sQuote("rseed"), " integer is requested!", call. = FALSE)
  }

  if (rseed > -1) {
    stop(sQuote("rseed"), " must be less than 0!", call. = FALSE)
  }

  if (is.null(n)) {
    stop(sQuote("n"), " must not be NULL!", call. = FALSE)
  }

  if (!is.wholenumber(n)) {
    stop(sQuote("n"), " integer or data frame is requested!", call. = FALSE)
  }
  
  if (!all(n > 0)) {
    stop("all ", sQuote("n"), " must be greater than 0!", call. = FALSE)
  }

  if (is.null(Theta)) {
    stop(sQuote("Theta"), " must not be NULL!", call. = FALSE)
  }

  if (!is.list(Theta)) {
    stop(sQuote("Theta"), " list is requested!", call. = FALSE)
  }
  
  Names <- names(Theta)
  
  c <- length(n)

  pdf <- as.character(unlist(Theta[grep("pdf", Names)]))
  
  pdf <- match.arg(pdf, .rebmix$pdf)
  
  theta1 <- as.numeric(unlist(Theta[grep("theta1", Names)]))
  
  theta1[is.na(theta1)] <- 0

  theta2 <- as.numeric(unlist(Theta[grep("theta2", Names)]))
  
  theta2[is.na(theta2)] <- 0

  d <- length(theta1) / c
  
  length(pdf) <- 1
  
  xmin <- rbind(rep(+Inf, d))
  xmax <- rbind(rep(-Inf, d))

  IDum <- rseed

  RNGMVNORM <- list()

  RNGMVNORM$Dataset <- list()
    
  for (i in 1:length(Dataset)) {
    message("Dataset = ", Dataset[i])

    flush.console()

    output <- .C("RRNGMVNORM",
      IDum = as.integer(IDum), 
      d = as.integer(d),
      c = as.integer(c),
      N = as.integer(n),
      length.pdf = as.integer(1),
      length.Theta = as.integer(3),
      length.theta = as.integer(c(d, length(theta2) / c, -length(theta2) / c)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = integer(1),
      Y = double(sum(n) * d),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RNGMVNORM!", call. = FALSE); return(NA)
    }

    dim(output$Y) <- c(output$n, d)

    xmin <- as.numeric(apply(rbind(xmin, output$Y), 2, min))
    xmax <- as.numeric(apply(rbind(xmax, output$Y), 2, max))

    RNGMVNORM$Dataset[[i]] <- as.data.frame(output$Y, stringsAsFactors = FALSE)

    IDum <- IDum - 1
  }
  
  names(RNGMVNORM$Dataset) <- Dataset
  
  RNGMVNORM$w <- n / output$n
  
  RNGMVNORM$Variables <- rep(.rebmix$Variables[1], d)

  RNGMVNORM$ymin <- xmin
  RNGMVNORM$ymax <- xmax
    
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("RNGMVNORM"))])

  class(RNGMVNORM) <- "RNGMVNORM"
  
  return(RNGMVNORM)
} ## RNGMVNORM
