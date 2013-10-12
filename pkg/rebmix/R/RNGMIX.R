RNGMIX <- function(Dataset = NULL,
  rseed = -1,
  n = NULL,
  Theta = NULL)
{
  digits <- getOption("digits"); options(digits = 15)

  message("RNGMIX Version 2.5.1");
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

  if ((!is.wholenumber(n)) && (!is.data.frame(n))) {
    stop(sQuote("n"), " integer or data frame is requested!", call. = FALSE)
  }
  
  if (!all(n > 0)) {
    stop("all ", sQuote("n"), " must be greater than 0!", call. = FALSE)
  }

  if (is.null(Theta)) {
    stop(sQuote("Theta"), " must not be NULL!", call. = FALSE)
  }

  if ((!is.character(Theta)) && (!is.data.frame(Theta))) {
    stop(sQuote("Theta"), " character matrix or data frame is requested!", call. = FALSE)
  }
  
  if (length(n) != ncol(Theta)) {
    stop("number of columns in ", sQuote("n"), " and ", sQuote("Theta"), " must match!", call. = FALSE)
  }

  nrow <- nrow(Theta)
  ncol <- ncol(Theta)

  c <- ncol

  n <- rbind(n)

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

  xmin <- rbind(rep(+Inf, d))
  xmax <- rbind(rep(-Inf, d))

  IDum <- rseed

  RNGMIX <- list()

  RNGMIX$Dataset <- list()
    
  for (i in 1:length(Dataset)) {
    message("Dataset = ", Dataset[i])

    flush.console()

    output <- .C("RRNGMIX",
      IDum = as.integer(IDum), 
      d = as.integer(d),
      c = as.integer(c),
      N = as.integer(n),
      ParFamType = as.character(pdf),
      Par0 = as.double(theta1),
      Par1 = as.double(theta2),
      n = integer(1),
      X = double(sum(n) * d),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RNGMIX!", call. = FALSE); return(NA)
    }

    dim(output$X) <- c(output$n, d)

    xmin <- as.numeric(apply(rbind(xmin, output$X), 2, min))
    xmax <- as.numeric(apply(rbind(xmax, output$X), 2, max))

    RNGMIX$Dataset[[i]] <- as.data.frame(output$X, stringsAsFactors = FALSE)

    IDum <- IDum - 1
  }
  
  names(RNGMIX$Dataset) <- Dataset
  
  RNGMIX$w <- as.data.frame(rbind(n / output$n), stringsAsFactors = FALSE)

  rownames(RNGMIX$w) <- "w"
  colnames(RNGMIX$w) <- paste("comp", if (c > 1) 1:c else "", sep = "")

  RNGMIX$Theta <- rbind(pdf, theta1, theta2)

  dim(RNGMIX$Theta) <- c(3 * d, c)

  if (d > 1) {
    rownames(RNGMIX$Theta) <- c(paste("pdf", 1:d, sep = ""),
      paste("theta1.", 1:d, sep = ""),
      paste("theta2.", 1:d, sep = ""))
  }
  else {
    rownames(RNGMIX$Theta) <- c("pdf", "theta1", "theta2")
  }

  Index <- NULL

  for (i in 1:d){
    Index <- c(Index, seq(from = i, to = i + 2 * d, by = d))
  }

  RNGMIX$Theta <- cbind(RNGMIX$Theta[Index, ])

  M <- match(RNGMIX$Theta[, 1], .rebmix$pdf)

  Index <- NULL

  for (i in 1:length(M)) {
    if (M[i] %in% which(.rebmix$pdf.nargs == 1)) {
      Index <- c(Index, i + 2)
    }
  }

  if (is.null(Index)) {
    RNGMIX$Theta <- as.data.frame(RNGMIX$Theta, stringsAsFactors = FALSE)
  }
  else {
    RNGMIX$Theta <- as.data.frame(RNGMIX$Theta[-Index, ], stringsAsFactors = FALSE)
  }

  colnames(RNGMIX$Theta) <- paste("comp", if (c > 1) 1:c else "", sep = "")

  RNGMIX$Variables <- rep(.rebmix$Variables[1], d)

  M <- na.omit(M)

  for (i in 1:length(M)) {
    if (M[i] %in% which(.rebmix$pdf.Variables == .rebmix$Variables[2])) {
      RNGMIX$Variables[i] <- .rebmix$Variables[2]
    }
  }

  RNGMIX$ymin <- xmin
  RNGMIX$ymax <- xmax
  
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("RNGMIX"))])

  class(RNGMIX) <- "RNGMIX"
  
  return(RNGMIX)
} ## RNGMIX
