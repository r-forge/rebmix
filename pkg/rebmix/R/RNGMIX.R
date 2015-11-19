setMethod("RNGMIX",
          signature(model = "missing"),
function(model,
  Dataset.name,
  rseed,
  n,
  Theta)
{
  digits <- getOption("digits"); options(digits = 15)

  message("RNGMIX Version 2.7.3")
  
  flush.console()
  
  if (is.null(Dataset.name)) {
    stop(sQuote("Dataset.name"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(Dataset.name)) {
    stop(sQuote("Dataset.name"), " character vector is requested!", call. = FALSE)
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
    stop(sQuote("n"), " integer is requested!", call. = FALSE)
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
  
  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)
  
  theta1 <- as.numeric(unlist(Theta[grep("theta1", Names)]))
  
  theta1[is.na(theta1)] <- 0

  theta2 <- as.numeric(unlist(Theta[grep("theta2", Names)]))
  
  theta2[is.na(theta2)] <- 0

  d <- length(theta1) / c
  
  length(pdf) <- d
  
  xmin <- rep(+Inf, d)
  xmax <- rep(-Inf, d)

  IDum <- rseed

  model <- new("RNGMIX")
  
  model@Dataset.name <- as.character(Dataset.name)
  model@rseed <- as.numeric(rseed)
  model@n <- as.numeric(n)
  model@Theta <- as.list(Theta)

  model@Dataset <- list()
    
  for (i in 1:length(Dataset.name)) {
    message("Dataset.name = ", Dataset.name[i])

    flush.console()

    output <- .C("RRNGMIX",
      IDum = as.integer(IDum), 
      d = as.integer(d),
      c = as.integer(c),
      N = as.integer(n),
      length.pdf = as.integer(d),
      length.Theta = as.integer(2),
      length.theta = as.integer(c(d, d)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = integer(1),
      Y = double(sum(n) * d),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RNGMIX!", call. = FALSE); return(NA)
    }

    dim(output$Y) <- c(output$n, d)

    xmin <- as.numeric(apply(rbind(xmin, output$Y), 2, min))
    xmax <- as.numeric(apply(rbind(xmax, output$Y), 2, max))

    model@Dataset[[i]] <- as.data.frame(output$Y, stringsAsFactors = FALSE)

    IDum <- IDum - 1
  }
  
  names(model@Dataset) <- Dataset.name
  
  model@w <- n / output$n
  
  M <- match(pdf, .rebmix$pdf)

  Index <- NULL

  for (i in 1:length(M)) {
    if (M[i] %in% which(.rebmix$pdf.nargs == 1)) {
      Index <- c(Index, i + 2)
    }
  }

  model@Variables <- rep(.rebmix$Variables[1], d)

  M <- na.omit(M)

  for (i in 1:length(M)) {
    if (M[i] %in% which(.rebmix$pdf.Variables == .rebmix$Variables[2])) {
      model@Variables[i] <- .rebmix$Variables[2]
    }
  }  

  model@ymin <- xmin
  model@ymax <- xmax
    
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("model"))])
  
  return(model)
}) ## RNGMIX

setMethod("RNGMIX",
          signature(model = "RNGMVNORM"),
function(model,
  Dataset.name,
  rseed,
  n,
  Theta)
{
  digits <- getOption("digits"); options(digits = 15)

  message("RNGMIX Version 2.7.3")
  
  flush.console()
  
  if (is.null(Dataset.name)) {
    stop(sQuote("Dataset.name"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(Dataset.name)) {
    stop(sQuote("Dataset.name"), " character vector is requested!", call. = FALSE)
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
    stop(sQuote("n"), " integer is requested!", call. = FALSE)
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
  
  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)
  
  theta1 <- as.numeric(unlist(Theta[grep("theta1", Names)]))
  
  theta1[is.na(theta1)] <- 0

  theta2 <- as.numeric(unlist(Theta[grep("theta2", Names)]))
  
  theta2[is.na(theta2)] <- 0

  d <- length(theta1) / c
  
  length(pdf) <- 1
  
  xmin <- rep(+Inf, d)
  xmax <- rep(-Inf, d)

  IDum <- rseed

  model@Dataset.name <- as.character(Dataset.name)
  model@rseed <- as.numeric(rseed)
  model@n <- as.numeric(n)
  model@Theta <- as.list(Theta)

  model@Dataset <- list()
    
  for (i in 1:length(Dataset.name)) {
    message("Dataset = ", Dataset.name[i])

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
      stop("in RNGMIX!", call. = FALSE); return(NA)
    }

    dim(output$Y) <- c(output$n, d)

    xmin <- as.numeric(apply(rbind(xmin, output$Y), 2, min))
    xmax <- as.numeric(apply(rbind(xmax, output$Y), 2, max))

    model@Dataset[[i]] <- as.data.frame(output$Y, stringsAsFactors = FALSE)

    IDum <- IDum - 1
  }
  
  names(model@Dataset) <- Dataset.name
  
  model@w <- n / output$n
  
  model@Variables <- rep(.rebmix$Variables[1], d)

  model@ymin <- xmin
  model@ymax <- xmax
    
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("model"))])
  
  return(model)
}) ## RNGMIX
