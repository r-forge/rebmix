setMethod("RCLRMIX",
          signature(model = "RCLRMIX"),
function(model, ...)
{
  Names <- names(model@x@Theta[[model@pos]])

  pdf <- unlist(model@x@Theta[[model@pos]][grep("pdf", Names)])

  theta1 <- unlist(model@x@Theta[[model@pos]][grep("theta1", Names)])

  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(model@x@Theta[[model@pos]][grep("theta2", Names)])

  theta2[is.na(theta2)] <- 0
  
  theta3 <- unlist(model@x@Theta[[model@pos]][grep("theta3", Names)])

  theta3[is.na(theta3)] <- 0  

  c <- length(model@x@w[[model@pos]])

  w <- model@x@w[[model@pos]]

  d <- length(pdf) / c

  dataset <- model@Dataset

  if (missing(dataset) || (length(dataset) == 0)) {
    dataset <- as.matrix(model@x@Dataset[[model@pos]])
  }
  else {
    dataset <- as.matrix(dataset)
  }  

  n <- nrow(dataset)

  if (sum(pdf %in% .rebmix$pdf[c(4, 6)]) == c * d) {
    model@p <- w

    model@pi <- list(); nlevels <- array()

    for (i in 1:d) {
      for (j in 1:c) {
        if (pdf[(j - 1) * d + i] == .rebmix$pdf[4]) {
          if (j == 1) {
            nlevels[i] <- as.integer(theta1[(j - 1) * d + i]) + 1

            model@pi[[i]] <- matrix(data = 0.0, nrow = c, ncol = nlevels[i])

            colnames(model@pi[[i]]) <- paste(0:(nlevels[i] - 1), sep = "")
            rownames(model@pi[[i]]) <- paste(1:c, sep = "")
          }

          for (ii in 1:nlevels[i]) {
            model@pi[[i]][j, ii] <- dbinom(ii - 1, size = as.integer(theta1[(j - 1) * d + i]), prob = as.numeric(theta2[(j - 1) * d + i]))
          }
        }
        else
        if (pdf[(j - 1) * d + i] == .rebmix$pdf[6]) {
          if (j == 1) {
            nlevels[i] <- length(unique(dataset[, i]))

            model@pi[[i]] <- matrix(data = 0.0, nrow = c, ncol = nlevels[i])

            colnames(model@pi[[i]]) <- paste(0:(nlevels[i] - 1), sep = "")
            rownames(model@pi[[i]]) <- paste(1:c, sep = "")
          }

          for (ii in 1:nlevels[i]) {
            model@pi[[i]][j, ii] <- ddirac(ii - 1, location = as.integer(theta1[(j - 1) * d + i]))
          }
        }
      }
    }

    Y <- dataset; y <- as.matrix(unique(Y)); Nt <- array(); Np <- array()

    for (j in 1:nrow(y)) {
      x <- array(); k <- 1

      for (l in 1:nrow(Y)) {
        if (all(y[j, ] == Y[l, ])) {
          x[k] <- l; k <- k + 1
        }
      }

      Nt[j] <- length(x); Y <- as.matrix(Y[-x, ]); Np[j] <- 0.0

      for (l in 1:c) {
        Pl <- 1.0

        for(i in 1:d) {
          for (ii in 1:length(model@pi[[i]][l, ])) {
            if (y[j, i] == ii - 1) {
              Pl <- Pl * model@pi[[i]][l, ii]
            }
          }
        }

        Np[j] <- Np[j] + model@p[l] * Pl * n
      }
    }

    model@P <- as.data.frame(cbind(y, Nt, Np))

    if (is.null(colnames(dataset))) {
      colnames(model@P) <- paste(c(1:d, "Nt", "Np"), sep = "")
      
      names(model@pi) <- 1:d
    }
    else {
      colnames(model@P) <- c(colnames(dataset), "Nt", "Np")
      
      names(model@pi) <- colnames(dataset)
    }

    rownames(model@pi[[i]]) <- paste(1:c, sep = "")
  }

  output <- .C(C_RCombineComponentsMIX,
    c = as.integer(c),
    w = as.double(model@x@w[[model@pos]]),
    length.pdf = as.integer(d),
    length.Theta = as.integer(3),
    length.theta = as.integer(c(d, d, d)),
    pdf = as.character(pdf),
    Theta = as.double(c(theta1, theta2, theta3)),
    n = as.integer(n),
    x = as.double(dataset),
    tau = double(n * c),
    F = integer(c),
    T = integer(c),
    EN = double(c),
    ED = double(c),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RCLRMIX!", call. = FALSE); return(NA)
  }

  model@tau <- matrix(data = output$tau, ncol = c, byrow = TRUE)

  colnames(model@tau) <- paste(1:c, sep = "")
  rownames(model@tau) <- paste(1:n, sep = "")

  if (output$c > 1) {
    model@from <- output$F
    model@to <- output$T
    model@EN <- output$EN
    model@ED <- output$ED
  }

  output <- .C(C_RCLRMIX,
    n = n,
    X = as.double(dataset),
    d = as.integer(d),
    c = as.integer(unlist(c)),
    w = as.double(unlist(w)),
    pdf = as.character(unlist(pdf)),
    theta1 = as.double(unlist(theta1)),
    theta2 = as.double(unlist(theta2)),
    theta3 = as.double(unlist(theta3)),
    Z = integer(n),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RCLRMIX!", call. = FALSE); return(NA)
  }

  unique.Z <- unique(output$Z)

  model@c <- length(unique.Z)

  i <- length(model@from)

  while (i > 0) {
    from.in.unique.Z <- model@from[i] %in% unique.Z
    to.in.unique.Z <- model@to[i] %in% unique.Z

    if (from.in.unique.Z && to.in.unique.Z) {
    }
    else
    if (from.in.unique.Z) {
      j <- 1

      while (j < i) {
        if (model@from[j] == model@to[i]) model@from[j] <- model@from[i]
        if (model@to[j] == model@to[i]) model@to[j] <- model@from[i]

        j <- j + 1
      }

      model@from <- model@from[-i]
      model@to <- model@to[-i]
      model@EN <- model@EN[-i]
      model@ED <- model@ED[-i]
    }
    else {
      model@from <- model@from[-i]
      model@to <- model@to[-i]
      model@EN <- model@EN[-i]
      model@ED <- model@ED[-i]
    }

    i <- i - 1
  }

  model@Zp <- as.factor(output$Z)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RCLRMIX

setMethod("RCLRMIX",
          signature(model = "RCLRMVNORM"),
function(model, ...)
{
  Names <- names(model@x@Theta[[model@pos]])

  pdf <- unlist(model@x@Theta[[model@pos]][grep("pdf", Names)])

  theta1 <- unlist(model@x@Theta[[model@pos]][grep("theta1", Names)])

  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(model@x@Theta[[model@pos]][grep("theta2", Names)])

  theta2[is.na(theta2)] <- 0

  c <- length(model@x@w[[model@pos]])

  w <- model@x@w[[model@pos]]

  d <- length(pdf) / c
  
  dataset <- model@Dataset

  if (missing(dataset) || (length(dataset) == 0)) {
    dataset <- as.matrix(model@x@Dataset[[model@pos]])
  }
  else {
    dataset <- as.matrix(dataset)
  } 

  n <- nrow(dataset)

  output <- .C(C_RCombineComponentsMVNORM,
    c = as.integer(c),
    w = as.double(model@x@w[[model@pos]]),
    length.pdf = as.integer(d),
    length.Theta = as.integer(4),
    length.theta = as.integer(c(d, d * d, -d * d, -1)),
    pdf = as.character(pdf),
    Theta = as.double(c(theta1, theta2)),
    n = as.integer(n),
    x = as.double(dataset),
    tau = double(n * c),
    F = integer(c),
    T = integer(c),
    EN = double(c),
    ED = double(c),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RCLRMIX!", call. = FALSE); return(NA)
  }

  model@tau <- matrix(data = output$tau, ncol = c, byrow = TRUE)

  colnames(model@tau) <- paste(1:c, sep = "")
  rownames(model@tau) <- paste(1:n, sep = "")

  if (output$c > 1) {
    model@from <- output$F
    model@to <- output$T
    model@EN <- output$EN
    model@ED <- output$ED
  }

  output <- .C(C_RCLRMVNORM,
    n = n,
    X = as.double(dataset),
    d = as.integer(d),
    c = as.integer(unlist(c)),
    w = as.double(unlist(w)),
    pdf = as.character(unlist(pdf)),
    theta1 = as.double(unlist(theta1)),
    theta2 = as.double(unlist(theta2)),
    Z = integer(n),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in RCLRMIX!", call. = FALSE); return(NA)
  }

  unique.Z <- unique(output$Z)

  model@c <- length(unique.Z)

  i <- length(model@from)

  while (i > 0) {
    from.in.unique.Z <- model@from[i] %in% unique.Z
    to.in.unique.Z <- model@to[i] %in% unique.Z

    if (from.in.unique.Z && to.in.unique.Z) {
    }
    else
    if (from.in.unique.Z) {
      j <- 1

      while (j < i) {
        if (model@from[j] == model@to[i]) model@from[j] <- model@from[i]
        if (model@to[j] == model@to[i]) model@to[j] <- model@from[i]

        j <- j + 1
      }

      model@from <- model@from[-i]
      model@to <- model@to[-i]
      model@EN <- model@EN[-i]
      model@ED <- model@ED[-i]
    }
    else {
      model@from <- model@from[-i]
      model@to <- model@to[-i]
      model@EN <- model@EN[-i]
      model@ED <- model@ED[-i]
    }

    i <- i - 1
  }

  model@Zp <- as.factor(output$Z)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RCLRMIX

setMethod("RCLRMIX",
          signature(model = "ANY"),
function(model,
  x,
  Dataset,
  pos,
  Zt, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  message("RCLRMIX Version 2.14.0")

  flush.console()

  model <- new(model,
    x = x,
    Dataset = Dataset,
    pos = pos,
    Zt = Zt)

  model <- RCLRMIX(model = model, ...)

  Zp <- as.numeric(levels(model@Zp))[model@Zp]
  Zt <- as.numeric(levels(model@Zt))[model@Zt]

  if (length(Zt) > 0) {
    prob <- array(data = 0.0, dim = model@c)

    for (i in model@c:1) {
      if (i < model@c) {
        Zp[Zp == model@from[i]] <- model@to[i]
      }

      error <- is.error(Zt, Zp)

      prob[i] <- length(error[error == 0]) / length(error)
    }

    model@prob <- as.numeric(prob)
  }

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RCLRMIX
