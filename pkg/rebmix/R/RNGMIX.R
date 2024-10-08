setMethod("RNGMIX",
          signature(model = "RNGMIX"),
function(model, ...)
{
  Names <- names(model@Theta)

  pdf <- unlist(model@Theta[grep("pdf", Names)])

  theta1 <- unlist(model@Theta[grep("theta1", Names)])

  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(model@Theta[grep("theta2", Names)])

  theta2[is.na(theta2)] <- 0
  
  theta3 <- unlist(model@Theta[grep("theta3", Names)])

  theta3[is.na(theta3)] <- 0  

  c <- length(model@n); d <- length(model@Variables)

  length(pdf) <- d

  xmin <- rep(+Inf, d)
  xmax <- rep(-Inf, d)

  IDum <- model@rseed

  for (i in 1:length(model@Dataset.name)) {
    message("Dataset = ", model@Dataset.name[i])

    flush.console()

    output <- .C(C_RRNGMIX,
      IDum = as.integer(IDum),
      d = as.integer(d),
      c = as.integer(c),
      N = as.integer(model@n),
      length.pdf = as.integer(d),
      length.Theta = as.integer(3),
      length.theta = as.integer(c(d, d, d)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2, theta3)),
      n = integer(1),
      Y = double(sum(model@n) * d),
      Z = integer(sum(model@n)),
      error = integer(9),
      PACKAGE = "rebmix")

    error <- error.to.string(output$error);
      
    if (error[1] != "") {
      stop(error[1], call. = FALSE); return(NA)
    }
    
    if (error[2] != "") {
      warning(error[2], call. = FALSE, immediate. = TRUE)
    }  
    
    if (error[3] != "") {
      warning(error[3], call. = FALSE, immediate. = TRUE)
    } 

    dim(output$Y) <- c(output$n, d)

    xmin <- as.numeric(apply(rbind(xmin, output$Y), 2, min))
    xmax <- as.numeric(apply(rbind(xmax, output$Y), 2, max))

    model@Dataset[[i]] <- as.data.frame(output$Y, stringsAsFactors = FALSE)
    
    colnames(model@Dataset[[i]]) <- paste(1:d, sep = "")

    if (i == 1) {
      model@Zt <- factor(output$Z)
    }

    IDum <- IDum - 1
  }

  names(model@Dataset) <- model@Dataset.name

  model@w <- model@n / output$n

  model@ymin <- xmin
  model@ymax <- xmax

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RNGMIX

setMethod("RNGMIX",
          signature(model = "RNGMVNORM"),
function(model, ...)
{
  Names <- names(model@Theta)

  pdf <- unlist(model@Theta[grep("pdf", Names)])

  theta1 <- unlist(model@Theta[grep("theta1", Names)])

  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(model@Theta[grep("theta2", Names)])

  theta2[is.na(theta2)] <- 0

  c <- length(model@n); d <- length(model@Variables)

  length(pdf) <- 1

  xmin <- rep(+Inf, d)
  xmax <- rep(-Inf, d)

  IDum <- model@rseed

  for (i in 1:length(model@Dataset.name)) {
    message("Dataset = ", model@Dataset.name[i])

    flush.console()

    output <- .C(C_RRNGMVNORM,
      IDum = as.integer(IDum),
      d = as.integer(d),
      c = as.integer(c),
      N = as.integer(model@n),
      length.pdf = as.integer(d),
      length.Theta = as.integer(4),
      length.theta = as.integer(c(d, d * d, -d * d, -1)),
      Theta = as.double(c(theta1, theta2)),
      n = integer(1),
      Y = double(sum(model@n) * d),
      Z = integer(sum(model@n)),
      error = integer(9),
      PACKAGE = "rebmix")

    error <- error.to.string(output$error);
      
    if (error[1] != "") {
      stop(error[1], call. = FALSE); return(NA)
    }
    
    if (error[2] != "") {
      warning(error[2], call. = FALSE, immediate. = TRUE)
    }  
    
    if (error[3] != "") {
      warning(error[3], call. = FALSE, immediate. = TRUE)
    } 

    dim(output$Y) <- c(output$n, d)

    xmin <- as.numeric(apply(rbind(xmin, output$Y), 2, min))
    xmax <- as.numeric(apply(rbind(xmax, output$Y), 2, max))

    model@Dataset[[i]] <- as.data.frame(output$Y, stringsAsFactors = FALSE)
    
    colnames(model@Dataset[[i]]) <- paste(1:d, sep = "")

    if (i == 1) {
      model@Zt <- factor(output$Z)
    }

    IDum <- IDum - 1
  }

  names(model@Dataset) <- model@Dataset.name

  model@w <- model@n / output$n

  model@ymin <- xmin
  model@ymax <- xmax

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RNGMIX

setMethod("RNGMIX",
          signature(model = "ANY"),
function(model,
  Dataset.name,
  rseed,
  n,
  Theta)
{
  digits <- getOption("digits"); options(digits = 15)

  message("RNGMIX Version 2.16.1")

  flush.console()

  model <- new(model,
    Dataset.name = Dataset.name,
    rseed = rseed,
    n = n,
    Theta = Theta)

  output <- RNGMIX(model = model)

  model@Dataset <- output@Dataset
  model@Zt <- output@Zt
  model@w <- output@w
  model@ymin <- output@ymin
  model@ymax <- output@ymax

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## RNGMIX
