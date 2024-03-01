### Panic Branislav.
setMethod("EMMIX",
          signature(model = "REBMIX"),
function(model, Theta, ...)
{
  summary <- list()

  summary.EM <- list()

  w <- Theta@w

  pdf <- Theta@pdf

  d <- Theta@d

  c <- Theta@c

  Names = names(Theta@Theta)

  theta1 <- unlist(Theta@Theta[grep("theta1", Names)])

  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(Theta@Theta[grep("theta2", Names)])

  theta2[is.na(theta2)] <- 0

  theta3 <- unlist(Theta@Theta[grep("theta3", Names)])

  theta3[is.na(theta3)] <- 0

  for (i in 1:length(model@Dataset)) {
    Dataset.name <- names(model@Dataset)[i]
    
    Dataset <- model@Dataset[[i]]
    
    if (as.character(class(Dataset)) == "data.frame") {
      Y.type <- 0
    
      X <- as.matrix(model@Dataset[[i]])

      n <- nrow(X)
      d <- ncol(X)    
    }
    else
    if (as.character(class(Dataset)) == "Histogram") {
      Y.type <- 1
    
      X <- as.matrix(model@Dataset[[i]]@Y)

      n <- nrow(X)
      d <- ncol(X) - 1 
    }     

    message("Dataset = ", Dataset.name)

    flush.console()

    output <- .C(C_REMMIX,
      d = as.integer(d),
      n = as.integer(n),
      Y = as.double(X),
      Y.type = as.integer(Y.type),
      pdf = as.character(model@pdf),
      c = as.integer(c),
      w = as.double(w),
      Theta1 = as.double(theta1),
      Theta2 = as.double(theta2),
      Theta3 = as.double(theta3),
      EMVariant = as.character(model@EMcontrol@variant),
      EMAcceleration = as.character(model@EMcontrol@acceleration),
      EMTolerance = as.double(model@EMcontrol@tolerance),
      EMAccelerationMul = as.double(model@EMcontrol@acceleration.multiplier),
      EMMaxIter = as.integer(model@EMcontrol@maximum.iterations),
      EMK = as.integer(model@EMcontrol@K),
      EMMerge = as.integer(model@EMcontrol@eliminate.zero.components),
      n_iter = integer(1),
      summary.logL = double(1),
      summary.M = integer(1),
      Error = character(1),
      PACKAGE = "rebmix")

    error <- unlist(strsplit(output$error, "\n"));
      
    if (error[1] != "") {
      stop(error[1], call. = FALSE); return(NA)
    }
   
    if (error[2] != "") {
      warning(error[2], call. = FALSE, immediate. = TRUE)
    }  
    
    if (error[3] != "") {
      warning(error[3], call. = FALSE, immediate. = TRUE)
    }  

    if (all(output$w == 0.0)) {
      warning("EM did not converge for ", Dataset.name, " dataset!")
    }

    c <- output$c

    length(output$w) <- c

    length(output$Theta1) <- c * d

    length(output$Theta2) <- c * d

    length(output$Theta3) <- c * d

    model@w[[i]] <- output$w

    model@Theta[[i]] <- list()

    length(model@Theta[[i]]) <- 4 * c

    names(model@Theta[[i]])[seq(1, 4 * c, 4)] <- paste("pdf", 1:c, sep = "")
    names(model@Theta[[i]])[seq(2, 4 * c, 4)] <- paste("theta1.", 1:c, sep = "")
    names(model@Theta[[i]])[seq(3, 4 * c, 4)] <- paste("theta2.", 1:c, sep = "")
    names(model@Theta[[i]])[seq(4, 4 * c, 4)] <- paste("theta3.", 1:c, sep = "")

    M1 <- which(model@pdf %in% .rebmix$pdf[.rebmix$pdf.nargs < 2])
    M2 <- which(model@pdf %in% .rebmix$pdf[.rebmix$pdf.nargs < 3])

    for (j in 1:c) {
	  model@Theta[[i]][[1 + (j - 1) * 4]] <- model@pdf
	  model@Theta[[i]][[2 + (j - 1) * 4]] <- output$Theta1[seq((j - 1) * d + 1, j * d, 1)]
  	  model@Theta[[i]][[3 + (j - 1) * 4]] <- output$Theta2[seq((j - 1) * d + 1, j * d, 1)]
  	  model@Theta[[i]][[4 + (j - 1) * 4]] <- output$Theta3[seq((j - 1) * d + 1, j * d, 1)]

      model@Theta[[i]][[3 + (j - 1) * 4]][M1] <- NA
      model@Theta[[i]][[4 + (j - 1) * 4]][M2] <- NA
    }

    summary[[i]] <- c(Dataset.name,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      output$c,
      NA,
      NA,
      rep(NA, d),
      rep(NA, d),
      rep(NA, d),
      rep(NA, d),
      NA,
      output$summary.logL,
      output$summary.M)

    summary.EM[[i]] <- c(Dataset.name,
      .rebmix$EMStrategy[4],
      output$EMVariant,
      output$EMAcceleration,
      output$EMAccelerationMul,
      output$EMTolerance,
      output$EMMaxIter,
      output$n_iter,
      output$n_iter)
    
    model@opt.c[[i]] <- NA
    model@opt.IC[[i]] <- NA
    model@opt.logL[[i]] <- NA
    model@opt.D[[i]] <- NA
    model@all.K[[i]] <- NA
    model@all.IC[[i]] <- NA
  }

  model@summary <- as.data.frame(do.call("rbind", summary), stringsAsFactors = FALSE)

  model@summary.EM <- as.data.frame(do.call("rbind", summary.EM), stringsAsFactors = FALSE)

  colnames(model@summary) <- c("Dataset",
    "Preprocessing",
    "cmax",
    "cmin",
    "Criterion",
    "ar",
    "Restraints",
    "c",
    "v/k",
    "K",
    paste("y0", if (d > 1) 1:d else "", sep = ""),
    paste("ymin", if (d > 1) 1:d else "", sep = ""),
    paste("ymax", if (d > 1) 1:d else "", sep = ""),
    paste("h", if (d > 1) 1:d else "", sep = ""),
    "IC",
    "logL",
    "M")
    
  colnames(model@summary.EM) <- c("Dataset",
    "strategy",
    "variant",
    "acceleration",
    "acceleration.multiplier",
    "tolerance",
    "maximum.iterations",
    "opt.iterations.nbr",
    "total.iterations.nbr")

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## EMMIX

setMethod("EMMIX",
          signature(model = "REBMVNORM"),
function(model, Theta, ...)
{
  summary <- list()

  summary.EM <- list()

  w <- Theta@w

  pdf <- Theta@pdf

  d <- Theta@d
  
  c <- Theta@c 

  Names <- names(Theta@Theta)

  theta1 <- unlist(Theta@Theta[grep("theta1", Names)])

  theta1[is.na(theta1)] <- 0

  theta2 <- unlist(Theta@Theta[grep("theta2", Names)])

  theta2[is.na(theta2)] <- 0 

  for (i in 1:length(model@Dataset)) {
    Dataset.name <- names(model@Dataset)[i]

    Dataset <- model@Dataset[[i]]
    
    if (as.character(class(Dataset)) == "data.frame") {
      Y.type <- 0
    
      X <- as.matrix(model@Dataset[[i]])

      n <- nrow(X)
      d <- ncol(X)    
    }
    else
    if (as.character(class(Dataset)) == "Histogram") {
      Y.type <- 1
    
      X <- as.matrix(model@Dataset[[i]]@Y)

      n <- nrow(X)
      d <- ncol(X) - 1 
    }  

    message("Dataset = ", Dataset.name)

    flush.console()

    output <- .C(C_REMMVNORM,
      d = as.integer(d),
      n = as.integer(n),
      Y = as.double(X),
      Y.type = as.integer(Y.type),
      pdf = as.character(model@pdf),
      c = as.integer(c),
      w = as.double(w),
      Theta1 = as.double(theta1),
      Theta2 = as.double(theta2),
      EMVariant = as.character(model@EMcontrol@variant),
      EMAcceleration = as.character(model@EMcontrol@acceleration),
      EMTolerance = as.double(model@EMcontrol@tolerance),
      EMAccelerationMul = as.double(model@EMcontrol@acceleration.multiplier),
      EMMaxIter = as.integer(model@EMcontrol@maximum.iterations),
      EMK = as.integer(model@EMcontrol@K),
      EMMerge = as.integer(model@EMcontrol@eliminate.zero.components),
      n_iter = integer(1),
      summary.logL = double(1),
      summary.M = integer(1),
      Error = character(1),
      PACKAGE = "rebmix")

    error <- unlist(strsplit(output$error, "\n"));
      
    if (error[1] != "") {
      stop(error[1], call. = FALSE); return(NA)
    }
   
    if (error[2] != "") {
      warning(error[2], call. = FALSE, immediate. = TRUE)
    }  
    
    if (error[3] != "") {
      warning(error[3], call. = FALSE, immediate. = TRUE)
    }

    c <- output$c

    length(output$w) <- c

    if (all(output$w == 0.0)) {
      warning("EM did not converge for ", Dataset.name, " dataset!")
    }

    length(output$Theta1) <- c * d

    length(output$Theta2) <- c * d * d

    model@w[[i]] <- output$w

    model@Theta[[i]] <- list()

    length(model@Theta[[i]]) <- 3 * c

    names(model@Theta[[i]])[seq(1, 3 * c, 3)] <- paste("pdf", 1:c, sep = "")
    names(model@Theta[[i]])[seq(2, 3 * c, 3)] <- paste("theta1.", 1:c, sep = "")
    names(model@Theta[[i]])[seq(3, 3 * c, 3)] <- paste("theta2.", 1:c, sep = "")

    M <- which(pdf %in% .rebmix$pdf[.rebmix$pdf.nargs == 1])

    for (j in 1:c) {
	  model@Theta[[i]][[1 + (j - 1) * 3]] <- model@pdf
	  model@Theta[[i]][[2 + (j - 1) * 3]] <- output$Theta1[seq((j - 1) * d + 1, j * d, 1)]
  	  model@Theta[[i]][[3 + (j - 1) * 3]] <- output$Theta2[seq((j - 1) * d * d + 1, j * d * d, 1)]

      model@Theta[[i]][[3 + (j - 1) * 3]][M] <- NA
    }

    summary[[i]] <- c(Dataset.name,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      output$c,
      NA,
      NA,
      rep(NA, d),
      rep(NA, d),
      rep(NA, d),
      rep(NA, d),
      NA,
      output$summary.logL,
      output$summary.M)

    summary.EM[[i]] <- c(Dataset.name,
      .rebmix$EMStrategy[4],
      output$EMVariant,
      output$EMAcceleration,
      output$EMAccelerationMul,
      output$EMTolerance,
      output$EMMaxIter,
      output$n_iter,
      output$n_iter)
    
    model@opt.c[[i]] <- NA
    model@opt.IC[[i]] <- NA
    model@opt.logL[[i]] <- NA
    model@opt.D[[i]] <- NA
    model@all.K[[i]] <- NA
    model@all.IC[[i]] <- NA
  }

  model@summary <- as.data.frame(do.call("rbind", summary), stringsAsFactors = FALSE)

  model@summary.EM <- as.data.frame(do.call("rbind", summary.EM), stringsAsFactors = FALSE)

  colnames(model@summary) <- c("Dataset",
    "Preprocessing",
    "cmax",
    "cmin",
    "Criterion",
    "ar",
    "Restraints",
    "c",
    "v/k",
    "K",
    paste("y0", if (d > 1) 1:d else "", sep = ""),
    paste("ymin", if (d > 1) 1:d else "", sep = ""),
    paste("ymax", if (d > 1) 1:d else "", sep = ""),
    paste("h", if (d > 1) 1:d else "", sep = ""),
    "IC",
    "logL",
    "M")

  colnames(model@summary.EM) <- c("Dataset",
    "strategy",
    "variant",
    "acceleration",
    "acceleration.multiplier",
    "tolerance",
    "maximum.iterations",
    "opt.iterations.nbr",
    "total.iterations.nbr")

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## EMMIX

setMethod("EMMIX",
          signature(model = "ANY"),
function(model,
  Dataset,
  Theta, 
  EMcontrol, ...)  
{
  digits <- getOption("digits"); options(digits = 15)

  Theta.model <- paste("EM", substr(model, 4, nchar(model)), ".Theta", sep = "")

  if (missing(Theta) || length(Theta) == 0) {
    stop(sQuote("Theta"), " must not be empty!", call. = FALSE)
  }

  if (as.character(class(Theta)) != Theta.model) {
    stop(sQuote("Theta"), " object of class ", Theta.model, " is requested!", call. = FALSE)
  }

  if (missing(EMcontrol) || length(EMcontrol) == 0) {
    EMcontrol <- new("EM.Control", strategy = "single")
  }
  else {
    EMcontrol@strategy <- "single"
  }

  model <- new(model,
     Dataset = Dataset,
     Preprocessing = "histogram",
     cmax = Theta@c,
     cmin = Theta@c,
     Criterion = "AIC",
     pdf = Theta@pdf,
     theta1 = numeric(),
     theta2 = numeric(),
     theta3 = rep(NA, Theta@d),
     K = "auto",
     y0 = numeric(),
     ymin = numeric(),
     ymax = numeric(),
     ar = 0.1,
     Restraints = "loose",
     EMcontrol = EMcontrol)

  output <- EMMIX(model = model, Theta = Theta, ...)

  for (k in (1:length(Dataset))) {
    model@w[[length(model@w) + 1]] <- output@w[[k]]
    model@Theta[[length(model@Theta) + 1]] <- output@Theta[[k]]
  }

  if (is.null(model@summary)) {
    model@summary <- output@summary
  }
  else {
    model@summary <- merge(model@summary, output@summary, all = TRUE, sort = FALSE)
  }
  
  if (is.null(model@summary.EM)) {
    model@summary.EM <- output@summary.EM 
  }
  else {
    model@summary.EM <- merge(model@summary.EM, output@summary.EM, all = TRUE, sort = FALSE)
  }

  for (k in (1:length(Dataset))) {
    model@opt.c[[length(model@opt.c) + 1]] <- output@opt.c[[k]]
    model@opt.IC[[length(model@opt.IC) + 1]] <- output@opt.IC[[k]]
    model@opt.logL[[length(model@opt.logL) + 1]] <- output@opt.logL[[k]]
    model@opt.D[[length(model@opt.D) + 1]] <- output@opt.D[[k]]
    model@all.K[[length(model@all.K) + 1]] <- output@all.K[[k]]
    model@all.IC[[length(model@all.IC) + 1]] <- output@all.IC[[k]]
  }

  model@pos <- which(as.numeric(model@summary[, "logL"]) == max(as.numeric(model@summary[, "logL"])))

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
})
### End
