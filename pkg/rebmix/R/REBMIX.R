setMethod("REBMIX",
          signature(model = "REBMIX"),
function(model, ...)
{
  summary <- list()

### Panic Branislav.

  summary.EM <- list()

### End  

  for (i in 1:length(model@Dataset)) {
    Dataset.name <- names(model@Dataset)[i]

    X <- as.matrix(model@Dataset[[i]])

    message("Dataset = ", Dataset.name)

    flush.console()

    n <- nrow(X)
    d <- ncol(X)
    
    length.pdf <- d

    if (is.character(model@K)) {
      if (model@Preprocessing == .rebmix$Preprocessing[3]) {
        K <- as.integer(sqrt(n) * 0.75):as.integer(sqrt(n) * 1.25)
      }
      else {
        Sturges <- ceiling(1.0 + log2(n)); RootN <- ceiling(2.0 * n^0.5)

        K <- kseq(Sturges, RootN, f = 0.25)
      }
      
      K <- unique(K) 
      
      dK <- max(K) - min(K) + 1; length.K <- length(K); K <- rep(K, d)
    }
    else {
      if (is.matrix(model@K)) {
        K <- model@K[i, ]
        
        dK <- 1; length.K <- 1
      }
      else {
        K <- sort(unique(model@K))
        
        dK <- max(K) - min(K) + 1; length.K <- length(K); K <- rep(K, d)
      }
    }

    if (length(model@theta1) > 0) {
      length.theta1 <- +d; theta1 <- model@theta1; theta1[is.na(theta1)] <- 0
    }
    else {
      length.theta1 <- -d; theta1 <- numeric()
    }

    if (length(model@theta2) > 0) {
      length.theta2 <- +d; theta2 <- model@theta2; theta2[is.na(theta2)] <- 0
    }
    else {
      length.theta2 <- -d; theta2 <- numeric()
    }

    output <- .C(C_RREBMIX,
      Preprocessing = as.character(model@Preprocessing),
      cmax = as.integer(model@cmax),
      cmin = as.integer(model@cmin),
      Criterion = as.character(model@Criterion),
      d = as.integer(d),
      Variables = as.character(model@Variables),
      length.pdf = as.integer(length.pdf),
      pdf = as.character(model@pdf),
      length.Theta = as.integer(2),
      length.theta = as.integer(c(length.theta1, length.theta2)),
      Theta = as.double(c(theta1, theta2)),
      length.K = as.integer(length.K),
      K = as.integer(K),
      length.y0 = as.integer(length(model@y0)),
      y0 = as.double(model@y0),
      length.ymin = as.integer(length(model@ymin)),
      ymin = as.double(model@ymin),
      length.ymax = as.integer(length(model@ymax)),
      ymax = as.double(model@ymax),
      ar = as.double(model@ar),
      Restraints = as.character(model@Restraints),
      n = as.integer(n),
      Y = as.double(X),
### Panic Branislav.      
      EMStrategy = as.character(model@EMcontrol@strategy),
      EMVariant = as.character(model@EMcontrol@variant),
      EMAcceleration = as.character(model@EMcontrol@acceleration),
      EMTolerance = as.double(model@EMcontrol@tolerance),
      EMAccelerationMul = as.double(model@EMcontrol@acceleration.multiplier),
      EMMaxIter = as.integer(model@EMcontrol@maximum.iterations),
      n_iter = integer(1),
      n_iter_all = integer(1),
### End            
      summary.k = integer(1),
      summary.h = double(d),
      summary.y0 = double(d),
      summary.ymin = double(d),
      summary.ymax = double(d),
      summary.IC = double(1),
      summary.logL = double(1),
      summary.M = integer(1),
      summary.c = integer(1),
      w = double(model@cmax),
      theta1 = double(model@cmax * d),
      theta2 = double(model@cmax * d),
      opt.length = integer(1),
      opt.c = integer(1000), # 1000 = ItMax see rebmixf.h
      opt.IC = double(1000),
      opt.logL = double(1000),
      opt.D = double(1000),
      all.length = integer(1),
      all.K = integer(dK),
      all.IC = double(dK),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RREBMIX!", call. = FALSE); return(NA) #
    }

    c <- output$summary.c

    length(output$summary.h) <- d
    length(output$w) <- c
    length(output$theta1) <- c * d
    length(output$theta2) <- c * d

    length(output$opt.c) <- output$opt.length
    length(output$opt.IC) <- output$opt.length
    length(output$opt.logL) <- output$opt.length
    length(output$opt.D) <- output$opt.length

    j <- order(output$opt.c, output$opt.logL)

    output$opt.c <- output$opt.c[j]
    output$opt.IC <- output$opt.IC[j]
    output$opt.logL <- output$opt.logL[j]
    output$opt.D <- output$opt.D[j]

    j <- !duplicated(output$opt.c, fromLast = TRUE)

    output$opt.c <- output$opt.c[j]
    output$opt.IC <- output$opt.IC[j]
    output$opt.logL <- output$opt.logL[j]
    output$opt.D <- output$opt.D[j]

    length(output$all.K) <- output$all.length
    length(output$all.IC) <- output$all.length

    model@w[[i]] <- output$w

    model@Theta[[i]] <- list()

    length(model@Theta[[i]]) <- 3 * c

    names(model@Theta[[i]])[seq(1, 3 * c, 3)] <- paste("pdf", 1:c, sep = "")
    names(model@Theta[[i]])[seq(2, 3 * c, 3)] <- paste("theta1.", 1:c, sep = "")
    names(model@Theta[[i]])[seq(3, 3 * c, 3)] <- paste("theta2.", 1:c, sep = "")

    M <- which(pdf %in% .rebmix$pdf[.rebmix$pdf.nargs == 1])

    for (j in 1:c) {
	  model@Theta[[i]][[1 + (j - 1) * 3]] <- model@pdf
	  model@Theta[[i]][[2 + (j - 1) * 3]] <- output$theta1[seq((j - 1) * d + 1, j * d, 1)]
  	  model@Theta[[i]][[3 + (j - 1) * 3]] <- output$theta2[seq((j - 1) * d + 1, j * d, 1)]

      model@Theta[[i]][[3 + (j - 1) * 3]][M] <- NA
    }

    output$K <- paste("c(", paste(K, collapse = ","), ")", sep = "")

    if (model@Preprocessing == .rebmix$Preprocessing[1]) {
      length(output$summary.y0) <- d

      summary[[i]] <- c(Dataset.name,
        output$Preprocessing,
        output$cmax,
        output$cmin,
        output$Criterion,
        output$ar,
        output$Restraints,
        output$summary.c,
        output$summary.k,
        output$K,
        output$summary.y0,
        output$summary.ymin,
        output$summary.ymax,
        output$summary.h,
        output$summary.IC,
        output$summary.logL,
        output$summary.M)
    }
    else
    if (model@Preprocessing == .rebmix$Preprocessing[2]) {
      summary[[i]] <- c(Dataset.name,
        output$Preprocessing,
        output$cmax,
        output$cmin,
        output$Criterion,
        output$ar,
        output$Restraints,
        output$summary.c,
        output$summary.k,
        output$K,
        rep(NA, d),
        output$summary.ymin,
        output$summary.ymax,
        output$summary.h,
        output$summary.IC,
        output$summary.logL,
        output$summary.M)
    }
    if (model@Preprocessing == .rebmix$Preprocessing[3]) {
      summary[[i]] <- c(Dataset.name,
        output$Preprocessing,
        output$cmax,
        output$cmin,
        output$Criterion,
        output$ar,
        output$Restraints,
        output$summary.c,
        output$summary.k,
        output$K,
        rep(NA, d),
        output$summary.ymin,
        output$summary.ymax,
        output$summary.h,
        output$summary.IC,
        output$summary.logL,
        output$summary.M)
    }
    
### Panic Branislav.

    summary.EM[[i]] <- c(Dataset.name,
      output$EMStrategy,
      output$EMVariant,
      output$EMAcceleration,
      output$EMAccelerationMul,
      output$EMTolerance,
      output$EMMaxIter,
      output$n_iter,
      output$n_iter_all)

### End     

    model@opt.c[[i]] <- output$opt.c
    model@opt.IC[[i]] <- output$opt.IC
    model@opt.logL[[i]] <- output$opt.logL
    model@opt.D[[i]] <- output$opt.D
    model@all.K[[i]] <- output$all.K
    model@all.IC[[i]] <- output$all.IC
  }

  model@summary <- as.data.frame(do.call("rbind", summary), stringsAsFactors = FALSE)
  
### Panic Branislav.

  model@summary.EM <- as.data.frame(do.call("rbind", summary.EM), stringsAsFactors = FALSE)
  
### End  

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
    
### Panic Branislav.

  colnames(model@summary.EM) <- c("Dataset",
    "strategy",
    "variant",
    "acceleration",
    "acceleration.multiplier",
    "tolerance",
    "maximum.iterations",
    "opt.iterations.nbr",
    "total.iterations.nbr")
    
### End   

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## REBMIX

setMethod("REBMIX",
          signature(model = "REBMVNORM"),
function(model, ...)
{
  summary <- list()
  
### Panic Branislav.

  summary.EM <- list()
  
### End  

  for (i in 1:length(model@Dataset)) {
    Dataset.name <- names(model@Dataset)[i]

    X <- as.matrix(model@Dataset[[i]])

    message("Dataset = ", Dataset.name)

    flush.console()

    n <- nrow(X)
    d <- ncol(X)
    
    length.pdf <- d

    if (is.character(model@K)) {
      if (model@Preprocessing == .rebmix$Preprocessing[3]) {
        K <- as.integer(sqrt(n) * 0.75):as.integer(sqrt(n) * 1.25)
      }
      else {
        Sturges <- ceiling(1.0 + log2(n)); RootN <- ceiling(2.0 * n^0.5)

        K <- kseq(Sturges, RootN, f = 0.25)
      }
      
      K <- unique(K) 
      
      dK <- max(K) - min(K) + 1; length.K <- length(K); K <- rep(K, d)
    }
    else {
      if (is.matrix(model@K)) {
        K <- model@K[i, ]
        
        dK <- 1; length.K <- 1
      }
      else {
        K <- sort(unique(model@K))
        
        dK <- max(K) - min(K) + 1; length.K <- length(K); K <- rep(K, d)
      }
    }

    if (length(model@theta1) > 0) {
      length.theta1 <- +d; theta1 <- model@theta1; theta1[is.na(theta1)] <- 0
    }
    else {
      length.theta1 <- -d; theta1 <- numeric()
    }

    if (length(model@theta2) > 0) {
      length.theta2 <- +d * d; length.theta3 <- +1; theta2 <- model@theta2; theta2[is.na(theta2)] <- 0
    }
    else {
      length.theta2 <- -d * d; length.theta3 <- -1; theta2 <- numeric()
    }
    
    output <- .C(C_RREBMVNORM,
      Preprocessing = as.character(model@Preprocessing),
      cmax = as.integer(model@cmax),
      cmin = as.integer(model@cmin),
      Criterion = as.character(model@Criterion),
      d = as.integer(d),
      Variables = as.character(model@Variables),
      length.pdf = as.integer(length.pdf),
      pdf = as.character(model@pdf),
      length.Theta = as.integer(4),
      length.theta = as.integer(c(length.theta1, length.theta2, length.theta2, length.theta3)),
      Theta = as.double(c(theta1, theta2)),
      length.K = as.integer(length.K),
      K = as.integer(K),
      length.y0 = as.integer(length(model@y0)),
      y0 = as.double(model@y0),
      length.ymin = as.integer(length(model@ymin)),
      ymin = as.double(model@ymin),
      length.ymax = as.integer(length(model@ymax)),
      ymax = as.double(model@ymax),
      ar = as.double(model@ar),
      Restraints = as.character(model@Restraints),
      n = as.integer(n),
      Y = as.double(X),
### Panic Branislav.      
      EMStrategy = as.character(model@EMcontrol@strategy),
      EMVariant = as.character(model@EMcontrol@variant),
      EMAcceleration = as.character(model@EMcontrol@acceleration),
      EMTolerance = as.double(model@EMcontrol@tolerance),
      EMAccelerationMul = as.double(model@EMcontrol@acceleration.multiplier),
      EMMaxIter = as.integer(model@EMcontrol@maximum.iterations),
      n_iter = integer(1),
      n_iter_all = integer(1),
### End         
      summary.k = integer(1),
      summary.h = double(d),
      summary.y0 = double(d),
      summary.ymin = double(d),
      summary.ymax = double(d),
      summary.IC = double(1),
      summary.logL = double(1),
      summary.M = integer(1),
      summary.c = integer(1),
      w = double(model@cmax),
      theta1 = double(model@cmax * d),
      theta2 = double(model@cmax * d * d),
      opt.length = integer(1),
      opt.c = integer(1000), # 1000 = ItMax see rebmixf.h
      opt.IC = double(1000),
      opt.logL = double(1000),
      opt.D = double(1000),
      all.length = integer(1),
      all.K = integer(dK),
      all.IC = double(dK),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RREBMVNORM!", call. = FALSE); return(NA)
    }

    c <- output$summary.c

    length(output$summary.h) <- d
    length(output$w) <- c
    length(output$theta1) <- c * d
    length(output$theta2) <- c * d * d

    length(output$opt.c) <- output$opt.length
    length(output$opt.IC) <- output$opt.length
    length(output$opt.logL) <- output$opt.length
    length(output$opt.D) <- output$opt.length

    j <- order(output$opt.c, output$opt.logL)

    output$opt.c <- output$opt.c[j]
    output$opt.IC <- output$opt.IC[j]
    output$opt.logL <- output$opt.logL[j]
    output$opt.D <- output$opt.D[j]

    j <- !duplicated(output$opt.c, fromLast = TRUE)

    output$opt.c <- output$opt.c[j]
    output$opt.IC <- output$opt.IC[j]
    output$opt.logL <- output$opt.logL[j]
    output$opt.D <- output$opt.D[j]

    length(output$all.K) <- output$all.length
    length(output$all.IC) <- output$all.length

    model@w[[i]] <- output$w

    model@Theta[[i]] <- list()

    length(model@Theta[[i]]) <- 3 * c

    names(model@Theta[[i]])[seq(1, 3 * c, 3)] <- paste("pdf", 1:c, sep = "")
    names(model@Theta[[i]])[seq(2, 3 * c, 3)] <- paste("theta1.", 1:c, sep = "")
    names(model@Theta[[i]])[seq(3, 3 * c, 3)] <- paste("theta2.", 1:c, sep = "")

    M <- which(pdf %in% .rebmix$pdf[.rebmix$pdf.nargs == 1])

    for (j in 1:c) {
	  model@Theta[[i]][[1 + (j - 1) * 3]] <- model@pdf
	  model@Theta[[i]][[2 + (j - 1) * 3]] <- output$theta1[seq((j - 1) * d + 1, j * d, 1)]
  	  model@Theta[[i]][[3 + (j - 1) * 3]] <- output$theta2[seq((j - 1) * d * d + 1, j * d * d, 1)]

      model@Theta[[i]][[3 + (j - 1) * 3]][M] <- NA
    }

    output$K <- paste("c(", paste(K, collapse = ","), ")", sep = "")

    if (model@Preprocessing == .rebmix$Preprocessing[1]) {
      length(output$summary.y0) <- d

      summary[[i]] <- c(Dataset.name,
        output$Preprocessing,
        output$cmax,
        output$cmin,
        output$Criterion,
        output$ar,
        output$Restraints,
        output$summary.c,
        output$summary.k,
        output$K,
        output$summary.y0,
        output$summary.ymin,
        output$summary.ymax,
        output$summary.h,
        output$summary.IC,
        output$summary.logL,
        output$summary.M)
    }
    else
    if (model@Preprocessing == .rebmix$Preprocessing[2]) {
      summary[[i]] <- c(Dataset.name,
        output$Preprocessing,
        output$cmax,
        output$cmin,
        output$Criterion,
        output$ar,
        output$Restraints,
        output$summary.c,
        output$summary.k,
        output$K,
        rep(NA, d),
        output$summary.ymin,
        output$summary.ymax,
        output$summary.h,
        output$summary.IC,
        output$summary.logL,
        output$summary.M)
    }
    if (model@Preprocessing == .rebmix$Preprocessing[3]) {
      summary[[i]] <- c(Dataset.name,
        output$Preprocessing,
        output$cmax,
        output$cmin,
        output$Criterion,
        output$ar,
        output$Restraints,
        output$summary.c,
        output$summary.k,
        output$K,
        rep(NA, d),
        output$summary.ymin,
        output$summary.ymax,
        output$summary.h,
        output$summary.IC,
        output$summary.logL,
        output$summary.M)
    }
    
### Panic Branislav.

    summary.EM[[i]] <- c(Dataset.name,
        output$EMStrategy,
        output$EMVariant,
        output$EMAcceleration,
        output$EMAccelerationMul,
        output$EMTolerance,
        output$EMMaxIter,
        output$n_iter,
        output$n_iter_all)
        
### End    

    model@opt.c[[i]] <- output$opt.c
    model@opt.IC[[i]] <- output$opt.IC
    model@opt.logL[[i]] <- output$opt.logL
    model@opt.D[[i]] <- output$opt.D
    model@all.K[[i]] <- output$all.K
    model@all.IC[[i]] <- output$all.IC
  }

  model@summary <- as.data.frame(do.call("rbind", summary), stringsAsFactors = FALSE)
  
### Panic Branislav.

  model@summary.EM <- as.data.frame(do.call("rbind", summary.EM), stringsAsFactors = FALSE)
  
### End  

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
    
### Panic Branislav.

  colnames(model@summary.EM) <- c("Dataset",
    "strategy",
    "variant",
    "acceleration",
    "acceleration.multiplier",
    "tolerance",
    "maximum.iterations",
    "opt.iterations.nbr",
    "total.iterations.nbr")
    
### End       

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## REBMIX

setMethod("REBMIX",
          signature(model = "ANY"),
function(model,
  Dataset,
  Preprocessing,
  cmax,
  cmin,
  Criterion,
  pdf,
  theta1,
  theta2,
  K,
  y0,
  ymin,
  ymax,
  ar,
  Restraints,
### Panic Branislav.
  EMcontrol, ...)
### End    
{
  digits <- getOption("digits"); options(digits = 15)

  message("REBMIX Version 2.12.1")

  flush.console()

  model <- new(model,
     Dataset = Dataset,
     Preprocessing = Preprocessing,
     cmax = cmax,
     cmin = cmin,
     Criterion = Criterion,
     pdf = pdf,
     theta1 = theta1,
     theta2 = theta2,
     K = K,
     y0 = y0,
     ymin = ymin,
     ymax = ymax,
     ar = ar,
     Restraints = Restraints,
### Panic Branislav.
     EMcontrol = EMcontrol)
### End     

  output <- REBMIX(model = model, ...)

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
  
### Panic Branislav.  
  
  if (is.null(model@summary.EM)) {
    model@summary.EM <- output@summary.EM 
  }
  else {
    model@summary.EM <- merge(model@summary.EM, output@summary.EM, all = TRUE, sort = FALSE)
  }
  
### End    

  for (k in (1:length(Dataset))) {
    model@opt.c[[length(model@opt.c) + 1]] <- output@opt.c[[k]]
    model@opt.IC[[length(model@opt.IC) + 1]] <- output@opt.IC[[k]]
    model@opt.logL[[length(model@opt.logL) + 1]] <- output@opt.logL[[k]]
    model@opt.D[[length(model@opt.D) + 1]] <- output@opt.D[[k]]
    model@all.K[[length(model@all.K) + 1]] <- output@all.K[[k]]
    model@all.IC[[length(model@all.IC) + 1]] <- output@all.IC[[k]]
  }

  model@pos <- which(as.numeric(model@summary[, "logL"]) == max(as.numeric(model@summary[, "logL"])))

  if (model@cmax <= max(as.numeric(model@summary[, "c"]))) {
    message("Note: c = ", cmax, "! Consider increasing ", sQuote("cmax"), "!")
  }

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## REBMIX
