setMethod("REBMIX",
          signature(model = "REBMIX"),
function(model, ...)
{
  summary <- list()

### Panic Branislav.
  summary.EM <- list()
### End  

  for (i in 1:length(model@Dataset)) {
    if (as.character(class(model@Dataset[[i]])) == "data.frame") {
      Y.type <- 0
      
      Dataset.name <- names(model@Dataset)[i]
      
      X <- as.matrix(model@Dataset[[i]])
      
      d <- ncol(X)
      
      h <- NULL
    }
    else
    if (as.character(class(model@Dataset[[i]])) == "Histogram") {
      Y.type <- 1
      
      Dataset.name <- names(model@Dataset)[i]
      
      X <- as.matrix(model@Dataset[[i]]@Y)
      
      d <- ncol(X) - 1
      
      h <- model@Dataset[[i]]@h
    }

    message("Dataset = ", Dataset.name)

    flush.console()

    n <- nrow(X)
    
    length.pdf <- d
    
    if (Y.type == 0) {
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
    }
    else
    if (Y.type == 1) {
      dK <- 1; length.K <- 0; K <- NULL 
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
    
    if (length(model@theta3) > 0) {
      length.theta3 <- +d; theta3 <- model@theta3; theta3[is.na(theta3)] <- 0
    }
    else {
      length.theta3 <- -d; theta3 <- numeric()
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
      length.Theta = as.integer(3),
      length.theta = as.integer(c(length.theta1, length.theta2, length.theta3)),
      Theta = as.double(c(theta1, theta2, theta3)),
      length.K = as.integer(length.K),
      K = as.integer(K),
      length.ymin = as.integer(length(model@ymin)),
      ymin = as.double(model@ymin),
      length.ymax = as.integer(length(model@ymax)),
      ymax = as.double(model@ymax),
      length.h = as.integer(length(h)),
      h = as.double(h),      
      ar = as.double(model@ar),
      Restraints = as.character(model@Restraints),
      Mode = as.character(model@Mode),
      n = as.integer(n),
      Y = as.double(X),
      Y.type = as.integer(Y.type),
### Panic Branislav.      
      EMStrategy = as.character(model@EMcontrol@strategy),
      EMVariant = as.character(model@EMcontrol@variant),
      EMAcceleration = as.character(model@EMcontrol@acceleration),
      EMTolerance = as.double(model@EMcontrol@tolerance),
      EMAccelerationMul = as.double(model@EMcontrol@acceleration.multiplier),
      EMMaxIter = as.integer(model@EMcontrol@maximum.iterations),
      EMK = as.integer(model@EMcontrol@K),
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
      theta3 = double(model@cmax * d),
      opt.length = integer(1),
      opt.c = integer(1000), # 1000 = ItMax see rebmixf.h
      opt.IC = double(1000),
      opt.logL = double(1000),
      opt.Dmin = double(1000),
      opt.D = double(1000),
      all.length = integer(1),
      all.K = integer(dK),
      all.IC = double(dK),
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

    c <- output$summary.c

    length(output$summary.h) <- d
    length(output$w) <- c
    length(output$theta1) <- c * d
    length(output$theta2) <- c * d
    length(output$theta3) <- c * d

    length(output$opt.c) <- output$opt.length
    length(output$opt.IC) <- output$opt.length
    length(output$opt.logL) <- output$opt.length
    length(output$opt.Dmin) <- output$opt.length
    length(output$opt.D) <- output$opt.length

    j <- order(output$opt.c, output$opt.logL)

    output$opt.c <- output$opt.c[j]
    output$opt.IC <- output$opt.IC[j]
    output$opt.logL <- output$opt.logL[j]
    output$opt.Dmin <- output$opt.Dmin[j]
    output$opt.D <- output$opt.D[j]

    j <- !duplicated(output$opt.c, fromLast = TRUE)

    output$opt.c <- output$opt.c[j]
    output$opt.IC <- output$opt.IC[j]
    output$opt.logL <- output$opt.logL[j]
    output$opt.Dmin <- output$opt.Dmin[j]
    output$opt.D <- output$opt.D[j]

    length(output$all.K) <- output$all.length
    length(output$all.IC) <- output$all.length

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
	  model@Theta[[i]][[2 + (j - 1) * 4]] <- output$theta1[seq((j - 1) * d + 1, j * d, 1)]
  	  model@Theta[[i]][[3 + (j - 1) * 4]] <- output$theta2[seq((j - 1) * d + 1, j * d, 1)]
  	  model@Theta[[i]][[4 + (j - 1) * 4]] <- output$theta3[seq((j - 1) * d + 1, j * d, 1)]

      model@Theta[[i]][[3 + (j - 1) * 4]][M1] <- NA
      model@Theta[[i]][[4 + (j - 1) * 4]][M2] <- NA
    }

    output$K <- paste("c(", paste(K, collapse = ","), ")", sep = "")

    if (Y.type == 0) {
      if (model@Preprocessing == .rebmix$Preprocessing[1]) {
        length(output$summary.y0) <- d

        summary[[i]] <- c(Dataset.name,
          output$Preprocessing,
          output$cmax,
          output$cmin,
          output$Criterion,
          output$ar,
          output$Restraints,
          output$Mode,
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
          output$Mode,
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
          output$Mode,
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
    }
    else
    if (Y.type == 1) {
      summary[[i]] <- c(Dataset.name,
        NA,
        output$cmax,
        output$cmin,
        output$Criterion,
        output$ar,
        output$Restraints,
        output$Mode,
        output$summary.c,
        NA,
        NA,
        rep(NA, d),
        rep(NA, d),
        rep(NA, d),
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
    model@opt.Dmin[[i]] <- output$opt.Dmin
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
    "Mode",
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
    if (as.character(class(model@Dataset[[i]])) == "data.frame") {
      Y.type <- 0
      
      Dataset.name <- names(model@Dataset)[i]
      
      X <- as.matrix(model@Dataset[[i]])
      
      d <- ncol(X)
      
      h <- NULL
    }
    else
    if (as.character(class(model@Dataset[[i]])) == "Histogram") {
      Y.type <- 1
      
      Dataset.name <- names(model@Dataset)[i]
      
      X <- as.matrix(model@Dataset[[i]]@Y)
      
      d <- ncol(X) - 1
      
      h <- model@Dataset[[i]]@h
    }

    message("Dataset = ", Dataset.name)

    flush.console()

    n <- nrow(X)
    
    length.pdf <- d

    if (Y.type == 0) {
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
    }
    else
    if (Y.type == 1) {
      dK <- 1; length.K <- 0; K <- NULL
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
      length.ymin = as.integer(length(model@ymin)),
      ymin = as.double(model@ymin),
      length.ymax = as.integer(length(model@ymax)),
      ymax = as.double(model@ymax),
      length.h = as.integer(length(h)),
      h = as.double(h),       
      ar = as.double(model@ar),
      Restraints = as.character(model@Restraints),
      Mode = as.character(model@Mode),
      n = as.integer(n),
      Y = as.double(X),
      Y.type = as.integer(Y.type),
### Panic Branislav.      
      EMStrategy = as.character(model@EMcontrol@strategy),
      EMVariant = as.character(model@EMcontrol@variant),
      EMAcceleration = as.character(model@EMcontrol@acceleration),
      EMTolerance = as.double(model@EMcontrol@tolerance),
      EMAccelerationMul = as.double(model@EMcontrol@acceleration.multiplier),
      EMMaxIter = as.integer(model@EMcontrol@maximum.iterations),
      EMK = as.integer(model@EMcontrol@K),      
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
      opt.Dmin = double(1000),
      opt.D = double(1000),
      all.length = integer(1),
      all.K = integer(dK),
      all.IC = double(dK),
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

    c <- output$summary.c

    length(output$summary.h) <- d
    length(output$w) <- c
    length(output$theta1) <- c * d
    length(output$theta2) <- c * d * d

    length(output$opt.c) <- output$opt.length
    length(output$opt.IC) <- output$opt.length
    length(output$opt.logL) <- output$opt.length
    length(output$opt.Dmin) <- output$opt.length
    length(output$opt.D) <- output$opt.length

    j <- order(output$opt.c, output$opt.logL)

    output$opt.c <- output$opt.c[j]
    output$opt.IC <- output$opt.IC[j]
    output$opt.logL <- output$opt.logL[j]
    output$opt.Dmin <- output$opt.Dmin[j]
    output$opt.D <- output$opt.D[j]

    j <- !duplicated(output$opt.c, fromLast = TRUE)

    output$opt.c <- output$opt.c[j]
    output$opt.IC <- output$opt.IC[j]
    output$opt.logL <- output$opt.logL[j]
    output$opt.Dmin <- output$opt.Dmin[j]
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

    if (Y.type == 0) {
      if (model@Preprocessing == .rebmix$Preprocessing[1]) {
        length(output$summary.y0) <- d

        summary[[i]] <- c(Dataset.name,
          output$Preprocessing,
          output$cmax,
          output$cmin,
          output$Criterion,
          output$ar,
          output$Restraints,
          output$Mode,
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
          output$Mode,
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
          output$Mode,
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
    }
    else
    if (Y.type == 1) {
      summary[[i]] <- c(Dataset.name,
        NA,
        output$cmax,
        output$cmin,
        output$Criterion,
        output$ar,
        output$Restraints,
        output$Mode,
        output$summary.c,
        NA,
        NA,
        rep(NA, d),
        rep(NA, d),
        rep(NA, d),
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
    model@opt.Dmin[[i]] <- output$opt.Dmin
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
    "Mode",
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
  theta3,  
  K,
  ymin,
  ymax,
  ar,
  Restraints,
  Mode,
### Panic Branislav.
  EMcontrol, ...)
### End    
{
  digits <- getOption("digits"); options(digits = 15)

  message("REBMIX Version 2.16.1")

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
     theta3 = theta3,
     K = K,
     ymin = ymin,
     ymax = ymax,
     ar = ar,
     Restraints = Restraints,
     Mode = Mode,
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
    model@opt.Dmin[[length(model@opt.Dmin) + 1]] <- output@opt.Dmin[[k]]
    model@opt.D[[length(model@opt.D) + 1]] <- output@opt.D[[k]]
    model@all.K[[length(model@all.K) + 1]] <- output@all.K[[k]]
    model@all.IC[[length(model@all.IC) + 1]] <- output@all.IC[[k]]
  }

  model@pos <- which(as.numeric(model@summary[, "logL"]) == max(as.numeric(model@summary[, "logL"])))

  if ((model@cmax > model@cmin) && (model@cmax <= max(as.numeric(model@summary[, "c"])))) {
    message("Note: c = ", cmax, "! Consider increasing ", sQuote("cmax"), "!")
  }

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## REBMIX
