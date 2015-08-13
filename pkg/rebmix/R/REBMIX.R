.REBMIX <- function(Dataset = NULL, 
  Preprocessing = NULL, 
  cmax = 15,
  Criterion = "AIC",
  Variables = NULL,
  pdf = NULL,
  Theta1 = NULL,
  Theta2 = NULL,
  K = NULL,
  y0 = NULL,
  ymin = NULL,
  ymax = NULL,
  ar = 0.1,
  Restraints = "loose", ...)
{
  REBMIX <- NULL
  REBMIX$Dataset <- Dataset
  REBMIX$w <- list()
  REBMIX$Theta <- list()
  REBMIX$summary <- list()
  REBMIX$pos <- 1
  REBMIX$opt.c <- list()
  REBMIX$opt.IC <- list()
  REBMIX$opt.logL <- list()
  REBMIX$opt.D <- list()
  REBMIX$all.length <- list()  
  REBMIX$all.K <- list()
  REBMIX$all.IC <- list()  

  for (i in 1:length(Dataset)) {
    DatasetName <- names(Dataset)[i]

    X <- as.matrix(Dataset[[i]])

    message("Dataset = ", DatasetName)
    flush.console()

    n <- nrow(X)
    d <- ncol(X)

    if (d < 1) {
      stop(sQuote("d"), " must be greater than 0!", call. = FALSE)
    }

    if (n < 1) {
      stop(sQuote("n"), " must be greater than 0!", call. = FALSE)
    }
    
    if (length(pdf) != d) {
      stop("lengths of ", sQuote("pdf"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }    

    if (!is.null(Theta1)) {
      Theta1[is.na(Theta1)] <- 0
      
      if (length(Theta1) != d) {
        stop("lengths of ", sQuote("Theta1"), " and ", sQuote("d"), " must match!", call. = FALSE)
      }
    }

    if (!is.null(Theta2)) {
      Theta2[is.na(Theta2)] <- 0
      
      if (length(Theta2) != d) {
        stop("lengths of ", sQuote("Theta2"), " and ", sQuote("d"), " must match!", call. = FALSE)
      }
    }
    
    if (!is.null(y0)) {
      y0[is.na(y0)] <- 0
      
      if (length(y0) != d) {
        stop("lengths of ", sQuote("y0"), " and ", sQuote("d"), " must match!", call. = FALSE)
      }      
    }    

    if (!is.null(ymin)) {
      ymin[is.na(ymin)] <- 0
      
      if (length(ymin) != d) {
        stop("lengths of ", sQuote("ymin"), " and ", sQuote("d"), " must match!", call. = FALSE)
      }      
    }

    if (!is.null(ymax)) {
      ymax[is.na(ymax)] <- 0
      
      if (length(ymax) != d) {
        stop("lengths of ", sQuote("ymax"), " and ", sQuote("d"), " must match!", call. = FALSE)
      }      
    }
    
    if (as.integer(length(pdf)) > 0) {
        length.pdf = +d;
    }
    else {
        length.pdf = -d;
    }
    
    if (as.integer(length(Theta1)) > 0) {
        length.Theta1 = +d;
    }
    else {
        length.Theta1 = -d;
    }    
    
    if (as.integer(length(Theta2)) > 0) {
        length.Theta2 = +d;
    }
    else {
        length.Theta2 = -d;
    }      
    
    output <- .C("RREBMIX",
      Preprocessing = as.character(Preprocessing), 
      cmax = as.integer(cmax),
      Criterion = as.character(Criterion),
      d = as.integer(d),
      Variables = as.character(Variables),
      length.pdf = length.pdf,
      pdf = as.character(pdf),
      length.Theta1 = length.Theta1,
      Theta1 = as.double(Theta1),
      length.Theta2 = length.Theta2,
      Theta2 = as.double(Theta2),
      length.K = as.integer(length(K)),
      K = as.integer(K),
      length.y0 = as.integer(length(y0)),
      y0 = as.double(y0),      
      length.ymin = as.integer(length(ymin)),
      ymin = as.double(ymin),
      length.ymax = as.integer(length(ymax)),
      ymax = as.double(ymax),
      ar = as.double(ar),
      Restraints = as.character(Restraints),
      n = as.integer(n),
      Y = as.double(X),
      summary.k = integer(1),
      summary.h = double(d),
      summary.y0 = double(d),      
      summary.IC = double(1),
      summary.logL = double(1),
      summary.M = integer(1), 
      summary.c = integer(1),
      W = double(cmax),        
      Theta.pdf = as.character(rep("THE_LONGEST_PARAMETRIC_FAMILY_TYPE", cmax * d)),
      Theta.Theta1 = double(cmax * d),        
      Theta.Theta2 = double(cmax * d),  
      opt.length = integer(1),
      opt.c = integer(1000), ## 1000 = ItMax see rebmixf.h
      opt.IC = double(1000),
      opt.logL = double(1000),
      opt.D = double(1000),
      all.length = integer(1),
      all.K = integer(max(K) - min(K) + 1),
      all.IC = double(max(K) - min(K) + 1),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in REBMIX!", call. = FALSE); return(NA)
    }

    c <- output$summary.c

    length(output$summary.h) <- d
    length(output$W) <- c
    length(output$Theta.pdf) <- c * d
    length(output$Theta.Theta1) <- c * d
    length(output$Theta.Theta2) <- c * d

    dim(output$Theta.pdf) <- c(d, c)
    dim(output$Theta.Theta1) <- c(d, c)
    dim(output$Theta.Theta2) <- c(d, c)
    
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

    REBMIX$w[[i]] <- as.data.frame(rbind(output$W), stringsAsFactors = FALSE)

    rownames(REBMIX$w[[i]]) <- "w"
    colnames(REBMIX$w[[i]]) <- paste("comp", if (c > 1) 1:c else "", sep = "")

    REBMIX$Theta[[i]] <- rbind(output$Theta.pdf, output$Theta.Theta1, output$Theta.Theta2)

    dim(REBMIX$Theta[[i]]) <- c(3 * d, c)

    if (d > 1) {
      rownames(REBMIX$Theta[[i]]) <- c(paste("pdf", 1:d, sep = ""),
        paste("theta1.", 1:d, sep = ""),
        paste("theta2.", 1:d, sep = ""))
    }
    else {
      rownames(REBMIX$Theta[[i]]) <- c("pdf", "theta1", "theta2")
    }

    Index <- NULL

    for (j in 1:d){
      Index <- c(Index, seq(from = j, to = j + 2 * d, by = d))
    }

    REBMIX$Theta[[i]] <- cbind(REBMIX$Theta[[i]][Index, ])

    M <- match(REBMIX$Theta[[i]][, 1], .rebmix$pdf)

    Index <- NULL

    for (j in 1:length(M)) {
      if (M[j] %in% which(.rebmix$pdf.nargs == 1)) {
        Index <- c(Index, j + 2)
      }
    }

    if (is.null(Index)) {
      REBMIX$Theta[[i]] <- as.data.frame(REBMIX$Theta[[i]], stringsAsFactors = FALSE)
    }
    else {
      REBMIX$Theta[[i]] <- as.data.frame(REBMIX$Theta[[i]][-Index, ], stringsAsFactors = FALSE)
    }

    colnames(REBMIX$Theta[[i]]) <- paste("comp", if (c > 1) 1:c else "", sep = "")
    
    output$K <- paste("c(", paste(K, collapse = ","), ")", sep = "")
    
    if (Preprocessing == .rebmix$Preprocessing[1]) {
      length(output$summary.y0) <- d    
    
      REBMIX$summary[[i]] <- c(DatasetName, 
        output$Preprocessing, 
        output$cmax, 
        output$Criterion, 
        output$ar,
        output$Restraints,
        output$summary.c,
        output$summary.k,
        output$K,
        output$summary.y0,
        output$summary.h,
        output$summary.IC,
        output$summary.logL,
        output$summary.M)
    }
    else
    if (Preprocessing == .rebmix$Preprocessing[2]) {
      REBMIX$summary[[i]] <- c(DatasetName, 
        output$Preprocessing, 
        output$cmax, 
        output$Criterion, 
        output$ar,
        output$Restraints,
        output$summary.c,
        output$summary.k,
        output$K,
        output$summary.h,
        output$summary.IC,
        output$summary.logL,
        output$summary.M)        
    }
    if (Preprocessing == .rebmix$Preprocessing[3]) {
      REBMIX$summary[[i]] <- c(DatasetName, 
        output$Preprocessing, 
        output$cmax, 
        output$Criterion, 
        output$ar,
        output$Restraints,
        output$summary.c,
        output$summary.k,
        output$K,
        output$summary.h,
        output$summary.IC,
        output$summary.logL,
        output$summary.M)        
    }
    
    REBMIX$opt.c[[i]] <- output$opt.c;
    REBMIX$opt.IC[[i]] <- output$opt.IC;
    REBMIX$opt.logL[[i]] <- output$opt.logL;
    REBMIX$opt.D[[i]] <- output$opt.D;    
    REBMIX$all.length[[i]] <- output$all.length;
    REBMIX$all.K[[i]] <- output$all.K;
    REBMIX$all.IC[[i]] <- output$all.IC;     
  }

  REBMIX$summary <- as.data.frame(do.call("rbind", REBMIX$summary), stringsAsFactors = FALSE)

  if (Preprocessing == .rebmix$Preprocessing[1]) {
    colnames(REBMIX$summary) <- c("Dataset", 
      "Preprocessing", 
      "cmax", 
      "Criterion", 
      "ar", 
      "Restraints", 
      "c", 
      "v/k", 
      "K",       
      paste("y0", if (d > 1) 1:d else "", sep = ""), 
      paste("h", if (d > 1) 1:d else "", sep = ""), 
      "IC", 
      "logL",
      "M")
  }
  else
  if (Preprocessing == .rebmix$Preprocessing[2]) {
    colnames(REBMIX$summary) <- c("Dataset", 
      "Preprocessing", 
      "cmax", 
      "Criterion", 
      "ar", 
      "Restraints", 
      "c", 
      "v/k",
      "K", 
      paste("h", if (d > 1) 1:d else "", sep = ""), 
      "IC", 
      "logL",
      "M")
  }
  if (Preprocessing == .rebmix$Preprocessing[3]) {
    colnames(REBMIX$summary) <- c("Dataset", 
      "Preprocessing", 
      "cmax", 
      "Criterion", 
      "ar", 
      "Restraints", 
      "c", 
      "v/k",
      "K",        
      paste("h", if (d > 1) 1:d else "", sep = ""),
      "IC", 
      "logL",
      "M")
  }

  rm(list = ls()[!(ls() %in% c("REBMIX"))])

  class(REBMIX) <- "REBMIX"
 
  return(REBMIX)
} ## .REBMIX 

REBMIX <- function(Dataset = NULL, 
  Preprocessing = NULL, 
  cmax = 15,
  Criterion = "AIC",
  Variables = NULL,
  pdf = NULL,
  Theta1 = NULL,
  Theta2 = NULL,
  K = NULL,
  y0 = NULL,
  ymin = NULL,
  ymax = NULL,
  ar = 0.1,
  Restraints = "loose", ...)
{
  digits <- getOption("digits"); options(digits = 15)

  message("REBMIX Version 2.7.2");
  flush.console()

  if (is.null(Dataset)) {
    stop(sQuote("Dataset"), " must not be NULL!", call. = FALSE)
  }

  if (!is.list(Dataset)) {
    stop(sQuote("Dataset"), " list of data frames is requested!", call. = FALSE)
  }
  
  if (is.null(names(Dataset))) {
    names(Dataset) <- paste("dataset", 1:length(Dataset), sep = "")
  }
    
  for (i in 1:length(Dataset)) {
    if (!is.data.frame(Dataset[[i]])) {
      stop(sQuote("Dataset"), " list of data frames or character vector is requested!", call. = FALSE)
    }
      
    if ((is.na(names(Dataset)[i])) || (names(Dataset)[i] == "")) {
      names(Dataset)[i] <- paste("dataset", i, sep = "")  
    }
  }
  
  if (is.null(Preprocessing)) {
    stop(sQuote("Preprocessing"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(Preprocessing)) {
    stop(sQuote("Preprocessing"), " character vector is requested!", call. = FALSE)
  }

  Preprocessing <- match.arg(Preprocessing, .rebmix$Preprocessing, several.ok = TRUE)

  if (!is.wholenumber(cmax)) {
    stop(sQuote("cmax"), " integer is requested!", call. = FALSE)
  }

  if (cmax < 1) {
    stop(sQuote("cmax"), " must be greater than 0!", call. = FALSE)
  }

  if (!is.character(Criterion)) {
    stop(sQuote("Criterion"), " character vector is requested!", call. = FALSE)
  }

  Criterion <- match.arg(Criterion, .rebmix$Criterion, several.ok = TRUE)

  if (is.null(Variables)) {
    stop(sQuote("Variables"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(Variables)) {
    stop(sQuote("Variables"), " character vector is requested!", call. = FALSE)
  }

  Variables <- match.arg(Variables, .rebmix$Variables, several.ok = TRUE)

  if (is.null(pdf)) {
    stop(sQuote("pdf"), " must not be NULL!", call. = FALSE)
  }

  if (!is.character(pdf)) {
    stop(sQuote("pdf"), " character vector is requested!", call. = FALSE)
  }

  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)

  if (is.null(K)) {
    stop(sQuote("K"), " must not be NULL!", call. = FALSE)
  }

  if (is.list(K)) {
    for (i in 1:length(K)) {
      if (!is.wholenumber(K[[i]])) {
        stop(sQuote("K"), " integer vector is requested!", call. = FALSE)
      }

      if (!all(K[[i]] > 0)) {
        stop("all ", sQuote("K"), " must be greater than 0!", call. = FALSE)
      }
    }

    if (length(K) != length(Preprocessing)) {
      stop("lengths of ", sQuote("Preprocessing"), " and ", sQuote("K"), " must match!", call. = FALSE)
    }
  }
  else {
    if (!is.wholenumber(K)) {
      stop(sQuote("K"), " integer vector is requested!", call. = FALSE)
    }

    if (!all(K > 0)) {
      stop("all ", sQuote("K"), " must be greater than 0!", call. = FALSE)
    }
  }

  if (!is.numeric(ar)) {
    stop(sQuote("ar"), " numeric is requested!", call. = FALSE)
  }

  if ((ar <= 0.0) || (ar > 1.0)) {
    stop(sQuote("ar"), " must be greater than 0.0 and less or equal than 1.0!", call. = FALSE)
  }

  if (!is.character(Restraints)) {
    stop(sQuote("Restraints"), " character is requested!", call. = FALSE)
  }

  Restraints <- match.arg(Restraints, .rebmix$Restraints, several.ok = FALSE)

  REBMIX <- NULL
  REBMIX$Dataset <- Dataset
  REBMIX$w <- list()
  REBMIX$Theta <- list()
  REBMIX$summary <- NULL
  REBMIX$pos <- 1
  REBMIX$opt.c <- list()
  REBMIX$opt.IC <- list()
  REBMIX$opt.logL <- list()
  REBMIX$opt.D <- list()
  REBMIX$all.length <- list()
  REBMIX$all.K <- list()
  REBMIX$all.IC <- list()  
  
  REBMIX$call <- list( 
    Preprocessing = Preprocessing, 
    cmax = cmax,
    Criterion = Criterion,
    Variables = Variables,
    pdf = pdf,
    Theta1 = Theta1,
    Theta2 = Theta2,
    K = K,
    y0 = y0,
    ymin = ymin,
    ymax = ymax,
    ar = ar,
    Restraints = Restraints)

  for (i in 1:length(Preprocessing)) {
    for (j in 1:length(Criterion)) {
      output <- .REBMIX(Dataset = Dataset, 
        Preprocessing = Preprocessing[i], 
        cmax = cmax,
        Criterion = Criterion[j],
        Variables = Variables,
        pdf = pdf,
        Theta1 = Theta1,
        Theta2 = Theta2,
        K = if (is.list(K)) K[[i]] else K,
        y0 = y0,        
        ymin = ymin,
        ymax = ymax,
        ar = ar,
        Restraints = Restraints, ...)

      for (k in (1:length(Dataset))) {
        REBMIX$w[[length(REBMIX$w) + 1]] <- output$w[[k]] 
        REBMIX$Theta[[length(REBMIX$Theta) + 1]] <- output$Theta[[k]]
      }

      if (is.null(REBMIX$summary)) {
        REBMIX$summary <- output$summary
      }
      else {
        REBMIX$summary <- merge(REBMIX$summary, output$summary, all = TRUE, sort = FALSE)
      }
      
      for (k in (1:length(Dataset))) {
        REBMIX$opt.c[[length(REBMIX$opt.c) + 1]] <- output$opt.c[[k]] 
        REBMIX$opt.IC[[length(REBMIX$opt.IC) + 1]] <- output$opt.IC[[k]] 
        REBMIX$opt.logL[[length(REBMIX$opt.logL) + 1]] <- output$opt.logL[[k]] 
        REBMIX$opt.D[[length(REBMIX$opt.D) + 1]] <- output$opt.D[[k]] 
        REBMIX$all.length[[length(REBMIX$all.length) + 1]] <- output$all.length[[k]] 
        REBMIX$all.K[[length(REBMIX$all.K) + 1]] <- output$all.K[[k]] 
        REBMIX$all.IC[[length(REBMIX$all.IC) + 1]] <- output$all.IC[[k]] 
      }      
    }
  }
  
  REBMIX$pos <- which(as.numeric(REBMIX$summary[, "logL"]) == max(as.numeric(REBMIX$summary[, "logL"])))  
  
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("REBMIX"))])

  class(REBMIX) <- "REBMIX"
  
  return(REBMIX)
} ## REBMIX
