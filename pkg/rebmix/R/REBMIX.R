.REBMIX <- function(Dataset = NULL, 
  Preprocessing = NULL, 
  D = 0.025, 
  cmax = 15,
  Criterion = "AIC",
  Variables = NULL,
  pdf = NULL,
  Theta1 = NULL,
  Theta2 = NULL,
  K = NULL,
  ymin = NULL,
  ymax = NULL,
  b = 1.0,
  ar = 0.1,
  Restraints = "loose")
{
  REBMIX <- NULL
  REBMIX$Dataset <- Dataset
  REBMIX$w <- list()
  REBMIX$Theta <- list()
  REBMIX$Variables <- Variables
  REBMIX$pdf <- pdf
  REBMIX$Theta1 <- Theta1
  REBMIX$Theta2 <- Theta2
  REBMIX$summary <- list()
  REBMIX$pos <- 1
  REBMIX$all.Imax <- list()
  REBMIX$all.c <- list()
  REBMIX$all.IC <- list()
  REBMIX$all.logL <- list()
  REBMIX$all.D <- list()

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

    output <- .C("RREBMIX",
      PreType = as.character(Preprocessing), 
      D = as.double(D),
      cmax = as.integer(cmax),
      ICType = as.character(Criterion),
      d = as.integer(d),
      VarType = as.character(Variables),
      IniFamType = as.character(pdf),
      length.Ini0 = as.integer(length(Theta1)),
      Ini0 = as.double(Theta1),
      length.Ini1 = as.integer(length(Theta2)),
      Ini1 = as.double(Theta2),
      kmax = as.integer(length(K)),
      K = as.integer(K),
      length.ymin = as.integer(length(ymin)),
      ymin = as.double(ymin),
      length.ymax = as.integer(length(ymax)),
      ymax = as.double(ymax),
      b = as.double(b),
      ar = as.double(ar),
      ResType = as.character(Restraints),
      n = as.integer(n),
      X = as.double(X),
      k = integer(1),
      h = double(d),
      y0 = double(d),
      IC = double(1),
      logL = double(1),
      df = integer(1), 
      c = integer(1),
      W = double(cmax),        
      ParFamType = as.character(rep("THE_LONGEST_PARAMETRIC_FAMILY_TYPE", cmax * d)),
      Par0 = double(cmax * d),        
      Par1 = double(cmax * d),  
      all.Imax = integer(1),
      all.c = integer(1000), ## 1000 = ItMax see rebmixf.h
      all.IC = double(1000),
      all.logL = double(1000),
      all.D = double(1000),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in REBMIX!", call. = FALSE); return(NA)
    }

    c <- output$c

    length(output$h) <- d
    length(output$y0) <- d
    length(output$W) <- c
    length(output$ParFamType) <- c * d
    length(output$Par0) <- c * d
    length(output$Par1) <- c * d

    dim(output$ParFamType) <- c(d, c)
    dim(output$Par0) <- c(d, c)
    dim(output$Par1) <- c(d, c)
    
    length(output$all.c) <- output$all.Imax 
    length(output$all.IC) <- output$all.Imax   
    length(output$all.logL) <- output$all.Imax 
    length(output$all.D) <- output$all.Imax 

    REBMIX$w[[i]] <- as.data.frame(rbind(output$W), stringsAsFactors = FALSE)

    rownames(REBMIX$w[[i]]) <- "w"
    colnames(REBMIX$w[[i]]) <- paste("comp", if (c > 1) 1:c else "", sep = "")

    REBMIX$Theta[[i]] <- rbind(output$ParFamType, output$Par0, output$Par1)

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

    if (Preprocessing == .rebmix$Preprocessing[1]) {
      REBMIX$summary[[i]] <- c(DatasetName, 
        output$PreType, 
        output$D,
        output$cmax, 
        output$ICType, 
        output$ar,
        output$ResType,
        output$c,
        output$b,
        output$k,
        output$y0,
        output$h,
        output$IC,
        output$logL,
        output$df)
    }
    else
    if (Preprocessing == .rebmix$Preprocessing[2]) {
      REBMIX$summary[[i]] <- c(DatasetName, 
        output$PreType, 
        output$D,
        output$cmax, 
        output$ICType, 
        output$ar,
        output$ResType,
        output$c,
        output$b,
        output$k,
        output$h,
        output$IC,
        output$logL,
        output$df)        
    }
    if (Preprocessing == .rebmix$Preprocessing[3]) {
      REBMIX$summary[[i]] <- c(DatasetName, 
        output$PreType, 
        output$D,
        output$cmax, 
        output$ICType, 
        output$ar,
        output$ResType,
        output$c,
        output$b,
        output$k,
        output$h,
        output$IC,
        output$logL,
        output$df)        
    }
    
    REBMIX$all.Imax[[i]] <- output$all.Imax;
    REBMIX$all.c[[i]] <- output$all.c;
    REBMIX$all.IC[[i]] <- output$all.IC;
    REBMIX$all.logL[[i]] <- output$all.logL;
    REBMIX$all.D[[i]] <- output$all.D;     
  }

  REBMIX$summary <- as.data.frame(do.call("rbind", REBMIX$summary), stringsAsFactors = FALSE)

  if (Preprocessing == .rebmix$Preprocessing[1]) {
    colnames(REBMIX$summary) <- c("Dataset", 
      "Preprocessing", 
      "D", 
      "cmax", 
      "Criterion", 
      "ar", 
      "Restraints", 
      "c", 
      "b",
      "v/k", 
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
      "D", 
      "cmax", 
      "Criterion", 
      "ar", 
      "Restraints", 
      "c", 
      "b", 
      "v/k", 
      paste("h", if (d > 1) 1:d else "", sep = ""), 
      "IC", 
      "logL",
      "M")
  }
  if (Preprocessing == .rebmix$Preprocessing[3]) {
    colnames(REBMIX$summary) <- c("Dataset", 
      "Preprocessing", 
      "D", 
      "cmax", 
      "Criterion", 
      "ar", 
      "Restraints", 
      "c", 
      "b",
      "v/k", 
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
  D = 0.025, 
  cmax = 15,
  Criterion = "AIC",
  Variables = NULL,
  pdf = NULL,
  Theta1 = NULL,
  Theta2 = NULL,
  K = NULL,
  ymin = NULL,
  ymax = NULL,
  b = 1.0,
  ar = 0.1,
  Restraints = "loose")
{
  digits <- getOption("digits"); options(digits = 15)

  message("REBMIX Version 2.5.0");
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

  if (!is.numeric(D)) {
    stop(sQuote("D"), " numeric is requested!", call. = FALSE)
  }

  if ((D < 0.0) || (D > 1.0)) {
    stop(sQuote("D"), " must be greater or equal than 0.0 and less or equal than 1.0!", call. = FALSE)
  }

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

  if (!is.numeric(b)) {
    stop(sQuote("b"), " numeric is requested!", call. = FALSE)
  }

  if ((b < 0.0) || (b > 1.0)) {
    stop(sQuote("b"), " must be greater or equal than 0.0 and less or equal than 1.0!", call. = FALSE)
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
  REBMIX$Variables <- Variables
  REBMIX$pdf <- pdf
  REBMIX$Theta1 <- Theta1
  REBMIX$Theta2 <- Theta2  
  REBMIX$summary <- NULL
  REBMIX$pos <- 1
  REBMIX$all.Imax <- list()
  REBMIX$all.c <- list()
  REBMIX$all.IC <- list()
  REBMIX$all.logL <- list()
  REBMIX$all.D <- list()

  for (i in 1:length(Preprocessing)) {
    for (j in 1:length(Criterion)) {
      output <- .REBMIX(Dataset = Dataset, 
        Preprocessing = Preprocessing[i], 
        D = D, 
        cmax = cmax,
        Criterion = Criterion[j],
        Variables = Variables,
        pdf = pdf,
        Theta1 = Theta1,
        Theta2 = Theta2,
        K = if (is.list(K)) K[[i]] else K,
        ymin = ymin,
        ymax = ymax,
        b = b,
        ar = ar,
        Restraints = Restraints)

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
        REBMIX$all.Imax[[length(REBMIX$all.Imax) + 1]] <- output$all.Imax[[k]] 
        REBMIX$all.c[[length(REBMIX$all.c) + 1]] <- output$all.c[[k]] 
        REBMIX$all.IC[[length(REBMIX$all.IC) + 1]] <- output$all.IC[[k]] 
        REBMIX$all.logL[[length(REBMIX$all.logL) + 1]] <- output$all.logL[[k]] 
        REBMIX$all.D[[length(REBMIX$all.D) + 1]] <- output$all.D[[k]] 
      }      
    }
  }
  
  REBMIX$pos <- which(as.numeric(REBMIX$summary[, "logL"]) == max(as.numeric(REBMIX$summary[, "logL"])))  
  
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("REBMIX"))])

  class(REBMIX) <- "REBMIX"
  
  return(REBMIX)
} ## REBMIX
