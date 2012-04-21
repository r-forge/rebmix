boot <- function(x,
  pos = 1,
  Bootstrap = "parametric",
  B = 100, 
  n = NULL,
  replace = TRUE, 
  prob = NULL, ...) 
UseMethod("boot")

boot.default <- function(x,
  pos = 1,
  Bootstrap = "parametric",
  B = 100, 
  n = NULL,
  replace = TRUE, 
  prob = NULL, ...)
{
  stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
} ## boot.default

boot.REBMIX <- function(x,
  pos = 1,
  Bootstrap = "parametric",
  B = 100, 
  n = NULL,
  replace = TRUE, 
  prob = NULL, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  if (missing(x) || (class(x) != "REBMIX")) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  if ((pos < 1) || (pos > nrow(x$summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x$summary), "!", call. = FALSE)
  }
  
  if (!is.character(Bootstrap)) {
    stop(sQuote("Bootstrap"), " character is requested!", call. = FALSE)
  }

  Bootstrap <- match.arg(Bootstrap, .rebmix.boot$Bootstrap, several.ok = FALSE)    
  
  if (!is.wholenumber(B)) {
    stop(sQuote("B"), " integer is requested!", call. = FALSE)
  }

  if (B < 1) {
    stop(sQuote("B"), " must be greater than 0!", call. = FALSE)
  }
  
  nmax <- nrow(as.matrix(x$Dataset[[which(names(x$Dataset) == x$summary[pos, "Dataset"])]]))
  
  if (is.null(n)) {
    n <- nmax
  }
  else {
    if (!is.wholenumber(n)) {
      stop(sQuote("n"), " integer is requested!", call. = FALSE)
    }
  
    if ((n < 1) || (n > nmax)) {
      stop(sQuote("n"), " must be greater than 0 and less or equal than ", nmax, "!", call. = FALSE)
    }
  }  
  
  if (Bootstrap == .rebmix.boot$Bootstrap[1]) {
    bsample <- RNGMIX(Dataset = paste("bsample_", 1:B, sep = ""),
      n = rbind(round(n * as.numeric(x$w[[pos]]))),
      Theta = x$Theta[[pos]], ...)
  }
  else
  if (Bootstrap == .rebmix.boot$Bootstrap[2]) {
    Dataset <- as.matrix(x$Dataset[[which(names(x$Dataset) == x$summary[pos, "Dataset"])]])
    
    bsample <- list()
    
    bsample$Dataset <- list()
    
    range <- NULL

    for (i in 1:B) {
      R1 <- sample.int(nmax, size = n, replace = replace, prob = prob)
      
      bsample$Dataset[[i]] <- as.data.frame(Dataset[R1, ], stringsAsFactors = FALSE)
      
      range <- apply(rbind(range, Dataset[R1, , drop = FALSE]), 2, range)
    }
    
    bsample$ymin <- range[1, ]; bsample$ymax <- range[2, ]
  }
    
  d <- length(x$pdf)
  
  C <- x$summary[pos, "Preprocessing"]
  
  if (C == .rebmix$Preprocessing[1]) {
    y0 <- as.numeric(x$summary[pos, paste("y0", if (d > 1) 1:d, sep = "")])
    h <- as.numeric(x$summary[pos, paste("h", if (d > 1) 1:d, sep = "")])    
    
    kmin <- ceiling((y0 - bsample$ymin) / h - 0.5) 
    kmax <- ceiling((bsample$ymax - y0) / h - 0.5)
  
    k <- max(kmin + kmax + 1)
  
    ymin <- y0 - h * (kmin + 0.5); ymax <- ymin + h * k
  }
  else
  if (C == .rebmix$Preprocessing[2]) {
    h <- as.numeric(x$summary[pos, paste("h", if (d > 1) 1:d, sep = "")])    

    k <- ceiling((bsample$ymax - bsample$ymin) / h)
    
    k <- max(k)
    
    ymin <- bsample$ymin; ymax <- ymin + h * k
  }
  else
  if (C == .rebmix$Preprocessing[3]) {
    k <- as.numeric(x$summary[pos, "v/k"])

    ymin <- bsample$ymin; ymax <- bsample$ymax
  }

  bsampleest <- REBMIX(Dataset = bsample$Dataset,
    Preprocessing = as.character(x$summary[pos, "Preprocessing"]),
    D = as.numeric(x$summary[pos, "D"]),
    cmax = as.numeric(x$summary[pos, "cmax"]),
    Criterion = as.character(x$summary[pos, "Criterion"]),
    Variables = as.character(x$Variables),
    pdf = as.character(x$pdf),
    Theta1 = if (is.null(x$Theta1)) NULL else as.numeric(x$Theta1),
    Theta2 = if (is.null(x$Theta2)) NULL else as.numeric(x$Theta2),
    K = k,
    ymin = if (is.null(ymin)) NULL else as.numeric(ymin),
    ymax = if (is.null(ymax)) NULL else as.numeric(ymax),
    b = as.numeric(x$summary[pos, "b"]),
    ar = as.numeric(x$summary[pos, "ar"]),
    Restraints = as.character(x$summary[pos, "Restraints"]))

  freq <- table(as.numeric(bsampleest$summary$c))
  c <- as.integer(names(freq)[which.max(freq)])
  
  w <- rbind(bsampleest$w[as.numeric(bsampleest$summary$c) == c])
  
  Theta <- rbind(bsampleest$Theta[as.numeric(bsampleest$summary$c) == c])
  
  output <- list()
  
  output$c <- as.data.frame(rbind(as.numeric(bsampleest$summary$c)), stringsAsFactors = FALSE)
  output$c.se <- sd(as.numeric(bsampleest$summary$c))
  output$c.cv <- output$c.se / mean(as.numeric(bsampleest$summary$c))
  
  rownames(output$c) <- "c"
  colnames(output$c) <- paste(1:B, sep = "")  
  
  output$c.mode <- c
  output$c.prob <- length(w) / B
  
  output$w <- NULL
  output$w.se <- NULL
  
  for (i in 1:length(w)) {
    output$w <- cbind(output$w, as.numeric(w[[i]]))
  }
  
  output$w.se <- apply(output$w, 1, sd)
  output$w.cv <- output$w.se / apply(output$w, 1, mean)
  
  output$w <- as.data.frame(rbind(output$w), stringsAsFactors = FALSE)
  
  rownames(output$w) <- paste("comp", if (c > 1) 1:c else "", sep = "")
  colnames(output$w) <- paste(which(bsampleest$summary$c == c), sep = "") 
  
  nrow <- nrow(Theta[[1]])
  ncol <- ncol(Theta[[1]])  
  
  R1 <- NULL; R2 <- NULL; D1 <- NULL; D2 <- NULL
  
  i <- 1; j <- 1
    
  while (i < nrow) {
    pdf <- match.arg(Theta[[1]][i, 1], .rebmix$pdf)

    if (pdf %in% .rebmix$pdf[.rebmix$pdf.nargs == 2]) {
      R1 <- c(R1, i + 1); R2 <- c(R2, i + 2); D1 <- c(D1, j); D2 <- c(D2, j)
      
      i <- i + 3; j <- j + 1
    }
    else
    if (pdf %in% .rebmix$pdf[.rebmix$pdf.nargs == 1]) {
      R1 <- c(R1, i + 1); D1 <- c(D1, j)
      
      i <- i + 2; j <- j + 1
    }
  }
  
  for (i in 1:length(R1)) {
    theta1 <- if (d > 1) paste("theta1.",  D1[i], sep = "") else "theta1"
    
    output[[theta1]] <- NULL

    for (j in 1:length(Theta)) {
      output[[theta1]] <- cbind(output[[theta1]], as.numeric(Theta[[j]][R1[i], ]))
    }
    
    se <- apply(output[[theta1]], 1, sd)

    output[[paste(theta1, ".se", sep = "")]] <- se

    cv <- se / apply(output[[theta1]], 1, mean)
    
    output[[paste(theta1, ".cv", sep = "")]] <- cv
    
    output[[theta1]] <- data.frame(rbind(output[[theta1]]), stringsAsFactors = FALSE)
    
    rownames(output[[theta1]]) <- paste("comp", if (c > 1) 1:c else "", sep = "")
    colnames(output[[theta1]]) <- paste(which(bsampleest$summary$c == c), sep = "") 
  }
  
  for (i in 1:length(R2)) {
    theta2 <- if (d > 1) paste("theta2.",  D2[i], sep = "") else "theta2"
    
    output[[theta2]] <- NULL

    for (j in 1:length(Theta)) {
      output[[theta2]] <- cbind(output[[theta2]], as.numeric(Theta[[j]][R2[i], ]))
    }
    
    se <- apply(output[[theta2]], 1, sd)

    output[[paste(theta2, ".se", sep = "")]] <- se

    cv <- se / apply(output[[theta2]], 1, mean)
    
    output[[paste(theta2, ".cv", sep = "")]] <- cv
    
    output[[theta2]] <- data.frame(rbind(output[[theta2]]), stringsAsFactors = FALSE)
    
    rownames(output[[theta2]]) <- paste("comp", if (c > 1) 1:c else "", sep = "")
    colnames(output[[theta2]]) <- paste(which(bsampleest$summary$c == c), sep = "") 
  }  
  
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("output"))]) 
  
  class(output) <- "boot.REBMIX" 
   
  return(output)  
} ## boot.REBMIX
