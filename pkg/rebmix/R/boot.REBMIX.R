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
      n = rbind(ceiling(n * as.numeric(x$w[[pos]]))),
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
    
  d <- length(x$call$pdf)
  
  bsampleest <- REBMIX(Dataset = bsample$Dataset,
    Preprocessing = as.character(x$summary[pos, "Preprocessing"]),
    cmax = as.numeric(x$summary[pos, "cmax"]),
    Criterion = as.character(x$summary[pos, "Criterion"]),
    Variables = x$call$Variables,
    pdf = x$call$pdf,
    theta1 = x$call$theta1,
    theta2 = x$call$theta2,
    K = eval(parse(text = as.character(x$summary[pos, "K"]))),
    y0 = x$call$y0,    
    ymin = x$call$ymin,
    ymax = x$call$ymax,
    ar = as.numeric(x$summary[pos, "ar"]),
    Restraints = as.character(x$summary[pos, "Restraints"]))

  freq <- table(as.numeric(bsampleest$summary$c))
  c <- as.integer(names(freq)[which.max(freq)])
  
  w <- bsampleest$w[as.numeric(bsampleest$summary$c) == c]
  
  Theta <- bsampleest$Theta[as.numeric(bsampleest$summary$c) == c]
  
  output <- list()
  
  output$c <- as.numeric(bsampleest$summary$c)
  output$c.se <- sd(as.numeric(bsampleest$summary$c))
  output$c.cv <- output$c.se / mean(as.numeric(bsampleest$summary$c))
  
  output$c.mode <- c
  output$c.prob <- length(w) / B
  
  output$w <- unlist(w); dim(output$w) <- c(length(w), c)
  
  colnames(output$w) <- paste("comp", if (c > 1) 1:c else "", sep = "")
  rownames(output$w) <- paste(which(bsampleest$summary$c == c), sep = "")

  output$w.se <- apply(output$w, 2, sd)
  output$w.cv <- output$w.se / apply(output$w, 2, mean)

  for (i in 1:output$c.mode) {
    theta1 <- paste("theta1.",  i, sep = "")

    output[[theta1]] <- NULL

    for (j in 1:length(Theta)) {
      output[[theta1]] <- c(output[[theta1]], Theta[[j]][[theta1]])
    }

    dim(output[[theta1]]) <- c(length(Theta), d)

    colnames(output[[theta1]]) <- x$call$pdf
    rownames(output[[theta1]]) <- paste(which(bsampleest$summary$c == c), sep = "")

    theta1.se <- paste("theta1.",  i, ".se", sep = "")
    theta1.cv <- paste("theta1.",  i, ".cv", sep = "")

    output[[theta1.se]] <- apply(output[[theta1]], 2, sd)
    output[[theta1.cv]] <- output[[theta1.se]] / apply(output[[theta1]], 2, mean)
  }

  for (i in 1:output$c.mode) {
    theta2 <- paste("theta2.",  i, sep = "")

    output[[theta2]] <- NULL

    for (j in 1:length(Theta)) {
      output[[theta2]] <- c(output[[theta2]], Theta[[j]][[theta2]])
    }

    dim(output[[theta2]]) <- c(length(Theta), d)

    colnames(output[[theta2]]) <- x$call$pdf
    rownames(output[[theta2]]) <- paste(which(bsampleest$summary$c == c), sep = "")

    theta2.se <- paste("theta2.",  i, ".se", sep = "")
    theta2.cv <- paste("theta2.",  i, ".cv", sep = "")

    output[[theta2.se]] <- apply(output[[theta2]], 2, sd)
    output[[theta2.cv]] <- output[[theta2.se]] / apply(output[[theta2]], 2, mean)
  }
  
  options(digits = digits)  

  rm(list = ls()[!(ls() %in% c("output"))]) 
  
  class(output) <- "boot.REBMIX" 
   
  return(output)  
} ## boot.REBMIX
