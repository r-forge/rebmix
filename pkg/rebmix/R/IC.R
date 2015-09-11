.IC <- function(x, 
  Criterion = "AIC",
  pos = 1, ...)
{
  if (missing(x) || (class(x) != "REBMIX")) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  if (!is.character(Criterion)) {
    stop(sQuote("Criterion"), " character vector is requested!", call. = FALSE)
  }

  Criterion <- match.arg(Criterion, .rebmix$Criterion, several.ok = TRUE)  

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  if ((pos < 1) || (pos > nrow(x$summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x$summary), "!", call. = FALSE)
  }

  Dataset <- as.character(x$summary[pos, "Dataset"])

  X <- as.matrix(x$Dataset[[which(names(x$Dataset) == x$summary[pos, "Dataset"])]])

  n <- nrow(X)
  d <- ncol(X)

  h <- as.double(x$summary[pos, paste("h", if (d > 1) 1:d else "", sep = "")])
  
  Names <- names(x$Theta[[pos]])
  
  c <- length(x$w[[pos]])

  pdf <- as.character(unlist(x$Theta[[pos]][grep("pdf", Names)]))
  
  pdf <- match.arg(pdf, .rebmix$pdf)
  
  theta1 <- as.numeric(unlist(x$Theta[[pos]][grep("theta1", Names)]))
  
  theta1[is.na(theta1)] <- 0

  theta2 <- as.numeric(unlist(x$Theta[[pos]][grep("theta2", Names)]))
  
  theta2[is.na(theta2)] <- 0

  length(pdf) <- d

  C <- x$summary[pos, "Preprocessing"]

  if (C == .rebmix$Preprocessing[1]) {
    y0 <- as.double(x$summary[pos, paste("y0", if (d > 1) 1:d else "", sep = "")])

    output <- .C("RInformationCriterionH",
      h = as.double(h),
      y0 = as.double(y0),
      k = as.integer(x$summary[pos, "v/k"]),
      d = as.integer(d),
      Criterion = as.character(Criterion),
      c = as.integer(c),
      W = as.double(x$w[[pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(2),
      length.theta = as.integer(c(d, d)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(X),
      IC = double(1),
      logL = double(1),
      M = integer(1),
      D = double(1),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in IC!", call. = FALSE); return(NA)
    }
  } 
  else 
  if (C == .rebmix$Preprocessing[2]) {
    output <- .C("RInformationCriterionPW",
      h = as.double(h),
      d = as.integer(d),
      Criterion = as.character(Criterion),
      c = as.integer(c),
      W = as.double(x$w[[pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(2),
      length.theta = as.integer(c(d, d)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(X),
      IC = double(1),
      logL = double(1),
      M = integer(1),
      D = double(1),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in IC!", call. = FALSE); return(NA)
    }
  } 
  else
  if (C == .rebmix$Preprocessing[3]) {
    k <- as.integer(x$summary[pos, "v/k"]) 

    output <- .C("RInformationCriterionKNN",
      h = as.double(h),
      k = as.integer(x$summary[pos, "v/k"]),
      d = as.integer(d),
      Criterion = as.character(Criterion),
      c = as.integer(c),
      W = as.double(x$w[[pos]]),
      length.pdf = as.integer(d),
      length.Theta = as.integer(2),
      length.theta = as.integer(c(d, d)),
      pdf = as.character(pdf),
      Theta = as.double(c(theta1, theta2)),
      n = as.integer(n),
      x = as.double(X),
      IC = double(1),
      logL = double(1),
      M = integer(1),
      D = double(1),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in IC!", call. = FALSE); return(NA)
    }
  }
  
  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)
} ## .IC

logL <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "AIC", pos = pos, ...)

  output <- output$logL
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## logL

AIC <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "AIC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## AIC

AIC3 <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "AIC3", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## AIC3

AIC4 <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "AIC4", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## AIC4

AICc <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "AICc", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## AICc

BIC <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "BIC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## BIC

CAIC <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "CAIC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## CAIC

HQC <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "HQC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## HQC

MDL2 <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "MDL2", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## MDL2

MDL5 <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "MDL5", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## MDL5

AWE <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "AWE", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## AWE

CLC <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "CLC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## CLC

ICL <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "ICL", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## ICL

ICLBIC <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "ICL-BIC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## ICLBIC

PRD <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "D", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## PRD

SSE <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "SSE", pos = pos, ...)

  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## SSE

PC <- function(x, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(x = x, Criterion = "PC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## PC




