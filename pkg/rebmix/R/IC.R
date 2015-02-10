.IC <- function(object, 
  Criterion = "AIC",
  pos = 1, ...)
{
  if (missing(object) || (class(object) != "REBMIX")) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }
  
  if (!is.character(Criterion)) {
    stop(sQuote("Criterion"), " character vector is requested!", call. = FALSE)
  }

  Criterion <- match.arg(Criterion, .rebmix$Criterion, several.ok = TRUE)  

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  if ((pos < 1) || (pos > nrow(object$summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(object$summary), "!", call. = FALSE)
  }

  Dataset <- as.character(object$summary[pos, "Dataset"])

  X <- as.matrix(object$Dataset[[which(names(object$Dataset) == object$summary[pos, "Dataset"])]])

  n <- nrow(X)
  d <- ncol(X)

  h <- as.double(object$summary[pos, paste("h", if (d > 1) 1:d else "", sep = "")])
  
  nrow <- nrow(object$Theta[[pos]])
  ncol <- ncol(object$Theta[[pos]])

  c <- ncol

  pdf <- array(data = NA, dim = c(nrow, ncol), dimnames = NULL)
  theta1 <- array(data = 0.0, dim = c(nrow, ncol), dimnames = NULL)
  theta2 <- array(data = 0.0, dim = c(nrow, ncol), dimnames = NULL)

  for (j in 1:ncol) {
    M <- match(object$Theta[[pos]][, j], .rebmix$pdf)

    d <- 1;

    for (l in 1:length(M)) {
      if (M[l] %in% which(.rebmix$pdf.nargs == 2)) {
        pdf[d, j] <- object$Theta[[pos]][l, j]
        theta1[d, j] <- as.numeric(object$Theta[[pos]][l + 1, j])
        theta2[d, j] <- as.numeric(object$Theta[[pos]][l + 2, j])

        d <- d + 1
      }
      else
      if (M[l] %in% which(.rebmix$pdf.nargs == 1)) {
        pdf[d, j] <- object$Theta[[pos]][l, j]
        theta1[d, j] <- as.numeric(object$Theta[[pos]][l + 1, j])

        d <- d + 1
      }
    }
  }

  d <- d - 1

  pdf <- pdf[1:d, ]; dim(pdf) <- c(d, ncol)
  
  theta1 <- theta1[1:d, ]; dim(theta1) <- c(d, ncol)
  theta2 <- theta2[1:d, ]; dim(theta2) <- c(d, ncol)

  C <- object$summary[pos, "Preprocessing"]

  if (C == .rebmix$Preprocessing[1]) {
    y0 <- as.double(object$summary[pos, paste("y0", if (d > 1) 1:d else "", sep = "")])

    output <- .C("RInformationCriterionH",
      h = as.double(h),
      y0 = as.double(y0),
      IniParFamType = as.character(object$call$pdf),
      k = as.integer(object$summary[pos, "v/k"]),
      n = as.integer(n),
      d = as.integer(d),
      object = as.double(X),
      ICType = as.character(Criterion),
      c = as.integer(c),
      W = as.double(object$w[[pos]]),
      ParFamType = as.character(pdf),
      Par0 = as.double(theta1),
      Par1 = as.double(theta2),      
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
      n = as.integer(n),
      d = as.integer(d),
      object = as.double(X),
      ICType = as.character(Criterion),
      c = as.integer(c),
      W = as.double(object$w[[pos]]),
      ParFamType = as.character(pdf),
      Par0 = as.double(theta1),
      Par1 = as.double(theta2),      
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
    k <- as.integer(object$summary[pos, "v/k"]) 

    output <- .C("RInformationCriterionKNN",
      k = as.integer(object$summary[pos, "v/k"]),
      h = as.double(h),
      n = as.integer(n),
      d = as.integer(d),
      object = as.double(X),
      ICType = as.character(Criterion),
      c = as.integer(c),
      W = as.double(object$w[[pos]]),
      ParFamType = as.character(pdf),
      Par0 = as.double(theta1),
      Par1 = as.double(theta2),      
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

logL <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "AIC", pos = pos, ...)

  output <- output$logL
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## logL

AIC <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "AIC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## AIC

AIC3 <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "AIC3", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## AIC3

AIC4 <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "AIC4", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## AIC4

AICc <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "AICc", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## AICc

BIC <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "BIC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## BIC

CAIC <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "CAIC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## CAIC

HQC <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "HQC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## HQC

MDL2 <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "MDL2", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## MDL2

MDL5 <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "MDL5", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## MDL5

AWE <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "AWE", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## AWE

CLC <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "CLC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## CLC

ICL <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "ICL", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## ICL

ICLBIC <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "ICL-BIC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## ICLBIC

PRD <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "D", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## PRD

SSE <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "SSE", pos = pos, ...)

  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## SSE

PC <- function(object, 
  pos = 1, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  output <- .IC(object = object, Criterion = "PC", pos = pos, ...)
  
  output <- output$IC
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
} ## PC




