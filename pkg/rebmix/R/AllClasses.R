### Panic Branislav.
setClass(Class = "EM.Control",
slots = c(strategy = "character",
  variant = "character",
  acceleration = "character",
  tolerance = "numeric",
  acceleration.multiplier = "numeric",
  maximum.iterations = "numeric",
  K = "numeric",
### Panic Branislav.
  eliminate.zero.components = "logical"))
### End  

setMethod("initialize", "EM.Control",
function(.Object, ...,
  strategy,
  variant,
  acceleration,
  tolerance,
  acceleration.multiplier,
  maximum.iterations,
  K,
### Panic Branislav.
  eliminate.zero.components)
### End 
{
  # strategy.

  if (missing(strategy) || length(strategy) == 0) {
    strategy <- .rebmix$EMStrategy[1]
  } 
  else {
    strategy <- match.arg(strategy, .rebmix$EMStrategy)
  }

  # variant.

  if (missing(variant) || length(variant) == 0) {
    variant <- .rebmix$EMVariant[1]
  } 
  else {
    variant <- match.arg(variant, .rebmix$EMVariant)
  }

  # acceleration.

  if (missing(acceleration) || length(acceleration) == 0) {
    acceleration <- .rebmix$EMAcceleration[1]
  } 
  else {
    acceleration <- match.arg(acceleration, .rebmix$EMAcceleration)
  }

  # tolerance.
  
  if (missing(tolerance) || (length(tolerance) == 0)) {
    tolerance <- 1.E-4
  }

  if (!is.numeric(tolerance)) {
    stop(sQuote("tolerance"), " numeric is requested!", call. = FALSE)
  }

  length(tolerance) <- 1

  if (tolerance <= 0.0) {
    stop(sQuote("tolerance"), " must be greater than 0.0!", call. = FALSE)
  }  

  # acceleration.multiplier.
  
  if (missing(acceleration.multiplier) || (length(acceleration.multiplier) == 0)) {
    acceleration.multiplier <- 1.0
  }

  if (!is.numeric(acceleration.multiplier)) {
    stop(sQuote("acceleration.multiplier"), " numeric is requested!", call. = FALSE)
  }

  length(acceleration.multiplier) <- 1

  if (acceleration.multiplier < 1.0 || acceleration.multiplier > 2.0) {
    stop(sQuote("acceleration.multiplier"), " must be greater or equal than 1.0 and less or equal than 2.0!", call. = FALSE)
  }  

  # maximum.iterations.
  
  if (missing(maximum.iterations) || (length(maximum.iterations) == 0)) {
    maximum.iterations <- as.integer(1000)
  }

  if (!is.wholenumber(maximum.iterations)) {
    stop(sQuote("maximum.iterations"), " integer is requested!", call. = FALSE)
  }

  length(maximum.iterations) <- 1

  if (maximum.iterations < 1) {
    stop(sQuote("maximum.iterations"), " must be greater than 0!", call. = FALSE)
  }  
  
  # K.
  
  if (missing(K) || (length(K) == 0)) {
    K <- as.integer(0)
  }

  if (!is.wholenumber(K)) {
    stop(sQuote("K"), " integer is requested!", call. = FALSE)
  }

  length(K) <- 1

  if (K < 0) {
    stop(sQuote("K"), " must be greater or equal than 0!", call. = FALSE)
  }    

  # eliminate.zero.components.

  if (missing(eliminate.zero.components) || (length(eliminate.zero.components) == 0)){
    eliminate.zero.components <- FALSE;
  }

  length(eliminate.zero.components) <- 1

  if (!is.logical(eliminate.zero.components)){
    stop(sQuote("eliminate.zero.components"), " logical is requested!", call. = FALSE)
  } 

  .Object@strategy <- strategy
  .Object@variant <- variant
  .Object@acceleration <- acceleration
  .Object@tolerance <- tolerance
  .Object@acceleration.multiplier <- acceleration.multiplier
  .Object@maximum.iterations <- maximum.iterations
  .Object@K <- K
  .Object@eliminate.zero.components <- eliminate.zero.components  

  rm(list = ls()[!(ls() %in% c(".Object"))])

  .Object
}) ## initialize

setMethod("show", 
  signature(object = "EM.Control"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class EM.Control is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"strategy\":", "\n", sep = "")

  print(object@strategy, quote = FALSE)

  cat("Slot \"variant\":", "\n", sep = "")

  print(object@variant, quote = FALSE)
    
  cat("Slot \"acceleration\":", "\n", sep = "")

  print(object@acceleration, quote = FALSE)

  cat("Slot \"tolerance\":", "\n", sep = "")

  print(object@tolerance, quote = FALSE)

  cat("Slot \"acceleration.multiplier\":", "\n", sep = "")

  print(object@acceleration.multiplier, quote = FALSE)

  cat("Slot \"maximum.iterations\":", "\n", sep = "")

  print(object@maximum.iterations, quote = FALSE)
  
  cat("Slot \"K\":", "\n", sep = "")

  print(object@K, quote = FALSE)  

  rm(list = ls())
}) ## show
### End

# Class RNGMIX.Theta

setClass(Class = "RNGMIX.Theta",
slots = c(c = "numeric",
  d = "numeric",
  pdf = "character",
  Theta = "list"))

setMethod("initialize", "RNGMIX.Theta",
function(.Object, ...,
  c,
  pdf)
{
  # c.

  if (missing(c) || (length(c) == 0)) {
    stop(sQuote("c"), " must not be empty!", call. = FALSE)
  }

  if (!is.wholenumber(c)) {
    stop(sQuote("c"), " integer is requested!", call. = FALSE)
  }

  length(c) <- 1

  if (c < 1) {
    stop(sQuote("c"), " must be greater than 0!", call. = FALSE)
  }

  # pdf.

  if (missing(pdf) || (length(pdf) == 0)) {
    stop(sQuote("pdf"), " must not be empty!", call. = FALSE)
  }

  if (!is.character(pdf)) {
    stop(sQuote("pdf"), " character vector is requested!", call. = FALSE)
  }

  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)

  # d.

  d <- length(pdf)

  # Theta.

  Theta <- list()

  length(Theta) <- 4 * c

  names(Theta)[seq(1, 4 * c, 4)] <- paste("pdf", 1:c, sep = "")
  names(Theta)[seq(2, 4 * c, 4)] <- paste("theta1.", 1:c, sep = "")
  names(Theta)[seq(3, 4 * c, 4)] <- paste("theta2.", 1:c, sep = "")
  names(Theta)[seq(4, 4 * c, 4)] <- paste("theta3.", 1:c, sep = "")  

  M1 <- which(pdf %in% .rebmix$pdf[.rebmix$pdf.nargs < 2])
  M2 <- which(pdf %in% .rebmix$pdf[.rebmix$pdf.nargs < 3])
  M3 <- which(pdf %in% .rebmix$pdf[8])

  for (i in 1:c) {
    Theta[[1 + (i - 1) * 4]] <- pdf
    Theta[[2 + (i - 1) * 4]] <- array(data = 0.0, dim = d)
    Theta[[3 + (i - 1) * 4]] <- array(data = 0.0, dim = d)
    Theta[[4 + (i - 1) * 4]] <- array(data = 0.0, dim = d)

    Theta[[3 + (i - 1) * 4]][M1] <- NA
    Theta[[4 + (i - 1) * 4]][M2] <- NA

    Theta[[2 + (i - 1) * 4]][M3] <- -Inf
    Theta[[3 + (i - 1) * 4]][M3] <- Inf
  }

  .Object@c <- c
  .Object@d <- d
  .Object@pdf <- pdf
  .Object@Theta <- Theta

  rm(list = ls()[!(ls() %in% c(".Object"))])

  .Object
}) ## initialize

setMethod("show",
          signature(object = "RNGMIX.Theta"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RNGMIX.Theta is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"c\":", "\n", sep = "")

  print(object@c, quote = FALSE)

  cat("Slot \"d\":", "\n", sep = "")

  print(object@d, quote = FALSE)

  cat("Slot \"pdf\":", "\n", sep = "")

  print(object@pdf, quote = FALSE)

  cat("Slot \"Theta\":", "\n", sep = "")

  print(object@Theta, quote = FALSE)

  rm(list = ls())
}) ## show

## Class RNGMVNORM.Theta

setClass("RNGMVNORM.Theta", contains = "RNGMIX.Theta")

setMethod("initialize", "RNGMVNORM.Theta",
function(.Object, ...,
  c,
  d)
{
  # c.

  if (missing(c) || (length(c) == 0)) {
    stop(sQuote("c"), " must not be empty!", call. = FALSE)
  }

  if (!is.wholenumber(c)) {
    stop(sQuote("c"), " integer is requested!", call. = FALSE)
  }

  length(c) <- 1

  if (c < 1) {
    stop(sQuote("c"), " must be greater than 0!", call. = FALSE)
  }

  # d.

  if (missing(d) || (length(d) == 0)) {
    stop(sQuote("d"), " must not be empty!", call. = FALSE)
  }

  if (!is.wholenumber(d)) {
    stop(sQuote("d"), " integer is requested!", call. = FALSE)
  }

  length(d) <- 1

  if (d < 1) {
    stop(sQuote("d"), " must be greater than 0!", call. = FALSE)
  }

  # pdf.

  pdf <- rep(.rebmix$pdf[1], d)

  # Theta.

  Theta <- list()

  length(Theta) <- 3 * c

  names(Theta)[seq(1, 3 * c, 3)] <- paste("pdf", 1:c, sep = "")
  names(Theta)[seq(2, 3 * c, 3)] <- paste("theta1.", 1:c, sep = "")
  names(Theta)[seq(3, 3 * c, 3)] <- paste("theta2.", 1:c, sep = "")

  for (i in 1:c) {
    Theta[[1 + (i - 1) * 3]] <- pdf
    Theta[[2 + (i - 1) * 3]] <- array(data = 0.0, dim = d)
    Theta[[3 + (i - 1) * 3]] <- array(data = 0.0, dim = d * d)
  }

  .Object@c <- c
  .Object@d <- d
  .Object@pdf <- pdf
  .Object@Theta <- Theta

  rm(list = ls()[!(ls() %in% c(".Object"))])

  .Object
}) ## initialize

### Panic Branislav.
# Class EMMIX.Theta

setClass("EMMIX.Theta", 
  slots = c(w = "numeric"),
  contains = "RNGMIX.Theta")

setMethod("initialize", "EMMIX.Theta",
function(.Object, ...,
  c,
  d,
  pdf)
{
  # c.

  if (missing(c) || (length(c) == 0)) {
    stop(sQuote("c"), " must not be empty!", call. = FALSE)
  }

  if (!is.wholenumber(c)) {
    stop(sQuote("c"), " integer is requested!", call. = FALSE)
  }

  length(c) <- 1

  if (c < 1) {
    stop(sQuote("c"), " must be greater than 0!", call. = FALSE)
  }

  # pdf.

  if (missing(pdf) || (length(pdf) == 0)) {
    stop(sQuote("pdf"), " must not be empty!", call. = FALSE)
  }

  if (!is.character(pdf)) {
    stop(sQuote("pdf"), " character vector is requested!", call. = FALSE)
  }

  pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)

  # d.

  if (missing(d) || (length(d) == 0)) {
    d <- length(pdf)
  }else {
    if (!is.wholenumber(d)) {
       stop(sQuote("d"), " integer is requested!", call. = FALSE)
    }

    if (d != length(pdf)){
      stop(sQuote("d"), " must be equal to ", length(pdf), " !", call. = FALSE)
    }
  }

  length(d) <- 1

  # w.

  w <- double(c)

  # Theta.

  Theta <- list()

  length(Theta) <- 4 * c

  names(Theta)[seq(1, 4 * c, 4)] <- paste("pdf", 1:c, sep = "")
  names(Theta)[seq(2, 4 * c, 4)] <- paste("theta1.", 1:c, sep = "")
  names(Theta)[seq(3, 4 * c, 4)] <- paste("theta2.", 1:c, sep = "")
  names(Theta)[seq(4, 4 * c, 4)] <- paste("theta3.", 1:c, sep = "") 

  M1 <- which(pdf %in% .rebmix$pdf[.rebmix$pdf.nargs < 2])
  M2 <- which(pdf %in% .rebmix$pdf[.rebmix$pdf.nargs < 3])

  for (i in 1:c) {
    Theta[[1 + (i - 1) * 4]] <- pdf
    Theta[[2 + (i - 1) * 4]] <- array(data = 0.0, dim = d)
    Theta[[3 + (i - 1) * 4]] <- array(data = 0.0, dim = d)
    Theta[[4 + (i - 1) * 4]] <- array(data = 0.0, dim = d)

    Theta[[3 + (i - 1) * 4]][M1] <- NA
    Theta[[4 + (i - 1) * 4]][M2] <- NA
  }

  .Object@c <- c
  .Object@d <- d
  .Object@pdf <- pdf
  .Object@w <- w
  .Object@Theta <- Theta

  rm(list = ls()[!(ls() %in% c(".Object"))])

  .Object
}) ## initialize

setMethod("show",
          signature(object = "EMMIX.Theta"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class EMMIX.Theta is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"c\":", "\n", sep = "")

  print(object@c, quote = FALSE)

  cat("Slot \"d\":", "\n", sep = "")

  print(object@d, quote = FALSE)

  cat("Slot \"pdf\":", "\n", sep = "")

  print(object@pdf, quote = FALSE)

  cat("Slot \"w\":", "\n", sep = "")

  print(object@w, quote = FALSE)

  cat("Slot \"Theta\":", "\n", sep = "")

  print(object@Theta, quote = FALSE)

  rm(list = ls())
}) ## show

# Class EMMVNORM.Theta

setClass("EMMVNORM.Theta", 
  slots = c(w = "numeric"),
  contains = "RNGMVNORM.Theta")

setMethod("initialize", "EMMVNORM.Theta",
function(.Object, ...,
  c,
  d,
  pdf)
{
  # c.

  if (missing(c) || (length(c) == 0)) {
    stop(sQuote("c"), " must not be empty!", call. = FALSE)
  }

  if (!is.wholenumber(c)) {
    stop(sQuote("c"), " integer is requested!", call. = FALSE)
  }

  length(c) <- 1

  if (c < 1) {
    stop(sQuote("c"), " must be greater than 0!", call. = FALSE)
  }

  # d.

  if (missing(d) || (length(d) == 0)) {
    stop(sQuote("d"), " must not be empty!", call. = FALSE)
  }

  if (!is.wholenumber(d)) {
    stop(sQuote("d"), " integer is requested!", call. = FALSE)
  }

  length(d) <- 1

  if (d < 1) {
    stop(sQuote("d"), " must be greater than 0!", call. = FALSE)
  }

  # pdf.

  if (missing(pdf) || (length(pdf) == 0)) {
    pdf <- rep(.rebmix$pdf[1], d)
  }
  else {
    pdf <- match.arg(pdf, .rebmix$pdf[1], several.ok = TRUE)

    if (length(pdf) != d){
      stop(sQuote("pdf"), " Length must be equal to the ", sQuote("d"), " parameter!", call. = FALSE)
    }
  }

  # w.

  w <- double(c)

  # Theta.

  Theta <- list()

  length(Theta) <- 3 * c

  names(Theta)[seq(1, 3 * c, 3)] <- paste("pdf", 1:c, sep = "")
  names(Theta)[seq(2, 3 * c, 3)] <- paste("theta1.", 1:c, sep = "")
  names(Theta)[seq(3, 3 * c, 3)] <- paste("theta2.", 1:c, sep = "")

  for (i in 1:c) {
    Theta[[1 + (i - 1) * 3]] <- pdf
    Theta[[2 + (i - 1) * 3]] <- array(data = 0.0, dim = d)
    Theta[[3 + (i - 1) * 3]] <- array(data = 0.0, dim = d * d)
  }

  .Object@c <- c
  .Object@d <- d
  .Object@pdf <- pdf
  .Object@w <- w
  .Object@Theta <- Theta

  rm(list = ls()[!(ls() %in% c(".Object"))])

  .Object
}) ## initialize

setMethod("show",
          signature(object = "EMMVNORM.Theta"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class EMMVNORM.Theta is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"c\":", "\n", sep = "")

  print(object@c, quote = FALSE)

  cat("Slot \"d\":", "\n", sep = "")

  print(object@d, quote = FALSE)

  cat("Slot \"pdf\":", "\n", sep = "")

  print(object@pdf, quote = FALSE)

  cat("Slot \"w\":", "\n", sep = "")

  print(object@w, quote = FALSE)

  cat("Slot \"Theta\":", "\n", sep = "")

  print(object@Theta, quote = FALSE)

  rm(list = ls())
}) ## show
### End

# Class RNGMIX

setClass(Class = "RNGMIX",
slots = c(Dataset.name = "character",
  rseed = "numeric",
  n = "numeric",
  Theta = "list",
  Dataset = "list",
  Zt = "factor",
  w = "numeric",
  Variables = "character",
  ymin = "numeric",
  ymax = "numeric"),
prototype = list(rseed = -1))

setMethod("initialize", "RNGMIX",
function(.Object, ...,
  Dataset.name,
  rseed,
  n,
  Theta)
{
  # Dataset.name.

  if (missing(Dataset.name) || (length(Dataset.name) == 0)) {
    stop(sQuote("Dataset.name"), " must not be empty!", call. = FALSE)
  }

  if (!is.character(Dataset.name)) {
    stop(sQuote("Dataset.name"), " character vector is requested!", call. = FALSE)
  }

  # rseed.

  if (missing(rseed) || (length(rseed) == 0)) rseed <- .Object@rseed

  if (!is.wholenumber(rseed)) {
    stop(sQuote("rseed"), " integer is requested!", call. = FALSE)
  }

  length(rseed) <- 1

  if (rseed > -1) {
    stop(sQuote("rseed"), " must be less than 0!", call. = FALSE)
  }

  # n.

  if (missing(n) || (length(n) == 0)) {
    stop(sQuote("n"), " must not be empty!", call. = FALSE)
  }

  if (!is.wholenumber(n)) {
    stop(sQuote("n"), " integer vector is requested!", call. = FALSE)
  }

  if (!all(n > 0)) {
    stop("all ", sQuote("n"), " must be greater than 0!", call. = FALSE)
  }

  c <- length(n)

  # Theta.

  if (missing(Theta) || (length(Theta) == 0)) {
    stop(sQuote("Theta"), " must not be empty!", call. = FALSE)
  }

  if (!is.list(Theta)) {
    stop(sQuote("Theta"), " list is requested!", call. = FALSE)
  }

  Names <- names(Theta)

  if (length(grep("pdf", Names)) == 0) {
    stop(sQuote("pdf1"), " in " , sQuote("Theta"), " character vector is requested!", call. = FALSE)
  }

  length(grep("pdf", Names))

  j <- 0; length.pdf <- length(Theta[[grep("pdf", Names)[1]]])

  for (i in grep("pdf", Names)) {
    pdf <- as.character(Theta[[i]])

    pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)

    if (length(pdf) != length.pdf) {
      stop("lengths of ", sQuote("pdfi"), " in " , sQuote("Theta"), " must be equal!", call. = FALSE)
    }

    Theta[[i]] <- pdf; j <- j + 1
  }
  
  if ((length.pdf > 1) && (j != c)) {
    stop(sQuote("pdfi"), " in " , sQuote("Theta"), " and ", sQuote("n"), " must match!", call. = FALSE)
  }
  
  if (length(grep("theta1", Names)) == 0) {
    stop(sQuote("theta1.1"), " in " , sQuote("Theta"), " numeric vector is requested!", call. = FALSE)
  }

  j <- 0; length.theta1 <- length(Theta[[grep("theta1", Names)[1]]])

  for (i in grep("theta1", Names)) {
    theta1 <- as.numeric(Theta[[i]])

    if (length(theta1) != length.theta1) {
      stop("lengths of ", sQuote("theta1.l"), " in " , sQuote("Theta"), " must be equal!", call. = FALSE)
    }

    j <- j + 1
  }

  if ((length.pdf > 1) && (j != c)) {
    stop(sQuote("theta1.l"), " in " , sQuote("Theta"), " and ", sQuote("n"), " must match!", call. = FALSE)
  }
  
  if (length(grep("theta2", Names)) == 0) {
    stop(sQuote("theta2.1"), " in " , sQuote("Theta"), " numeric vector is requested!", call. = FALSE)
  }
  
  j <- 0; length.theta2 <- length(Theta[[grep("theta2", Names)[1]]])

  for (i in grep("theta2", Names)) {
    theta2 <- as.numeric(Theta[[i]])

    if (length(theta2) != length.theta2) {
      stop("lengths of ", sQuote("theta2.l"), " in " , sQuote("Theta"), " must be equal!", call. = FALSE)
    }

    j <- j + 1
  }

  if ((length.pdf > 1) && (j != c)) {
    stop(sQuote("theta2.l"), " in " , sQuote("Theta"), " and ", sQuote("n"), " must match!", call. = FALSE)
  }
  
  if (as.character(class(.Object)) == "RNGMIX") {
    if (length(grep("theta3", Names)) == 0) {
      stop(sQuote("theta3.1"), " in " , sQuote("Theta"), " numeric vector is requested!", call. = FALSE)
    }   
  
    j <- 0; length.theta3 <- length(Theta[[grep("theta3", Names)[1]]])

    for (i in grep("theta3", Names)) {
      theta3 <- as.numeric(Theta[[i]])

      if (length(theta3) != length.theta3) {
        stop("lengths of ", sQuote("theta3.l"), " in " , sQuote("Theta"), " must be equal!", call. = FALSE)
      }

      j <- j + 1
    }

    if ((length.pdf > 1) && (j != c)) {
      stop(sQuote("theta3.l"), " in " , sQuote("Theta"), " and ", sQuote("n"), " must match!", call. = FALSE)
    }
  }

  # Variables.

  for (i in 1:length(.rebmix$pdf)) {
    .Object@Variables[which(pdf == .rebmix$pdf[i])] <- .rebmix$pdf.Variables[i]
  }

  .Object@Dataset.name <- Dataset.name
  .Object@rseed <- rseed
  .Object@n <- n
  .Object@Theta <- Theta

  rm(list = ls()[!(ls() %in% c(".Object"))])

  .Object
}) ## initialize

setMethod("show",
          signature(object = "RNGMIX"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RNGMIX is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"w\":", "\n", sep = "")

  print(object@w, quote = FALSE)

  cat("Slot \"ymin\":", "\n", sep = "")

  print(object@ymin, quote = FALSE)

  cat("Slot \"ymax\":", "\n", sep = "")

  print(object@ymax, quote = FALSE)

  rm(list = ls())
}) ## show

# Class RNGMVNORM

setClass("RNGMVNORM", contains = "RNGMIX")

setMethod("show",
          signature(object = "RNGMVNORM"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RNGMVNORM is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"w\":", "\n", sep = "")

  print(object@w, quote = FALSE)

  cat("Slot \"ymin\":", "\n", sep = "")

  print(object@ymin, quote = FALSE)

  cat("Slot \"ymax\":", "\n", sep = "")

  print(object@ymax, quote = FALSE)

  rm(list = ls())
}) ## show

# Class REBMIX

setClass(Class = "REBMIX",
slots = c(Dataset = "list",
  Preprocessing = "character",
  cmax = "numeric",
  cmin = "numeric",
  Criterion = "character",
  Variables = "character",
  pdf = "character",
  theta1 = "numeric",
  theta2 = "numeric",
  theta3 = "numeric",  
  K = "ANY",
  ymin = "numeric",
  ymax = "numeric",
  ar = "numeric",
  Restraints = "character",
### Panic Branislav.  
  EMcontrol = "ANY",
### End    
  w = "list",
  Theta = "list",
  summary = "ANY",
### Panic Branislav.   
  summary.EM = "ANY",
### End   
  pos = "numeric",
  opt.c = "list",
  opt.IC = "list",
  opt.logL = "list",
  opt.D = "list",
  all.K = "list",
  all.IC = "list"),
prototype = list(cmax = 15,
  cmin = 1,
  Criterion = "AIC",
  ar = 0.1,
  Restraints = "loose",
  pos = 1))

setMethod("initialize", "REBMIX",
function(.Object, ...,
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
### Panic Branislav.  
  EMcontrol)
### End  
{
  # Dataset.

  if (missing(Dataset) || (length(Dataset) == 0)) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }

  if (!is.list(Dataset)) {
    stop(sQuote("Dataset"), " list is requested!", call. = FALSE)
  }

  if (is.null(names(Dataset))) {
    names(Dataset) <- paste("dataset", 1:length(Dataset), sep = "")
  }

  if (!all(unlist(lapply(Dataset, function(x) as.character(class(x)) == "Histogram" || as.character(class(x)) == "data.frame")))) {
    stop(sQuote("Dataset"), " list of data frames or objects of class ", "Histogram", " is requested!", call. = FALSE)
  }

  d <- NULL;
  
  for (i in 1:length(Dataset)) {
    if (as.character(class(Dataset[[i]])) == "data.frame") {
      di <- ncol(Dataset[[i]])
    
      d <- c(d, di)
    
      if (di < 1) {
        stop(sQuote("Dataset"), " numbers of columns in data frames must be greater than 0!", call. = FALSE)
      }

      Dataset[[i]] <- as.data.frame(Dataset[[i]][complete.cases(Dataset[[i]]), ])

      if (nrow(Dataset[[i]]) < 2) {
        stop(sQuote("Dataset"), " numbers of rows in data frames must be greater than 1!", call. = FALSE)
      }
      
      # Preprocessing.

      if (missing(Preprocessing) || (length(Preprocessing) == 0)) {
        stop(sQuote("Preprocessing"), " must not be empty!", call. = FALSE)
      }

      if (!is.character(Preprocessing)) {
        stop(sQuote("Preprocessing"), " character vector is requested!", call. = FALSE)
      }

      Preprocessing <- match.arg(Preprocessing, .rebmix$Preprocessing, several.ok = FALSE)      
    }
    else {
      di <- ncol(Dataset[[i]]@Y) - 1
      
      d <- c(d, di)
    }
  }
  
  d <- unique(d)

  if (length(d) != 1) {
    stop(sQuote("Dataset"), " numbers of variables in data frames must be equal!", call. = FALSE)
  }  

  # cmax.

  if (missing(cmax) || (length(cmax) == 0)) cmax <- .Object@cmax

  if (!is.wholenumber(cmax)) {
    stop(sQuote("cmax"), " integer is requested!", call. = FALSE)
  }

  if (cmax < 1) {
    stop(sQuote("cmax"), " must be greater than 0!", call. = FALSE)
  }

  # cmin.

  if (missing(cmin) || (length(cmin) == 0)) cmin <- .Object@cmin

  if (!is.wholenumber(cmin)) {
    stop(sQuote("cmin"), " integer is requested!", call. = FALSE)
  }

  if (cmin < 1) {
    stop(sQuote("cmin"), " must be greater than 0!", call. = FALSE)
  }

  if (cmin > cmax) {
    stop(sQuote("cmax"), " must be greater or equal than ", cmin, "!", call. = FALSE)
  }

  # Criterion.

  if (missing(Criterion) || (length(Criterion) == 0)) Criterion <- .Object@Criterion

  if (!is.character(Criterion)) {
    stop(sQuote("Criterion"), " character is requested!", call. = FALSE)
  }

  Criterion <- match.arg(Criterion, .rebmix$Criterion, several.ok = FALSE)

  # pdf.

  if (.Object@class[1] == "REBMVNORM") {
    pdf <- rep(.rebmix$pdf[1], d)
  }
  else {
    if (missing(pdf) || (length(pdf) == 0)) {
      stop(sQuote("pdf"), " must not be empty!", call. = FALSE)
    }

    if (!is.character(pdf)) {
      stop(sQuote("pdf"), " character vector is requested!", call. = FALSE)
    }

    pdf <- match.arg(pdf, .rebmix$pdf, several.ok = TRUE)

    if (length(pdf) != d) {
      stop("lengths of ", sQuote("pdf"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
  }

  # theta1.

  if (missing(theta1) || (length(theta1) == 0)) {
    if (any(pdf == .rebmix$pdf[4])) {
      stop(sQuote("theta1"), " must not be empty for ", dQuote(.rebmix$pdf[4]), "!", call. = FALSE)
    }
    else {
      theta1 <- .Object@theta1
    }
  }
  else {
    if (length(theta1) != d) {
      stop("lengths of ", sQuote("theta1"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
    
    for (i in 1:d) {
      if (pdf[i] == .rebmix$pdf[4]) {
        if (!is.wholenumber(theta1[i])) {
          stop(sQuote("theta1"), " integer is requested for ", dQuote(.rebmix$pdf[4]), "!", call. = FALSE)
        }

        if (theta1[i] < 0.0) {
          stop(sQuote("theta1"), " for ", dQuote(.rebmix$pdf[4]), " must be greater or equal than 0!", call. = FALSE)
        }
      }
    }

    class(theta1) <- "numeric"
  }

  # theta2.

  if (missing(theta2) || (length(theta2) == 0)) {
    theta2 <- .Object@theta2
  }
  else {
    if (length(theta2) != d) {
      stop("lengths of ", sQuote("theta2"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }

    class(theta2) <- "numeric"
  }
  
  # theta3.

  if (missing(theta3) || (length(theta3) == 0)) {
    if (any(pdf == .rebmix$pdf[10])) {
      stop(sQuote("theta3"), " must not be empty for ", dQuote(.rebmix$pdf[10]), "!", call. = FALSE)
    }
    else {
      theta3 <- .Object@theta3
    }
  }
  else {
    if (length(theta3) != d) {
      stop("lengths of ", sQuote("theta3"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
    
    for (i in 1:d) {
      if (pdf[i] == .rebmix$pdf[10]) {
        if (is.na(theta3[i])) {
        }
        else
        if (is.character(theta3[i])) {
          stop(sQuote("theta3"), " numeric is requested!", call. = FALSE)
        }
        else
        if (abs(abs(theta3[i]) - 1.0) > .Machine$double.eps^0.5) {
          stop(sQuote("theta3"), " for ", dQuote(.rebmix$pdf[10]), " must be -1.0, NA or 1.0!", call. = FALSE)
        }
      }
    }

    class(theta3) <- "numeric"
  }

  # K.

  if (missing(K) || (length(K) == 0)) {
    stop(sQuote("K"), " must not be empty!", call. = FALSE)
  }

  if (is.matrix(K)) {
    if (!all(unlist(lapply(K, is.wholenumber) == TRUE))) {
      stop(sQuote("K"), " matrix of integer values is requested!", call. = FALSE)
    }

    if (!all(unlist(lapply(K, function(x) all(x > 0)))) == TRUE) {
      stop("all ", sQuote("K"), " must be greater than 0!", call. = FALSE)
    }
    
    if(ncol(K) != d) {
      stop(sQuote("K"), " number of columns in matrix must equal ", d, "!", call. = FALSE)    
    }
    
    if(nrow(K) != length(Dataset)) {
      stop(sQuote("K"), " number of rows in matrix must equal ", length(Dataset), "!", call. = FALSE)    
    }
  }
  else
  if (is.numeric(K)) {
    if (!is.wholenumber(K)) {
      stop(sQuote("K"), " integer vector is requested!", call. = FALSE)
    }

    if (!all(K > 0)) {
      stop("all ", sQuote("K"), " must be greater than 0!", call. = FALSE)
    }
  }
  else {
    K <- "auto"
  }

  # ymin.

  if (missing(ymin) || (length(ymin) == 0)) {
    ymin <- .Object@ymin
  }
  else {
    if (length(ymin) != d) {
      stop("lengths of ", sQuote("ymin"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
  }

  # ymax.

  if (missing(ymax) || (length(ymax) == 0)) {
    ymax <- .Object@ymax
  }
  else {
    if (length(ymax) != d) {
      stop("lengths of ", sQuote("ymax"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
    
    if (any(ymax <= ymin)) {
      stop(sQuote("ymax"), " must be greater than ", sQuote("ymin"), "!", call. = FALSE)
    }    
  }

  # ar.

  if (missing(ar) || (length(ar) == 0)) ar <- .Object@ar

  if (!is.numeric(ar)) {
    stop(sQuote("ar"), " numeric is requested!", call. = FALSE)
  }

  length(ar) <- 1

  if ((ar <= 0.0) || (ar > 1.0)) {
    stop(sQuote("ar"), " must be greater than 0.0 and less or equal than 1.0!", call. = FALSE)
  }

  # Restraints.

  if (missing(Restraints) || (length(Restraints) == 0)) Restraints <- .Object@Restraints

  if (!is.character(Restraints)) {
    stop(sQuote("Restraints"), " character is requested!", call. = FALSE)
  }

  Restraints <- match.arg(Restraints, .rebmix$Restraints, several.ok = FALSE)

  # Variables.

  for (i in 1:length(.rebmix$pdf)) {
    .Object@Variables[which(pdf == .rebmix$pdf[i])] <- .rebmix$pdf.Variables[i]
  }

  # Dataset.

  for (i in 1:d) {
    if (.Object@Variables[i] == .rebmix$Variables[2]) {
      for (j in 1:length(Dataset)) {
         if (all(sapply(Dataset[[j]][, i], is.wholenumber)) == FALSE) {
           stop(sQuote("Dataset"), " all values in column ", i, " must be integers!", call. = FALSE)
         }

         if (any(sapply(Dataset[[j]][, i], function(x) x < 0)) == TRUE) {
           stop(sQuote("Dataset"), " all values in column ", i, " must be greater or equal than 0!", call. = FALSE)
         }
      }
    }
  }
 
### Panic Branislav.
  if (missing(EMcontrol) || length(EMcontrol) == 0) {
    EMcontrol <- new("EM.Control")
  }

  if (as.character(class(EMcontrol)) != "EM.Control") {
    stop(sQuote("EMcontrol"), " object of class EM.Control is requested!", call. = FALSE)
  }
### End  

  .Object@Dataset <- Dataset
  .Object@Preprocessing <- Preprocessing
  .Object@cmax <- cmax
  .Object@cmin <- cmin
  .Object@Criterion <- Criterion
  .Object@pdf <- pdf
  .Object@theta1 <- theta1
  .Object@theta2 <- theta2
  .Object@theta3 <- theta3
  .Object@K <- K
  .Object@ymin <- ymin
  .Object@ymax <- ymax
  .Object@ar <- ar
  .Object@Restraints <- Restraints
### Panic Branislav.  
  .Object@EMcontrol <- EMcontrol
### End

  rm(list = ls()[!(ls() %in% c(".Object"))])

  .Object
}) ## initialize

setMethod("show",
          signature(object = "REBMIX"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(object@pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  length(object@pos) <- 1

  if ((object@pos < 1) || (object@pos > nrow(object@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(object@summary), "!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"w\":", "\n", sep = "")

  print(object@w[[object@pos]], quote = FALSE)

  cat("Slot \"Theta\":", "\n", sep = "")

  print(object@Theta[[object@pos]], quote = FALSE)

  cat("Slot \"summary\":", "\n", sep = "")

  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object@summary), nomatch = 0)

  DF <- object@summary[object@pos, p]

  is.num <- sapply(DF, is.number); DF[is.num] <- lapply(DF[is.num], as.number)

  print(DF, quote = FALSE)

  rm(list = ls())
}) ## show

## Class REBMVNORM

setClass("REBMVNORM", contains = "REBMIX")

setMethod("show",
          signature(object = "REBMVNORM"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMVNORM is requested!", call. = FALSE)
  }

  if (!is.wholenumber(object@pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  length(object@pos) <- 1

  if ((object@pos < 1) || (object@pos > nrow(object@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(object@summary), "!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"w\":", "\n", sep = "")

  print(object@w[[object@pos]], quote = FALSE)

  cat("Slot \"Theta\":", "\n", sep = "")

  print(object@Theta[[object@pos]], quote = FALSE)

  cat("Slot \"summary\":", "\n", sep = "")

  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object@summary), nomatch = 0)

  DF <- object@summary[object@pos, p]

  is.num <- sapply(DF, is.number); DF[is.num] <- lapply(DF[is.num], as.number)

  print(DF, quote = FALSE)

  rm(list = ls())
}) ## show

## Class REBMIX.boot

setClass("REBMIX.boot",
slots = c(x = "ANY",
  rseed = "numeric",
  pos = "numeric",
  Bootstrap = "character",
  B = "numeric",
  n = "numeric",
  replace = "logical",
  prob = "numeric",
  c = "numeric",
  c.se = "numeric",
  c.cv = "numeric",
  c.mode = "numeric",
  c.prob = "numeric",
  w = "matrix",
  w.se = "numeric",
  w.cv = "numeric",
  Theta = "list",
  Theta.se = "list",
  Theta.cv = "list"),
prototype = list(rseed = -1,
  pos = 1,
  Bootstrap = "parametric",
  B = 100,
  replace = TRUE))

setMethod("initialize", "REBMIX.boot",
function(.Object, ...,
  x,
  rseed,
  pos,
  Bootstrap,
  B,
  n,
  replace,
  prob)
{
  model <- gsub("\\.boot", "", .Object@class[1])

  # x.

  if (missing(x) || (length(x) == 0)) {
    stop(sQuote("x"), " must not be empty!", call. = FALSE)
  }

  if (as.character(class(x)) != model) {
    stop(sQuote("x"), " object of class ", model, " is requested!", call. = FALSE)
  }

  # rseed.

  if (missing(rseed) || (length(rseed) == 0)) rseed <- .Object@rseed

  if (!is.wholenumber(rseed)) {
    stop(sQuote("rseed"), " integer is requested!", call. = FALSE)
  }

  length(rseed) <- 1

  if (rseed > -1) {
    stop(sQuote("rseed"), " must be less than 0!", call. = FALSE)
  }

  # pos.

  if (missing(pos) || (length(pos) == 0)) pos <- .Object@pos

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }

  # Bootstrap.

  if (missing(Bootstrap) || (length(Bootstrap) == 0)) Bootstrap <- .Object@Bootstrap

  if (!is.character(Bootstrap)) {
    stop(sQuote("Bootstrap"), " character is requested!", call. = FALSE)
  }

  Bootstrap <- match.arg(Bootstrap, .rebmix.boot$Bootstrap, several.ok = FALSE)

  # B.

  if (missing(B) || (length(B) == 0)) B <- .Object@B

  if (!is.wholenumber(B)) {
    stop(sQuote("B"), " integer is requested!", call. = FALSE)
  }

  length(B) <- 1

  if (B < 1) {
    stop(sQuote("B"), " must be greater than 0!", call. = FALSE)
  }

  # n.

  if (missing(n) || (length(n) == 0)) n <- .Object@n

  nmax <- nrow(as.matrix(x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]]))

  if (length(n) == 0) {
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

  # replace.

  if (missing(replace) || (length(replace) == 0)) replace <- .Object@replace

  if (!is.logical(replace)) {
    stop(sQuote("replace"), " logical is requested!", call. = FALSE)
  }

  # prob.

  if (missing(prob) || (length(prob) == 0)) {
    prob <- .Object@prob
  }
  else {
    if (!is.numeric(prob)) {
      stop(sQuote("prob"), " numeric vector is requested!", call. = FALSE)
    }

    if (length(prob) != length(n)) {
      stop("lengths of ", sQuote("prob"), " and ", sQuote("n"), " must match!", call. = FALSE)
    }
  }

  .Object@x <- x
  .Object@rseed <- rseed
  .Object@pos <- pos
  .Object@Bootstrap <- Bootstrap
  .Object@B <- B
  .Object@n <- n
  .Object@replace <- replace
  .Object@prob <- prob

  rm(list = ls()[!(ls() %in% c(".Object"))])

  .Object
}) ## initialize

setMethod("show",
          signature(object = "REBMIX.boot"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX.boot is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"c\":", "\n", sep = "")

  print(object@c, quote = FALSE)

  cat("Slot \"c.se\":", "\n", sep = "")

  print(object@c.se, quote = FALSE)

  cat("Slot \"c.cv\":", "\n", sep = "")

  print(object@c.cv, quote = FALSE)

  cat("Slot \"c.mode\":", "\n", sep = "")

  print(object@c.mode, quote = FALSE)

  cat("Slot \"c.prob\":", "\n", sep = "")

  print(object@c.prob, quote = FALSE)

  rm(list = ls())
}) ## show

## Class REBMVNORM.boot

setClass("REBMVNORM.boot", contains = "REBMIX.boot")

setMethod("show",
          signature(object = "REBMVNORM.boot"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMVNORM.boot is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"c\":", "\n", sep = "")

  print(object@c, quote = FALSE)

  cat("Slot \"c.se\":", "\n", sep = "")

  print(object@c.se, quote = FALSE)

  cat("Slot \"c.cv\":", "\n", sep = "")

  print(object@c.cv, quote = FALSE)

  cat("Slot \"c.mode\":", "\n", sep = "")

  print(object@c.mode, quote = FALSE)

  cat("Slot \"c.prob\":", "\n", sep = "")

  print(object@c.prob, quote = FALSE)

  rm(list = ls())
}) ## show

## Class RCLRMIX

setClass("RCLRMIX",
slots = c(x = "ANY",
  Dataset = "ANY",
  pos = "numeric",
  Zt = "factor",
  Zp = "factor",
  c = "numeric",
  p = "numeric",
  pi = "list",
  P = "data.frame",
  tau = "matrix",
  prob = "numeric",
### Panic Branislav.
  Rule = "character",
### End
  from = "numeric",
  to = "numeric",
  EN = "numeric",
  ED = "numeric"),
prototype = list(pos = 1))

setMethod("initialize", "RCLRMIX",
function(.Object, ...,
  x,
  Dataset,
  pos,
  Zt,
  Rule)  
{
  model <- gsub("RCLR", "REB", .Object@class[1])

  # x.

  if (missing(x) || (length(x) == 0)) {
    stop(sQuote("x"), " must not be empty!", call. = FALSE)
  }
  
  if (as.character(class(x)) != model) {
    stop(sQuote("x"), " object of class ", model, " is requested!", call. = FALSE)
  }

  # pos.

  if (missing(pos) || (length(pos) == 0)) pos <- .Object@pos

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  # Dataset.
  
  if (missing(Dataset) || (length(Dataset) == 0)) {
    Dataset <- x@Dataset[[pos]]
    
    if (as.character(class(Dataset)) == "data.frame") {
      n <- nrow(Dataset)
    }
    else
    if (as.character(class(Dataset)) == "Histogram") {
      n <- nrow(Dataset@Y)
    }     
  }
  else {
    if ((as.character(class(Dataset)) != "Histogram") && (as.character(class(Dataset)) != "data.frame")) {
      stop(sQuote("Dataset"), " data frame or object of class Histogram is requested!", call. = FALSE)
    }
    
    d <- length(x@Variables)

    if (as.character(class(Dataset)) == "data.frame") {
      if (ncol(Dataset) != d) {
        stop(sQuote("Dataset"), " number of columns in data frame must equal ", d, "!", call. = FALSE)
      }
      
      n <- nrow(Dataset)

      if (n < 1) {
        stop(sQuote("Dataset"), " number of rows in data frame must be greater than 0!", call. = FALSE)
      }      
    }
    else
    if (as.character(class(Dataset)) == "Histogram") {
      if (length(Dataset@h) != d) {
        stop(sQuote("Dataset"), " number of variables in object of class Histogram must equal ", d, "!", call. = FALSE)
      }
      
      n <- nrow(Dataset@Y)

      if (n < 1) {
        stop(sQuote("Dataset"), " number of bins in object of class Histogram must be greater than 0!", call. = FALSE)
      }       
    }      
  }

  # Zt.

  if (missing(Zt) || (length(Zt) == 0)) {
    Zt <- .Object@Zt
  }
  else {  
    if (!is.factor(Zt)) {
      stop(sQuote("Zt"), " factor is requested!", call. = FALSE)
    }
  
    if (n != length(Zt)) {
      stop("length of ", sQuote("Zt"), " must equal ", n, "!", call. = FALSE)  
    }
  }

  levels(Zt) <- 1:length(levels(Zt))
  
  # Rule

  if (missing(Rule) || (length(Rule) == 0)) {
    Rule <- .rclrmix$Rule[[1]]
  }
  else {
    if (!is.character(Rule)) {
      stop(sQuote("Rule"), " character is requested!", call. = FALSE)
    }
  }

  Rule <- match.arg(Rule, .rclrmix$Rule, several.ok = FALSE)  

  .Object@x <- x
  .Object@Dataset <- Dataset
  .Object@pos <- pos
  .Object@Zt <- Zt
  .Object@Rule <- Rule  

  rm(list = ls()[!(ls() %in% c(".Object"))])

  .Object
}) ## initialize

setMethod("show",
          signature(object = "RCLRMIX"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLRMIX is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"prob\":", "\n", sep = "")

  print(object@prob, quote = FALSE)

  cat("Slot \"Zp\":", "\n", sep = "")

  print(object@Zp, quote = FALSE)

  cat("Slot \"from\":", "\n", sep = "")

  print(object@from, quote = FALSE)

  cat("Slot \"to\":", "\n", sep = "")

  print(object@to, quote = FALSE)

  cat("Slot \"EN\":", "\n", sep = "")

  print(object@EN, quote = FALSE)

  cat("Slot \"ED\":", "\n", sep = "")

  print(object@ED, quote = FALSE)
  
  cat("Slot \"Rule\":", "\n", sep = "")

  print(object@Rule, quote = FALSE)  

  rm(list = ls())
}) ## show

## Class RCLRMVNORM

setClass("RCLRMVNORM", contains = "RCLRMIX")

setMethod("show",
          signature(object = "RCLRMVNORM"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLRMVNORM is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"prob\":", "\n", sep = "")

  print(object@prob, quote = FALSE)

  cat("Slot \"Zp\":", "\n", sep = "")

  print(object@Zp, quote = FALSE)

  cat("Slot \"from\":", "\n", sep = "")

  print(object@from, quote = FALSE)

  cat("Slot \"to\":", "\n", sep = "")

  print(object@to, quote = FALSE)

  cat("Slot \"EN\":", "\n", sep = "")

  print(object@EN, quote = FALSE)

  cat("Slot \"ED\":", "\n", sep = "")

  print(object@ED, quote = FALSE)
  
  cat("Slot \"Rule\":", "\n", sep = "")

  print(object@Rule, quote = FALSE)  

  rm(list = ls())
}) ## show

## Class RCLS.chunk

setClass("RCLS.chunk",
slots = c(s = "numeric",
  levels = "character",
  ntrain = "numeric",
  train = "list",
  Zr = "list",
  ntest = "numeric",
  test = "data.frame",
  Zt = "factor"))

setMethod("show",
          signature(object = "RCLS.chunk"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLS.chunk is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"s\":", "\n", sep = "")

  print(object@s, quote = FALSE)

  cat("Slot \"levels\":", "\n", sep = "")

  print(object@levels, quote = FALSE)

  cat("Slot \"ntrain\":", "\n", sep = "")

  print(object@ntrain, quote = FALSE)

  cat("Slot \"ntest\":", "\n", sep = "")

  print(object@ntest, quote = FALSE)

  rm(list = ls())
}) ## show

## Class RCLSMIX

setClass("RCLSMIX",
slots = c(x = "list",
  Dataset = "data.frame",
  o = "numeric",
  s = "numeric",
  ntrain = "numeric",
  P = "numeric",
  ntest = "numeric",
  Zt = "factor",
  Zp = "factor",
  CM = "table",
  Accuracy = "numeric",
  Error = "numeric",
  Precision = "numeric",
  Sensitivity = "numeric",
  Specificity = "numeric",
  Chunks = "numeric"),
prototype = list(CM = table(0)))

setMethod("initialize", "RCLSMIX",
function(.Object, ...,
  x,
  Dataset,
  Zt)
{
  model <- gsub("RCLS", "REB", .Object@class[1])

  # x.

  if (missing(x) || (length(x) == 0)) {
    stop(sQuote("x"), " must not be empty!", call. = FALSE)
  }

  if (!all(unlist(lapply(x, class)) == model)) {
    stop(sQuote("x"), " list of ", model, " objects is requested!", call. = FALSE)
  }

  # o.

  .Object@o <- length(x)

  # s.

  s <- unique(unlist(lapply(x, function(x) length(x@Dataset))))

  if (length(s) != 1) {
    stop("lengths of ", sQuote("Dataset"), " in ", sQuote("x"), " must be equal!", call. = FALSE)
  }

  if (s == 1) {
    stop(sQuote("s"), " must be greater than 1!", call. = FALSE)
  }

  .Object@s <- s

  # ntrain.

  ntrain <- matrix(unlist(lapply(x, function(x) lapply(x@Dataset, nrow))), ncol = s, byrow = TRUE)

  ntrain <- ntrain[!duplicated(ntrain),]

  if (length(ntrain) != s) {
    stop(sQuote("Dataset"), " in ", sQuote("x"), " numbers of rows in data frames must be equal!", call. = FALSE)
  }

  .Object@ntrain <- ntrain

  # P.

  .Object@P <- ntrain / sum(ntrain)

  # Dataset.

  if (missing(Dataset) || (length(Dataset) == 0)) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }

  if (!is.data.frame(Dataset)) {
    stop(sQuote("Dataset"), " data frame is requested!", call. = FALSE)
  }
  
  n <- nrow(Dataset)

  # ntest.

  .Object@ntest <- nrow(Dataset)

  # Zt.

  if (missing(Zt) || (length(Zt) == 0)) Zt <- .Object@Zt

  if (!is.factor(Zt)) {
    stop(sQuote("Zt"), " factor is requested!", call. = FALSE)
  }
  
  if (n != length(Zt)) {
    stop("length of ", sQuote("Zt"), " must equal ", n, "!", call. = FALSE)  
  }

  levels(Zt) <- 1:length(levels(Zt))

  .Object@x <- x
  .Object@Dataset <- Dataset
  .Object@Zt <- Zt

  rm(list = ls()[!(ls() %in% c(".Object"))])

  .Object
}) ## initialize

setMethod("show",
          signature(object = "RCLSMIX"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLSMIX is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"CM\":", "\n", sep = "")

  print(object@CM, quote = FALSE)

  cat("Slot \"Error\":", "\n", sep = "")

  print(object@Error, quote = FALSE)

  cat("Slot \"Precision\":", "\n", sep = "")

  names(object@Precision) <- NULL

  print(object@Precision, quote = FALSE)

  cat("Slot \"Sensitivity\":", "\n", sep = "")

  names(object@Sensitivity) <- NULL

  print(object@Sensitivity, quote = FALSE)

  cat("Slot \"Specificity\":", "\n", sep = "")

  names(object@Specificity) <- NULL

  print(object@Specificity, quote = FALSE)

  cat("Slot \"Chunks\":", "\n", sep = "")

  names(object@Chunks) <- NULL

  print(object@Chunks, quote = FALSE)

  rm(list = ls())
}) ## show

## Class RCLSMVNORM

setClass("RCLSMVNORM", contains = "RCLSMIX")

setMethod("show",
          signature(object = "RCLSMVNORM"),
function(object)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLSMVNORM is requested!", call. = FALSE)
  }

  cat("An object of class ", "\"", as.character(class(object)), "\"", "\n", sep = "")

  cat("Slot \"CM\":", "\n", sep = "")

  print(object@CM, quote = FALSE)

  cat("Slot \"Error\":", "\n", sep = "")

  print(object@Error, quote = FALSE)

  cat("Slot \"Precision\":", "\n", sep = "")

  names(object@Precision) <- NULL

  print(object@Precision, quote = FALSE)

  cat("Slot \"Sensitivity\":", "\n", sep = "")

  names(object@Sensitivity) <- NULL

  print(object@Sensitivity, quote = FALSE)

  cat("Slot \"Specificity\":", "\n", sep = "")

  names(object@Specificity) <- NULL

  print(object@Specificity, quote = FALSE)

  cat("Slot \"Chunks\":", "\n", sep = "")

  names(object@Chunks) <- NULL

  print(object@Chunks, quote = FALSE)

  rm(list = ls())
}) ## show

# Class Histogram

setClass(Class = "Histogram",
slots = c(Y = "data.frame",
  K = "numeric",
  ymin = "numeric",
  ymax = "numeric",
  y0 = "numeric",  
  h = "numeric",
  n = "numeric",
  ns = "numeric"))

setMethod("initialize", "Histogram",
function(.Object, ...,
  Y,
  K,
  ymin,
  ymax,
  y0,
  h,
  n,
  ns)
{
  # Y.

  if (missing(Y) || (length(Y) == 0)) {
    stop(sQuote("Y"), " must not be empty!", call. = FALSE)
  }

  if (!is.data.frame(Y)) {
    stop(sQuote("Y"), " data frame is requested!", call. = FALSE)
  }

  d <- ncol(Y) - 1

  if (d < 1) {
    stop(sQuote("Y"), " number of columns in data frame must be greater than 1!", call. = FALSE)
  }

  Y <- as.data.frame(Y[complete.cases(Y), ])

  if (nrow(Y) < 1) {
    stop(sQuote("Y"), " number of rows in data frame must be greater than 0!", call. = FALSE)
  }
  
  colnames(Y) <- c(paste("y", if (d > 1) 1:d else "", sep = ""), "k")
  
  # K.

  if (missing(K) || (length(K) == 0)) {
    stop(sQuote("K"), " must not be empty!", call. = FALSE)
  }

  if (!is.wholenumber(K)) {
    stop(sQuote("K"), " integer or integer vector is requested!", call. = FALSE)
  }

  if (!all(K > 0)) {
    stop("all ", sQuote("K"), " must be greater than 0!", call. = FALSE)
  }
  
  if (length(K) == 1) {
    K <- rep(K, d)
  }
  else
  if (length(K) != d) {
    stop("lengths of ", sQuote("K"), " and ", sQuote("d"), " must match!", call. = FALSE)
  }  
  
  # ymin.

  if (missing(ymin) || (length(ymin) == 0)) {
    stop(sQuote("ymin"), " must not be empty!", call. = FALSE)
  }
  else {
    if (length(ymin) != d) {
      stop("lengths of ", sQuote("ymin"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
  }
  
  # ymax.

  if (missing(ymax) || (length(ymax) == 0)) {
    stop(sQuote("ymax"), " must not be empty!", call. = FALSE)
  }
  else {
    if (length(ymax) != d) {
      stop("lengths of ", sQuote("ymax"), " and ", sQuote("d"), " must match!", call. = FALSE)
    }
  }
  
  h <- (ymax - ymin) / K; y0 <- ymin + 0.5 * h
  
  n <- 0
  
  ns <- 0
  
  .Object@Y <- Y
  .Object@K <- K
  .Object@ymin <- ymin
  .Object@ymax <- ymax
  .Object@y0 <- y0
  .Object@h <- h
  .Object@n <- n
  .Object@ns <- ns

  rm(list = ls()[!(ls() %in% c(".Object"))])

  .Object
}) ## initialize
