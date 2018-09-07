setMethod("a.c", signature(x = "RNGMIX.Theta"), function(x) x@c)
setMethod("a.d", signature(x = "RNGMIX.Theta"), function(x) x@d)
setMethod("a.pdf", signature(x = "RNGMIX.Theta"), function(x) x@pdf)
setMethod("a.Theta", signature(x = "RNGMIX.Theta"), function(x) x@Theta)

setMethod("a.theta1<-", 
  signature(x = "RNGMIX.Theta", l = "missing"),
function(x, value)
{
  x@d <- 1; length(x@pdf) <- 1

  # value.

  if (missing(value) || (length(value) == 0)) {
    stop(sQuote("value"), " must not be empty!", call. = FALSE)
  }

  if (!is.number(value)) {
    stop(sQuote("value"), " numeric vector is requested!", call. = FALSE)
  } 

  if (length(value) != x@c) {
    stop("length of ", sQuote("value"), " must equal " , x@c, "!", call. = FALSE)
  }
  
  for (l in 1:x@c) {
    if (x@pdf == .rebmix$pdf[3]) {
      if (value[l] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[3]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
    else
    if (x@pdf == .rebmix$pdf[4]) {
      if (!is.wholenumber(value[l])) {
        stop(sQuote("value"), " integer is requested for ", dQuote(.rebmix$pdf[4]), "!", call. = FALSE)
      }

      if (value[l] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[4]), " must be greater or equal than 0!", call. = FALSE)
      }
    }
    else
    if (x@pdf == .rebmix$pdf[5]) {
      if (value[l] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[5]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
    else
    if (x@pdf == .rebmix$pdf[7]) {
      if (value[l] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[7]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
  }
  
  for (l in 1:x@c) {
    length(x@Theta[[1 + (l - 1) * 3]]) <- 1
    x@Theta[[2 + (l - 1) * 3]] <- value[l]
    length(x@Theta[[3 + (l - 1) * 3]]) <- 1
  }
  
  rm(list = ls()[!(ls() %in% c("x"))])
  
  x
}) ## a.theta1<-

setMethod("a.theta1<-", 
  signature(x = "RNGMIX.Theta"),
function(x, l, value)
{
  # l.

  if (missing(l) || (length(l) == 0)) {
    stop(sQuote("l"), " must not be empty!", call. = FALSE)
  }
  
  if (!is.wholenumber(l)) {
    stop(sQuote("l"), " integer is requested!", call. = FALSE)
  }

  length(l) <- 1

  if ((l < 1) || (l > x@c)) {
    stop(sQuote("l"), " must be greater than 0 and less or equal than ", x@c, "!", call. = FALSE)
  }
  
  # value.

  if (missing(value) || (length(value) == 0)) {
    stop(sQuote("value"), " must not be empty!", call. = FALSE)
  }

  if (!is.number(value)) {
    stop(sQuote("value"), " numeric vector is requested!", call. = FALSE)
  } 

  if (length(value) != x@d) {
    stop("length of ", sQuote("value"), " must equal " , x@d, "!", call. = FALSE)
  }
  
  for (i in 1:x@d) {
    if (x@pdf[i] == .rebmix$pdf[3]) {
      if (value[i] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[3]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
    else
    if (x@pdf[i] == .rebmix$pdf[4]) {
      if (!is.wholenumber(value[i])) {
        stop(sQuote("value"), " integer is requested for ", dQuote(.rebmix$pdf[4]), "!", call. = FALSE)
      }

      if (value[i] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[4]), " must be greater or equal than 0!", call. = FALSE)
      }
    }
    else
    if (x@pdf[i] == .rebmix$pdf[5]) {
      if (value[i] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[5]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
    else
    if (x@pdf[i] == .rebmix$pdf[7]) {
      if (value[i] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[7]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
  }

  x@Theta[[2 + (l - 1) * 3]] <- value
  
  rm(list = ls()[!(ls() %in% c("x"))])
  
  x
}) ## a.theta1<-

setMethod("a.theta2<-", 
  signature(x = "RNGMIX.Theta", l = "missing"),
function(x, value)
{
  x@d <- 1; length(x@pdf) <- 1

  # value.

  if (missing(value) || (length(value) == 0)) {
    stop(sQuote("value"), " must not be empty!", call. = FALSE)
  }

  if (!is.numeric(value)) {
    stop(sQuote("value"), " numeric vector is requested!", call. = FALSE)
  } 

  if (length(value) != x@c) {
    stop("length of ", sQuote("value"), " must equal " , x@c, "!", call. = FALSE)
  }
  
  for (l in 1:x@c) {
    if (x@pdf == .rebmix$pdf[1]) {
      if (value[l] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[1]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
    else
    if (x@pdf == .rebmix$pdf[2]) {
      if (value[l] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[2]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
    else
    if (x@pdf == .rebmix$pdf[3]) {
      if (value[l] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[3]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
    else
    if (x@pdf == .rebmix$pdf[4]) {
      if ((value[l] < 0.0) || (value[l] > 1.0)) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[4]), " must be greater or equal than 0.0 and less or equal than 1.0!", call. = FALSE)
      }    
    }
    else
    if (x@pdf == .rebmix$pdf[5]) {
      value[l] <- NA
    }
    else
    if (x@pdf == .rebmix$pdf[6]) {
      value[l] <- NA
    }
    else    
    if (x@pdf == .rebmix$pdf[7]) {
      if (value[l] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[7]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
    else    
    if (x@pdf == .rebmix$pdf[9]) {
      if (value[l] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[9]), " must be greater than 0.0!", call. = FALSE)
      }    
    }    
  }
  
  for (l in 1:x@c) {
    length(x@Theta[[1 + (l - 1) * 3]]) <- 1
    length(x@Theta[[2 + (l - 1) * 3]]) <- 1
    x@Theta[[3 + (l - 1) * 3]] <- value[l]
  }
  
  rm(list = ls()[!(ls() %in% c("x"))])  
  
  x
}) ## a.theta2<-

setMethod("a.theta2<-", 
  signature(x = "RNGMIX.Theta"),
function(x, l, value)
{
  # l.

  if (missing(l) || (length(l) == 0)) {
    stop(sQuote("l"), " must not be empty!", call. = FALSE)
  }
  
  if (!is.wholenumber(l)) {
    stop(sQuote("l"), " integer is requested!", call. = FALSE)
  }

  length(l) <- 1

  if ((l < 1) || (l > x@c)) {
    stop(sQuote("l"), " must be greater than 0 and less or equal than ", x@c, "!", call. = FALSE)
  }
  
  # value.

  if (missing(value) || (length(value) == 0)) {
    stop(sQuote("value"), " must not be empty!", call. = FALSE)
  }

  if (!is.numeric(value)) {
    stop(sQuote("value"), " numeric vector is requested!", call. = FALSE)
  } 

  if (length(value) != x@d) {
    stop("length of ", sQuote("value"), " must equal " , x@d, "!", call. = FALSE)
  }
  
  for (i in 1:x@d) {
    if (x@pdf[i] == .rebmix$pdf[1]) {
      if (value[i] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[1]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
    else
    if (x@pdf[i] == .rebmix$pdf[2]) {
      if (value[i] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[2]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
    else
    if (x@pdf[i] == .rebmix$pdf[3]) {
      if (value[i] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[3]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
    else
    if (x@pdf[i] == .rebmix$pdf[4]) {
      if ((value[i] < 0.0) || (value[i] > 1.0)) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[4]), " must be greater or equal than 0.0 and less or equal than 1.0!", call. = FALSE)
      }    
    }
    else
    if (x@pdf[i] == .rebmix$pdf[5]) {
      value[i] <- NA
    }
    else
    if (x@pdf[i] == .rebmix$pdf[6]) {
      value[i] <- NA
    }
    else    
    if (x@pdf[i] == .rebmix$pdf[7]) {
      if (value[i] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[7]), " must be greater than 0.0!", call. = FALSE)
      }    
    }
    else    
    if (x@pdf[i] == .rebmix$pdf[9]) {
      if (value[i] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[9]), " must be greater than 0.0!", call. = FALSE)
      }    
    }    
  }  

  x@Theta[[3 + (l - 1) * 3]] <- value
  
  rm(list = ls()[!(ls() %in% c("x"))])  
  
  x
}) ## a.theta2<-

setMethod("a.theta1<-", 
  signature(x = "RNGMVNORM.Theta"),
function(x, l, value)
{
  # l.

  if (missing(l) || (length(l) == 0)) {
    stop(sQuote("l"), " must not be empty!", call. = FALSE)
  }
  
  if (!is.wholenumber(l)) {
    stop(sQuote("l"), " integer is requested!", call. = FALSE)
  }

  length(l) <- 1

  if ((l < 1) || (l > x@c)) {
    stop(sQuote("l"), " must be greater than 0 and less or equal than ", x@c, "!", call. = FALSE)
  }
  
  # value.

  if (missing(value) || (length(value) == 0)) {
    stop(sQuote("value"), " must not be empty!", call. = FALSE)
  }

  if (!is.number(value)) {
    stop(sQuote("value"), " numeric vector is requested!", call. = FALSE)
  } 

  if (length(value) != x@d) {
    stop("length of ", sQuote("value"), " must equal " , x@d, "!", call. = FALSE)
  }    

  x@Theta[[2 + (l - 1) * 3]] <- value
  
  rm(list = ls()[!(ls() %in% c("x"))])  
  
  x
}) ## a.theta1<-

setMethod("a.theta2<-", 
  signature(x = "RNGMVNORM.Theta"),
function(x, l, value)
{
  # l.

  if (missing(l) || (length(l) == 0)) {
    stop(sQuote("l"), " must not be empty!", call. = FALSE)
  }
  
  if (!is.wholenumber(l)) {
    stop(sQuote("l"), " integer is requested!", call. = FALSE)
  }

  length(l) <- 1

  if ((l < 1) || (l > x@c)) {
    stop(sQuote("l"), " must be greater than 0 and less or equal than ", x@c, "!", call. = FALSE)
  }
  
  # value.

  if (missing(value) || (length(value) == 0)) {
    stop(sQuote("value"), " must not be empty!", call. = FALSE)
  }

  if (!is.number(value)) {
    stop(sQuote("value"), " numeric vector is requested!", call. = FALSE)
  } 

  if (length(value) != x@d * x@d) {
    stop("length of ", sQuote("value"), " must equal " , x@d * x@d, "!", call. = FALSE)
  }

  x@Theta[[3 + (l - 1) * 3]] <- value
  
  rm(list = ls()[!(ls() %in% c("x"))])  
  
  x
}) ## a.theta2<-

setMethod("a.Dataset.name", signature(x = "RNGMIX"), function(x) x@Dataset.name)
setMethod("a.rseed", signature(x = "RNGMIX"), function(x) x@rseed)
setMethod("a.n", signature(x = "RNGMIX"), function(x) x@n)
setMethod("a.Theta", signature(x = "RNGMIX"), function(x) x@Theta)

setMethod("a.Dataset", 
          signature(x = "RNGMIX"), 
function(x, pos)
{
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > length(x@Dataset))) {
    output <- x@Dataset
  }
  else {
    output <- x@Dataset[[pos]]
  }
  
  rm(list = ls()[!(ls() %in% c("output"))])  
  
  output
}) ## Dataset
          
setMethod("a.Zt", signature(x = "RNGMIX"), function(x) x@Zt)
setMethod("a.w", signature(x = "RNGMIX"), function(x) x@w)
setMethod("a.Variables", signature(x = "RNGMIX"), function(x) x@Variables)
setMethod("a.ymin", signature(x = "RNGMIX"), function(x) x@ymin)
setMethod("a.ymax", signature(x = "RNGMIX"), function(x) x@ymax)

setMethod("a.Dataset", 
          signature(x = "REBMIX"), 
function(x, pos)
{
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > length(x@Dataset))) {
    output <- x@Dataset
  }
  else {
    output <- x@Dataset[[pos]]
  }
  
  rm(list = ls()[!(ls() %in% c("output"))])  
  
  output
}) ## Dataset

setMethod("a.Preprocessing", signature(x = "REBMIX"), function(x) x@Preprocessing)
setMethod("a.cmax", signature(x = "REBMIX"), function(x) x@cmax)
setMethod("a.cmin", signature(x = "REBMIX"), function(x) x@cmin)
setMethod("a.Criterion", signature(x = "REBMIX"), function(x) x@Criterion)
setMethod("a.Variables", signature(x = "REBMIX"), function(x) x@Variables)
setMethod("a.pdf", signature(x = "REBMIX"), function(x) x@pdf)
setMethod("a.theta1", signature(x = "REBMIX"), function(x) x@theta1)
setMethod("a.theta2", signature(x = "REBMIX"), function(x) x@theta2)

setMethod("a.theta1.all", 
  signature(x = "REBMIX"), 
function(x, pos) 
{
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1
  
  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  c <- length(x@w[[pos]]); d <- length(x@pdf)  

  output <- matrix(data = 0.0, nrow = c, ncol = d, byrow = TRUE)
  
  for (l in 1:c) {
    output[l, ] <- x@Theta[[pos]][[2 + (l - 1) * 3]]
  }
  
  colnames(output) <- NULL
  rownames(output) <- paste("theta2.", 1:c, sep = "")
  
  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## a.theta1.all

setMethod("a.theta2.all", 
  signature(x = "REBMIX"), 
function(x, pos) 
{
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1
  
  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  c <- length(x@w[[pos]]); d <- length(x@pdf)  

  output <- matrix(data = 0.0, nrow = c, ncol = d, byrow = TRUE)
  
  for (l in 1:c) {
    output[l, ] <- x@Theta[[pos]][[3 + (l - 1) * 3]]
  }
  
  colnames(output) <- NULL
  rownames(output) <- paste("theta2.", 1:c, sep = "")
  
  rm(list = ls()[!(ls() %in% c("output"))])  

  output
}) ## a.theta2.all

setMethod("a.theta2.all", 
  signature(x = "REBMVNORM"), 
function(x, pos) 
{
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1
  
  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  c <- length(x@w[[pos]]); d <- length(x@pdf)  

  output <- matrix(data = 0.0, nrow = c, ncol = d * d, byrow = TRUE)
  
  for (l in 1:c) {
    output[l, ] <- x@Theta[[pos]][[3 + (l - 1) * 3]]
  }
  
  colnames(output) <- NULL
  rownames(output) <- paste("theta2.", 1:c, sep = "")
  
  rm(list = ls()[!(ls() %in% c("output"))])  

  output
}) ## a.theta2.all

setMethod("a.K", signature(x = "REBMIX"), function(x) x@K)
setMethod("a.y0", signature(x = "REBMIX"), function(x) x@y0)
setMethod("a.ymin", signature(x = "REBMIX"), function(x) x@ymin)
setMethod("a.ymax", signature(x = "REBMIX"), function(x) x@ymax)
setMethod("a.ar", signature(x = "REBMIX"), function(x) x@ar)
setMethod("a.Restraints", signature(x = "REBMIX"), function(x) x@Restraints)

setMethod("a.w", 
          signature(x = "REBMIX"), 
function(x, pos)
{
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1
  
  if ((pos < 1) || (pos > length(x@w))) {
    output <- x@w
  }
  else {
    output <- x@w[[pos]]
  }
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
}) ## w

setMethod("a.Theta", 
          signature(x = "REBMIX"), 
function(x, pos)
{
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1
  
  if ((pos < 1) || (pos > length(x@Theta))) {
    output <- x@Theta
  }
  else {
    output <- x@Theta[[pos]]
  }
  
  rm(list = ls()[!(ls() %in% c("output"))])
  
  output
}) ## Theta

setMethod("a.summary", signature(x = "REBMIX"), function(x) x@summary)
setMethod("a.pos", signature(x = "REBMIX"), function(x) x@pos)
setMethod("a.opt.c", signature(x = "REBMIX"), function(x) x@opt.c)
setMethod("a.opt.IC", signature(x = "REBMIX"), function(x) x@opt.IC)
setMethod("a.opt.logL", signature(x = "REBMIX"), function(x) x@opt.logL)
setMethod("a.opt.D", signature(x = "REBMIX"), function(x) x@opt.D)
setMethod("a.all.K", signature(x = "REBMIX"), function(x) x@all.K)
setMethod("a.all.IC", signature(x = "REBMIX"), function(x) x@all.IC)

setMethod("a.rseed", signature(x = "REBMIX.boot"), function(x) x@rseed)
setMethod("a.pos", signature(x = "REBMIX.boot"), function(x) x@pos)
setMethod("a.Bootstrap", signature(x = "REBMIX.boot"), function(x) x@Bootstrap)
setMethod("a.B", signature(x = "REBMIX.boot"), function(x) x@B)
setMethod("a.n", signature(x = "REBMIX.boot"), function(x) x@n)
setMethod("a.replace", signature(x = "REBMIX.boot"), function(x) x@replace)
setMethod("a.prob", signature(x = "REBMIX.boot"), function(x) x@prob)
setMethod("a.c", signature(x = "REBMIX.boot"), function(x) x@c)
setMethod("a.c.se", signature(x = "REBMIX.boot"), function(x) x@c.se)
setMethod("a.c.cv", signature(x = "REBMIX.boot"), function(x) x@c.cv)
setMethod("a.c.mode", signature(x = "REBMIX.boot"), function(x) x@c.mode)
setMethod("a.c.prob", signature(x = "REBMIX.boot"), function(x) x@c.prob)
setMethod("a.w", signature(x = "REBMIX.boot"), function(x) x@w)
setMethod("a.w.se", signature(x = "REBMIX.boot"), function(x) x@w.se)
setMethod("a.w.cv", signature(x = "REBMIX.boot"), function(x) x@w.cv)
setMethod("a.Theta", signature(x = "REBMIX.boot"), function(x) x@Theta)
setMethod("a.Theta.se", signature(x = "REBMIX.boot"), function(x) x@Theta.se)
setMethod("a.Theta.cv", signature(x = "REBMIX.boot"), function(x) x@Theta.cv)

setMethod("a.pos", signature(x = "RCLRMIX"), function(x) x@pos)
setMethod("a.Zt", signature(x = "RCLRMIX"), function(x) x@Zt)
setMethod("a.Zp", signature(x = "RCLRMIX"), function(x) x@Zp)
setMethod("a.c", signature(x = "RCLRMIX"), function(x) x@c)
setMethod("a.p", signature(x = "RCLRMIX"), function(x) x@p)
setMethod("a.pi", signature(x = "RCLRMIX"), function(x) x@pi)
setMethod("a.P", signature(x = "RCLRMIX"), function(x) x@P)
setMethod("a.tau", signature(x = "RCLRMIX"), function(x) x@tau)
setMethod("a.prob", signature(x = "RCLRMIX"), function(x) x@prob)
setMethod("a.from", signature(x = "RCLRMIX"), function(x) x@from)
setMethod("a.to", signature(x = "RCLRMIX"), function(x) x@to)
setMethod("a.EN", signature(x = "RCLRMIX"), function(x) x@EN)
setMethod("a.ED", signature(x = "RCLRMIX"), function(x) x@ED)

setMethod("a.s", signature(x = "RCLS.chunk"), function(x) x@s)
setMethod("a.levels", signature(x = "RCLS.chunk"), function(x) x@levels)
setMethod("a.ntrain", signature(x = "RCLS.chunk"), function(x) x@ntrain)
setMethod("a.train", signature(x = "RCLS.chunk"), function(x) x@train)
setMethod("a.ntest", signature(x = "RCLS.chunk"), function(x) x@ntest)
setMethod("a.test", signature(x = "RCLS.chunk"), function(x) x@test)
setMethod("a.Zt", signature(x = "RCLS.chunk"), function(x) x@Zt)

setMethod("a.o", signature(x = "RCLSMIX"), function(x) x@o)
setMethod("a.Dataset", signature(x = "RCLSMIX"), function(x) x@Dataset)
setMethod("a.s", signature(x = "RCLSMIX"), function(x) x@s)
setMethod("a.ntrain", signature(x = "RCLSMIX"), function(x) x@ntrain)
setMethod("a.P", signature(x = "RCLSMIX"), function(x) x@P)
setMethod("a.ntest", signature(x = "RCLSMIX"), function(x) x@ntest)
setMethod("a.Zt", signature(x = "RCLSMIX"), function(x) x@Zt)
setMethod("a.Zp", signature(x = "RCLSMIX"), function(x) x@Zp)
setMethod("a.CM", signature(x = "RCLSMIX"), function(x) x@CM)
setMethod("a.Accuracy", signature(x = "RCLSMIX"), function(x) x@Accuracy)
setMethod("a.Error", signature(x = "RCLSMIX"), function(x) x@Error)
setMethod("a.Precision", signature(x = "RCLSMIX"), function(x) x@Precision)
setMethod("a.Sensitivity", signature(x = "RCLSMIX"), function(x) x@Sensitivity)
setMethod("a.Specificity", signature(x = "RCLSMIX"), function(x) x@Specificity)
setMethod("a.Chunks", signature(x = "RCLSMIX"), function(x) x@Chunks)
