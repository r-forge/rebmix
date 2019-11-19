### Panic Branislav.

setMethod("a.strategy", signature(x = "EM.Control"), function(x) x@strategy)
setMethod("a.variant", signature(x = "EM.Control"), function(x) x@variant)
setMethod("a.acceleration", signature(x = "EM.Control"), function(x) x@acceleration)
setMethod("a.tolerance", signature(x = "EM.Control"), function(x) x@tolerance)
setMethod("a.acceleration.multiplier", signature(x = "EM.Control"), function(x) x@acceleration.multiplier)
setMethod("a.maximum.iterations", signature(x = "EM.Control"), function(x) x@maximum.iterations)

setMethod("a.strategy<-",
          signature = (x = "EM.Control"),
function(x, value)
{
  # value.

  if (missing(value) || (length(value) == 0)) {
    stop(sQuote("strategy"), " must not be empty!", call. = FALSE)
  }

  if (!is.character(value)) {
    stop(sQuote("strategy"), " must be character!", call. = FALSE)
  }

  x@strategy <- match.arg(value, .rebmix$EMStrategy)

  rm(list = ls()[!(ls() %in% c("x"))])

  x
}) ## a.strategy<-

setMethod("a.variant<-",
          signature = (x = "EM.Control"),
function(x, value)
{
  # value.

  if (missing(value) || (length(value) == 0)) {
    stop(sQuote("variant"), " must not be empty!", call. = FALSE)
  }

  if (!is.character(value)) {
    stop(sQuote("variant"), " must be character!", call. = FALSE)
  }

  x@variant <- match.arg(value, .rebmix$EMVariant)

  rm(list = ls()[!(ls() %in% c("x"))])

  x
}) ## a.variant<-

setMethod("a.acceleration<-",
          signature = (x = "EM.Control"),
function(x, value)
{
  # value.

  if (missing(value) || (length(value) == 0)) {
    stop(sQuote("acceleration"), " must not be empty!", call. = FALSE)
  }

  if (!is.character(value)){
    stop(sQuote("acceleration"), " must be character!", call. = FALSE)
  }

  x@acceleration <- match.arg(value, .rebmix$EMAcceleration)

  rm(list = ls()[!(ls() %in% c("x"))])

  x
}) ## a.acceleration<-

setMethod("a.tolerance<-",
          signature = (x = "EM.Control"),
function(x, value)
{
  # value.

  if (missing(value) || (length(value) == 0)) {
    stop(sQuote("tolerance"), " must not be empty!", call. = FALSE)
  }

  if (!is.numeric(value)) {
    stop(sQuote("tolerance"), " numeric is requested!", call. = FALSE)
  }

  length(value) <- 1

  if (value <= 0.0) {
    stop(sQuote("tolerance"), " must be greater than 0.0!", call. = FALSE)
  }   

  x@tolerance <- value

  rm(list = ls()[!(ls() %in% c("x"))])

  x
}) ## a.tolerance<-

setMethod("a.acceleration.multiplier<-",
          signature = (x = "EM.Control"),
function(x, value)
{
  # value.
  
  if (missing(value) || (length(value) == 0)) {
    stop(sQuote("acceleration.multiplier"), " must not be empty!", call. = FALSE)
  }

  if (!is.numeric(value)) {
    stop(sQuote("acceleration.multiplier"), " numeric is requested!", call. = FALSE)
  }

  length(value) <- 1

  if (value < 1.0 || value > 2.0) {
    stop(sQuote("acceleration.multiplier"), " must be greater or equal than 1.0 and less or equal than 2.0!", call. = FALSE)
  }  

  x@acceleration.multiplier <- value

  rm(list = ls()[!(ls() %in% c("x"))])

  x
}) ## a.acceleration.multiplier<-

setMethod("a.maximum.iterations<-",
          signature = (x = "EM.Control"),
function(x, value)
{
  # value.
  
  if (missing(value) || (length(value) == 0)) {
    stop(sQuote("maximum.iterations"), " must not be empty!", call. = FALSE)
  }

  if (!is.wholenumber(value)) {
    stop(sQuote("maximum.iterations"), " integer is requested!", call. = FALSE)
  }

  length(value) <- 1

  if (value < 1) {
    stop(sQuote("maximum.iterations"), " must be greater than 0!", call. = FALSE)
  }  

  x@maximum.iterations <- value

  rm(list = ls()[!(ls() %in% c("x"))])

  x
}) ## a.maximum.iterations<-

### End

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
    else
    if (x@pdf == .rebmix$pdf[10]) {
      if (value[l] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[10]), " must be greater than 0.0!", call. = FALSE)
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
    else
    if (x@pdf[i] == .rebmix$pdf[10]) {
      if (value[i] < 0.0) {
        stop(sQuote("value"), " for ", dQuote(.rebmix$pdf[10]), " must be greater than 0.0!", call. = FALSE)
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
}) ## a.Dataset

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
}) ## a.Dataset

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
}) ## a.w

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
}) ## a.Theta

setMethod("a.summary",
          signature(x = "REBMIX"),
function(x, pos, col.name)
{
  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }

  length(pos) <- 1
  
  if ((pos < 1) || (pos > nrow(x@summary))) {
    output <- x@summary
  }
  else {
    output <- x@summary[pos, ]
  }

  if (!missing(col.name) && (length(col.name) > 0)) {
    if (!is.character(col.name)) {
      stop(sQuote("col.name"), " character is requested!", call. = FALSE)
    }

    col.name <- match.arg(col.name, colnames(output), several.ok = FALSE)

    output <- output[, col.name]

    if (is.number(output) == TRUE) {
      output <- as.numeric(output)
    }
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  output
}) ## a.summary

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

setMethod("a.Zp",
          signature(x = "RCLRMIX"),
function(x, s)
{
  Zp <- as.numeric(levels(x@Zp))[x@Zp]

  c <- x@c; s <- eval(s)

  if (!is.wholenumber(s)) {
    stop(sQuote("s"), " integer is requested!", call. = FALSE)
  }

  length(s) <- 1

  if ((s < 1) || (s > c)) {
    stop(sQuote("s"), " must be greater than 0 and less or equal than ", c, "!", call. = FALSE)
  }

  l <- c - 1

  while (s < length(unique(Zp))) {
    Zp[Zp == x@from[l]] <- x@to[l]

    l <- l - 1
  }

  rm(list = ls()[!(ls() %in% c("Zp"))])

  as.factor(Zp)
}) ## a.Zp

setMethod("a.c", signature(x = "RCLRMIX"), function(x) x@c)

setMethod("a.p",
          signature(x = "RCLRMIX"),
function(x, s)
{
  p <- x@p

  c <- x@c; s <- eval(s)

  if (!is.wholenumber(s)) {
    stop(sQuote("s"), " integer is requested!", call. = FALSE)
  }

  length(s) <- 1

  if ((s < 1) || (s > c)) {
    stop(sQuote("s"), " must be greater than 0 and less or equal than ", c, "!", call. = FALSE)
  }

  l <- c - 1; C <- numeric()

  while (l + 1 > s) {
    p[x@to[l]] <- p[x@to[l]] + p[x@from[l]]

    C <- c(C, x@from[l]); l <- l - 1
  }

  if (length(C) > 0) {
    p <- p[-C]
  }

  rm(list = ls()[!(ls() %in% c("p"))])

  p
}) ## a.p

setMethod("a.pi",
          signature(x = "RCLRMIX"),
function(x, s)
{
  p <- x@p; pi <- x@pi

  c <- x@c; d <- length(x@pi); s <- eval(s)

  if (!is.wholenumber(s)) {
    stop(sQuote("s"), " integer is requested!", call. = FALSE)
  }

  length(s) <- 1

  if ((s < 1) || (s > c)) {
    stop(sQuote("s"), " must be greater than 0 and less or equal than ", c, "!", call. = FALSE)
  }

  l <- c - 1; C <- numeric()

  while (l + 1 > s) {
    for (i in 1:d) {
      pi[[i]][x@to[l], ] <- (p[x@to[l]] * pi[[i]][x@to[l], ] + p[x@from[l]] * pi[[i]][x@from[l], ]) / (p[x@to[l]] + p[x@from[l]])
    }

    p[x@to[l]] <- p[x@to[l]] + p[x@from[l]]

    C <- c(C, x@from[l]); l <- l - 1
  }

  if (length(C) > 0) {
    p <- p[-C]; for (i in 1:d) pi[[i]] <- pi[[i]][-C, ]
  }

  rm(list = ls()[!(ls() %in% c("pi"))])

  pi
}) ## a.pi

setMethod("a.P",
          signature(x = "RCLRMIX"),
function(x, s)
{
  p <- x@p; pi <- x@pi; P <- x@P

  c <- x@c; d <- length(x@pi); s <- eval(s)

  if (!is.wholenumber(s)) {
    stop(sQuote("s"), " integer is requested!", call. = FALSE)
  }

  length(s) <- 1

  if ((s < 1) || (s > c)) {
    stop(sQuote("s"), " must be greater than 0 and less or equal than ", c, "!", call. = FALSE)
  }

  l <- c - 1; C <- numeric()

  while (l + 1 > s) {
    for (i in 1:d) {
      pi[[i]][x@to[l], ] <- (p[x@to[l]] * pi[[i]][x@to[l], ] + p[x@from[l]] * pi[[i]][x@from[l], ]) / (p[x@to[l]] + p[x@from[l]])
    }

    p[x@to[l]] <- p[x@to[l]] + p[x@from[l]]

    C <- c(C, x@from[l]); l <- l - 1
  }

  if (length(C) > 0) {
    p <- p[-C]; for (i in 1:d) pi[[i]] <- as.matrix(pi[[i]][-C, ])
  }

  dataset <- as.matrix(x@x@Dataset[[x@pos]])

  n <- nrow(dataset)

  Y <- dataset; y <- as.matrix(x@P[, 1:d]); Np <- array()

  for (j in 1:nrow(y)) {
    Np[j] <- 0.0

    for (l in 1:s) {
      Pl <- 1.0

      for(i in 1:d) {
        for (ii in 1:length(pi[[i]][l, ])) {
          if (y[j, i] == ii - 1) {
            Pl <- Pl * pi[[i]][l, ii]
          }
        }
      }

      Np[j] <- Np[j] + p[l] * Pl * n
    }
  }

  P[, d + 2] <- Np

  if (is.null(colnames(dataset))) {
    colnames(P) <- paste(c(1:d, "Nt", "Np"), sep = "")
  }
  else {
    colnames(P) <- c(colnames(dataset), "Nt", "Np")
  }

  rm(list = ls()[!(ls() %in% c("P"))])

  P
}) ## a.P

setMethod("a.tau",
          signature(x = "RCLRMIX"),
function(x, s)
{
  tau <- x@tau

  colnames <- colnames(x@tau)

  c <- x@c; s <- eval(s)

  if (!is.wholenumber(s)) {
    stop(sQuote("s"), " integer is requested!", call. = FALSE)
  }

  length(s) <- 1

  if ((s < 1) || (s > c)) {
    stop(sQuote("s"), " must be greater than 0 and less or equal than ", c, "!", call. = FALSE)
  }

  l <- c - 1; C <- numeric()

  while (l + 1 > s) {
    tau[, x@to[l]] <- tau[, x@to[l]] + tau[, x@from[l]]

    C <- c(C, x@from[l]); l <- l - 1
  }

  if (length(C) > 0) {
    tau <- as.matrix(tau[, -C]); colnames(tau) <- colnames[-C]
  }

  rm(list = ls()[!(ls() %in% c("tau"))])

  tau
}) ## a.tau

setMethod("a.prob", signature(x = "RCLRMIX"), function(x) x@prob)
setMethod("a.from", signature(x = "RCLRMIX"), function(x) x@from)
setMethod("a.to", signature(x = "RCLRMIX"), function(x) x@to)
setMethod("a.EN", signature(x = "RCLRMIX"), function(x) x@EN)
setMethod("a.ED", signature(x = "RCLRMIX"), function(x) x@ED)

setMethod("a.s", signature(x = "RCLS.chunk"), function(x) x@s)
setMethod("a.levels", signature(x = "RCLS.chunk"), function(x) x@levels)
setMethod("a.ntrain", signature(x = "RCLS.chunk"), function(x) x@ntrain)
setMethod("a.train", signature(x = "RCLS.chunk"), function(x) x@train)
setMethod("a.Zr", signature(x = "RCLS.chunk"), function(x) x@Zr)
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
