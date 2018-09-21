setMethod("summary",
          signature(object = "REBMIX"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }

  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object@summary), nomatch = 0)

  DF <- object@summary[p]

  is.num <- sapply(DF, is.number); DF[is.num] <- lapply(DF[is.num], as.number)

  print(DF, quote = FALSE, ...)

  cat(paste("Maximum logL = ", DF[object@pos, "logL"], " at pos = ", object@pos, ".\n", sep = "", collapse = ""))

  rm(list = ls())
}) ## summary

setMethod("summary",
          signature(object = "REBMVNORM"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMVNORM is requested!", call. = FALSE)
  }

  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object@summary), nomatch = 0)

  DF <- object@summary[p]

  is.num <- sapply(DF, is.number); DF[is.num] <- lapply(DF[is.num], as.number)

  print(DF, quote = FALSE, ...)

  cat(paste("Maximum logL = ", DF[object@pos, "logL"], " at pos = ", object@pos, ".\n", sep = "", collapse = ""))

  rm(list = ls())
}) ## summary

setMethod("summary",
          signature(object = "REBMIX.boot"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMIX.boot is requested!", call. = FALSE)
  }

  w.cv <- object@w.cv

  cat("w.cv", "\n", sep = "")

  print(w.cv, quote = FALSE, ...)

  d <- length(object@x@Variables)

  names <- names(object@Theta.cv)

  Names <- names[grep("theta1", names, fixed = TRUE)]

  theta1.cv <- matrix(unlist(object@Theta.cv[Names]), ncol = d, byrow = TRUE)

  rownames(theta1.cv) <- Names
  colnames(theta1.cv) <- NULL

  print(theta1.cv, quote = FALSE, ...)

  Names <- names[grep("theta2", names, fixed = TRUE)]

  theta2.cv <- matrix(unlist(object@Theta.cv[Names]), ncol = d, byrow = TRUE)

  rownames(theta2.cv) <- Names
  colnames(theta2.cv) <- NULL

  print(theta2.cv, quote = FALSE, ...)

  cat(paste("Mode probability = ", as.number(object@c.prob), " at c = ", object@c.mode, " components.\n", sep = "", collapse = ""))

  rm(list = ls())
}) ## summary

setMethod("summary",
          signature(object = "REBMVNORM.boot"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class REBMVNORM.boot is requested!", call. = FALSE)
  }

  w.cv <- object@w.cv
  
  cat("w.cv", "\n", sep = "")

  print(w.cv, quote = FALSE, ...)

  d <- length(object@x@Variables)

  names <- names(object@Theta.cv)

  Names <- names[grep("theta1", names, fixed = TRUE)]

  theta1.cv <- matrix(unlist(object@Theta.cv[Names]), ncol = d, byrow = TRUE)

  rownames(theta1.cv) <- Names
  colnames(theta1.cv) <- NULL

  print(theta1.cv, quote = FALSE, ...)

  Names <- names[grep("theta2", names, fixed = TRUE)]

  theta2.cv <- matrix(unlist(object@Theta.cv[Names]), ncol = d * d, byrow = TRUE)

  rownames(theta2.cv) <- Names
  colnames(theta2.cv) <- NULL

  print(theta2.cv, quote = FALSE, ...)

  cat(paste("Mode probability = ", as.number(object@c.prob), " at c = ", object@c.mode, " components.\n", sep = "", collapse = ""))

  rm(list = ls())
}) ## summary

setMethod("summary",
          signature(object = "RCLRMIX"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLRMIX is requested!", call. = FALSE)
  }

  c <- object@c - 1

  if (c > 0) {
    combine <- rbind(as.number(1:c), as.number(object@from), as.number(object@to), as.number(object@EN), as.number(object@ED))

    combine <- gsub(" ", "", combine, fixed = TRUE)

    rownames(combine) <- c("Number of clusters", "From cluster", "To cluster", "Entropy", "Entropy decrease")
    colnames(combine) <- rep("", c)

    combine <- capture.output(combine)[-1]
    combine <- paste(combine, "\n", sep = "")
    combine <- gsub("\"", " ", combine, fixed = TRUE)

    cat(combine, sep = "")
  }
  else {
    combine <- rbind(1)

    rownames(combine) <- "Number of clusters"
    colnames(combine) <- ""

    combine <- capture.output(combine)[-1]
    combine <- paste(combine, "\n", sep = "")
    combine <- gsub("\"", " ", combine, fixed = TRUE)

    cat(combine, sep = "")
  }

  rm(list = ls())
}) ## summary

setMethod("summary",
          signature(object = "RCLRMVNORM"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLRMVNORM is requested!", call. = FALSE)
  }

  c <- object@c - 1

  if (c > 0) {
    combine <- rbind(as.number(1:c), as.number(object@from), as.number(object@to), as.number(object@EN), as.number(object@ED))

    combine <- gsub(" ", "", combine, fixed = TRUE)

    rownames(combine) <- c("Number of clusters", "From cluster", "To cluster", "Entropy", "Entropy decrease")
    colnames(combine) <- rep("", c)

    combine <- capture.output(combine)[-1]
    combine <- paste(combine, "\n", sep = "")
    combine <- gsub("\"", " ", combine, fixed = TRUE)

    cat(combine, sep = "")
  }
  else {
    combine <- rbind(1)

    rownames(combine) <- "Number of clusters"
    colnames(combine) <- ""

    combine <- capture.output(combine)[-1]
    combine <- paste(combine, "\n", sep = "")
    combine <- gsub("\"", " ", combine, fixed = TRUE)

    cat(combine, sep = "")
  }

  rm(list = ls())
}) ## summary

setMethod("summary",
          signature(object = "RCLSMIX"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLSMIX is requested!", call. = FALSE)
  }

  CM <- as.data.frame(object@CM)

  colnames(CM) <- c("Test", "Predictive", "Frequency")

  print(CM, quote = FALSE, ...)

  cat(paste("Error = ", as.number(object@Error), ".\n", sep = "", collapse = ""))

  rm(list = ls())
}) ## summary

setMethod("summary",
          signature(object = "RCLSMVNORM"),
function(object, ...)
{
  if (missing(object)) {
    stop(sQuote("object"), " object of class RCLSMVNORM is requested!", call. = FALSE)
  }

  CM <- as.data.frame(object@CM)

  colnames(CM) <- c("Test", "Predictive", "Frequency")

  print(CM, quote = FALSE, ...)

  cat(paste("Error = ", as.number(object@Error), ".\n", sep = "", collapse = ""))

  rm(list = ls())
}) ## summary
