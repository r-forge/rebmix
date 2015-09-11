RCLSMIX <- function(x,
  P = NULL,
  Dataset = NULL, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  if (missing(x)) {
    stop(sQuote("x"), " object or list of classes REBMIX is requested!", call. = FALSE)
  }

  if (class(x) == "list") {
    o <- length(x)

    for (i in 1:o) {
      if (class(x[[i]]) != "REBMIX") {
        stop(sQuote("x"), " list of classes REBMIX is requested!", call. = FALSE)
      }
    }
  }
  else
  if (class(x) == "REBMIX") {
    o <- 1;

    x[[1]] <- x
  }
  else {
    stop(sQuote("x"), " object or list of classes REBMIX is requested!", call. = FALSE)
  }

  if (is.null(Dataset)) {
    stop(sQuote("Dataset"), " must not be NULL!", call. = FALSE)
  }

  Dataset <- as.data.frame(Dataset)

  s <- nrow(x[[1]]$summary)

  if (o > 1) {
    for (i in 2:o) {
      if (s != nrow(x[[i]]$summary)) {
        stop(sQuote("x"), " list of classes REBMIX with equal number of classes is requested!", call. = FALSE)
      }
    }
  }

  if (s > 1) {
    if (is.null(P)) {
      stop(sQuote("P"), " must not be NULL!", call. = FALSE)
    }

    if (s != length(P)) {
      stop(sQuote("x"), " and ", sQuote("P"), " must be of the same length!", call. = FALSE)
    }
  }

  d <- array(data = NA, dim = o, dimnames = NULL)

  c <- array(data = list(NULL), dim = c(o, s), dimnames = NULL)

  w <- array(data = list(NULL), dim = c(o, s), dimnames = NULL)

  pdf <- array(data = list(NULL), dim = c(o, s), dimnames = NULL)

  theta1 <- array(data = list(NULL), dim = c(o, s), dimnames = NULL)

  theta2 <- array(data = list(NULL), dim = c(o, s), dimnames = NULL)

  for (io in 1:o) {
    for (is in 1:s) {
      Names <- names(x[[io]]$Theta[[is]])
    
      pdf[[io, is]] <- as.character(unlist(x[[io]]$Theta[[is]][grep("pdf", Names)]))
    
      theta1[[io, is]] <- as.numeric(unlist(x[[io]]$Theta[[is]][grep("theta1", Names)]))
      
      theta1[[io, is]][is.na(theta1[[io, is]])] <- 0

      theta2[[io, is]] <- as.numeric(unlist(x[[io]]$Theta[[is]][grep("theta2", Names)]))
      
      theta2[[io, is]][is.na(theta2[[io, is]])] <- 0

      c[[io, is]] <- length(x[[io]]$w[[is]])

      w[[io, is]] <- x[[io]]$w[[is]]
      
      d[io] <- length(pdf[[io, is]]) / c[[io, is]]
    }
  }

  if (s > 1) {
    message("RCLSMIX Version 2.7.3");
    flush.console()

    output <- .C("RCLSMIX",
      n = as.integer(nrow(Dataset)),
      X = as.double(unlist(Dataset)),
      s = as.integer(s),
      o = as.integer(o),
      d = as.integer(d),
      c = as.integer(unlist(c)),
      W = as.double(unlist(w)),
      pdf = as.character(unlist(pdf)),
      theta1 = as.double(unlist(theta1)),
      theta2 = as.double(unlist(theta2)),
      P = as.double(unlist(P)),
      Z = integer(nrow(Dataset)),
      error = integer(1),
      PACKAGE = "rebmix")

    if (output$error == 1) {
      stop("in RCLSMIX!", call. = FALSE); return(NA)
    }
  }

  output <- as.factor(output$Z)
  levels(output) <- 0:(length(P) -1)

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## RCLSMIX
