predict.list <- function(object,
  P = NULL,
  Dataset = NULL, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  if (missing(object)) {
    stop(sQuote("object"), " object or list of classes REBMIX is requested!", call. = FALSE)
  }

  if (class(object) == "list") {
    o <- length(object)

    for (i in 1:o) {
      if (class(object[[i]]) != "REBMIX") {
        stop(sQuote("object"), " list of classes REBMIX is requested!", call. = FALSE)
      }
    }
  }
  else
  if (class(object) == "REBMIX") {
    o <- 1;

    object[[1]] <- object
  }
  else {
    stop(sQuote("object"), " object or list of classes REBMIX is requested!", call. = FALSE)
  }

  if (is.null(Dataset)) {
    stop(sQuote("Dataset"), " must not be NULL!", call. = FALSE)
  }

  Dataset <- as.data.frame(Dataset)

  s <- nrow(object[[1]]$summary)

  if (o > 1) {
    for (i in 2:o) {
      if (s != nrow(object[[i]]$summary)) {
        stop(sQuote("object"), " list of classes REBMIX with equal number of classes is requested!", call. = FALSE)
      }
    }
  }

  if (s > 1) {
    if (is.null(P)) {
      stop(sQuote("P"), " must not be NULL!", call. = FALSE)
    }

    if (s != length(P)) {
      stop(sQuote("object"), " and ", sQuote("P"), " must be of the same length!", call. = FALSE)
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
      nrow <- nrow(object[[io]]$Theta[[is]])
      ncol <- ncol(object[[io]]$Theta[[is]])

      c[[io, is]] <- ncol

      w[[io, is]] <- as.numeric(object[[io]]$w[[is]])

      pdf[[io, is]] <- array(data = NA, dim = c(nrow, ncol), dimnames = NULL)
      theta1[[io, is]] <- array(data = 0.0, dim = c(nrow, ncol), dimnames = NULL)
      theta2[[io, is]] <- array(data = 0.0, dim = c(nrow, ncol), dimnames = NULL)

      for (j in 1:ncol) {
        M <- match(object[[io]]$Theta[[is]][, j], .rebmix$pdf)

        d[io] <- 1;

        for (l in 1:length(M)) {
          if (M[l] %in% which(.rebmix$pdf.nargs == 2)) {
            pdf[[io, is]][d[io], j] <- object[[io]]$Theta[[is]][l, j]
            theta1[[io, is]][d[io], j] <- as.numeric(object[[io]]$Theta[[is]][l + 1, j])
            theta2[[io, is]][d[io], j] <- as.numeric(object[[io]]$Theta[[is]][l + 2, j])

            d[io] <- d[io] + 1
          }
          else
          if (M[l] %in% which(.rebmix$pdf.nargs == 1)) {
            pdf[[io, is]][d[io], j] <- object[[io]]$Theta[[is]][l, j]
            theta1[[io, is]][d[io], j] <- as.numeric(object[[io]]$Theta[[is]][l + 1, j])

            d[io] <- d[io] + 1
          }
        }
      }

      d[io] <- d[io] - 1

      pdf[[io, is]] <- pdf[[io, is]][1:d[io], ]; dim(pdf[[io, is]]) <- c(d[io], ncol)
      theta1[[io, is]] <- theta1[[io, is]][1:d[io], ]; dim(theta1[[io, is]]) <- c(d[io], ncol)
      theta2[[io, is]] <- theta2[[io, is]][1:d[io], ]; dim(theta2[[io, is]]) <- c(d[io], ncol)
    }
  }

  if (s > 1) {
    message("RCLSMIX Version 2.4.1");
    flush.console()

    output <- .C("RCLSMIX",
      n = as.integer(nrow(Dataset)),
      X = as.double(unlist(Dataset)),
      s = as.integer(s),
      o = as.integer(o),
      d = as.integer(d),
      c = as.integer(unlist(c)),
      W = as.double(unlist(w)),
      ParFamType = as.character(unlist(pdf)),
      Par0 = as.double(unlist(theta1)),
      Par1 = as.double(unlist(theta2)),
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
} ## predict.list
