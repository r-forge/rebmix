setMethod("boot",
          signature(x = "REBMIX"),
function(x,
  rseed,
  pos,
  Bootstrap,
  B,
  n,
  replace,
  prob, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  model <- new("REBMIX.boot",
    x = x,
    rseed = rseed,
    pos = pos,
    Bootstrap = Bootstrap,
    B = B,
    n = n,
    replace = replace,
    prob = prob)

  if (model@Bootstrap == .rebmix.boot$Bootstrap[1]) {
    bsample <- RNGMIX(Dataset.name = paste("bsample_", 1:model@B, sep = ""), rseed = rseed,
      n = ceiling(model@n * as.numeric(model@x@w[[model@pos]])),
      Theta = model@x@Theta[[model@pos]])

    Dataset <- bsample@Dataset
  }
  else
  if (model@Bootstrap == .rebmix.boot$Bootstrap[2]) {
    dataset <- as.matrix(model@x@Dataset[[model@pos]])

    Dataset <- list()

    for (i in 1:model@B) {
      R1 <- sample.int(n = nrow(dataset), size = model@n, replace = replace, prob = if (length(prob) == 0) NULL else prob)

      Dataset[[i]] <- as.data.frame(dataset[R1,], stringsAsFactors = FALSE)
    }
  }

  if (length(model@x@theta1) > 0) {
    theta1 <- model@x@theta1; theta1[is.na(theta1)] <- 0
  }
  else {
    theta1 <- numeric()
  }

  if (length(model@x@theta2) > 0) {
    theta2 <- model@x@theta2; theta2[is.na(theta2)] <- 0
  }
  else {
    theta2 <- numeric()
  }
  
  if (length(model@x@theta3) > 0) {
    theta3 <- model@x@theta3; theta3[is.na(theta3)] <- 0
  }
  else {
    theta3 <- numeric()
  }  

  d <- length(model@x@Variables)
  
  # cmax.
  
  cmax <- list(...)$cmax

  if (is.null(cmax) || (length(cmax) == 0)) cmax <- as.numeric(model@x@summary[model@pos, "cmax"])

  if (!is.wholenumber(cmax)) {
    stop(sQuote("cmax"), " integer is requested!", call. = FALSE)
  }

  if (cmax < 1) {
    stop(sQuote("cmax"), " must be greater than 0!", call. = FALSE)
  }
  
  # cmin.
  
  cmin <- list(...)$cmin

  if (is.null(cmin) || (length(cmin) == 0)) cmin <- as.numeric(model@x@summary[model@pos, "cmin"])

  if (!is.wholenumber(cmin)) {
    stop(sQuote("cmin"), " integer is requested!", call. = FALSE)
  }

  if (cmin < 1) {
    stop(sQuote("cmin"), " must be greater than 0!", call. = FALSE)
  }

  if (cmin >= cmax) {
    stop(sQuote("cmax"), " must be greater than ", cmin, "!", call. = FALSE)
  }  

  bsampleest <- REBMIX(Dataset = Dataset,
    Preprocessing = as.character(model@x@summary[model@pos, "Preprocessing"]),
    cmax = cmax,
    cmin = cmin,
    Criterion = as.character(model@x@summary[model@pos, "Criterion"]),
    Variables = model@x@Variables,
    pdf = model@x@pdf,
    theta1 = theta1,
    theta2 = theta2,
    theta3 = theta3,
    K = eval(parse(text = as.character(model@x@summary[model@pos, "K"]))),
    y0 = model@x@y0,
    ymin = model@x@ymin,
    ymax = model@x@ymax,
    ar = as.numeric(model@x@summary[model@pos, "ar"]),
    Restraints = as.character(model@x@summary[model@pos, "Restraints"]))

  freq <- table(as.numeric(bsampleest@summary$c))

  c <- as.integer(names(freq)[which.max(freq)])

  model@c <- as.numeric(bsampleest@summary$c)
  model@c.se <- sd(model@c)
  model@c.cv <- model@c.se / mean(model@c)

  w <- bsampleest@w[model@c == c]

  Theta <- bsampleest@Theta[model@c == c]

  model@c.mode <- c
  model@c.prob <- length(w) / model@B

  model@w <- matrix(unlist(w), ncol = c, byrow = TRUE)

  colnames(model@w) <- paste(1:c, sep = "")
  rownames(model@w) <- paste(which(bsampleest@summary$c == c), sep = "")

  model@w.se <- as.vector(apply(model@w, 2, sd))
  model@w.cv <- as.vector(model@w.se / apply(model@w, 2, mean))

  for (i in 1:model@c.mode) {
    theta1 <- paste("theta1.",  i, sep = "")

    model@Theta[[theta1]] <- NULL

    for (j in 1:length(Theta)) {
      model@Theta[[theta1]] <- c(model@Theta[[theta1]], Theta[[j]][[theta1]])
    }

    model@Theta[[theta1]] <- matrix(model@Theta[[theta1]], ncol = d, byrow = TRUE)

    colnames(model@Theta[[theta1]]) <- NULL
    rownames(model@Theta[[theta1]]) <- paste(which(bsampleest@summary$c == c), sep = "")

    theta1.se <- paste("theta1.",  i, ".se", sep = "")
    theta1.cv <- paste("theta1.",  i, ".cv", sep = "")

    model@Theta.se[[theta1.se]] <- as.vector(apply(model@Theta[[theta1]], 2, sd))
    model@Theta.cv[[theta1.cv]] <- as.vector(model@Theta.se[[theta1.se]] / apply(model@Theta[[theta1]], 2, mean))
  }

  for (i in 1:model@c.mode) {
    theta2 <- paste("theta2.",  i, sep = "")

    model@Theta[[theta2]] <- NULL

    for (j in 1:length(Theta)) {
      model@Theta[[theta2]] <- c(model@Theta[[theta2]], Theta[[j]][[theta2]])
    }

    model@Theta[[theta2]] <- matrix(model@Theta[[theta2]], ncol = d, byrow = TRUE)

    colnames(model@Theta[[theta2]]) <- NULL
    rownames(model@Theta[[theta2]]) <- paste(which(bsampleest@summary$c == c), sep = "")

    theta2.se <- paste("theta2.",  i, ".se", sep = "")
    theta2.cv <- paste("theta2.",  i, ".cv", sep = "")

    model@Theta.se[[theta2.se]] <- as.vector(apply(model@Theta[[theta2]], 2, sd))
    model@Theta.cv[[theta2.cv]] <- as.vector(model@Theta.se[[theta2.se]] / apply(model@Theta[[theta2]], 2, mean))
  }
  
  for (i in 1:model@c.mode) {
    theta3 <- paste("theta3.",  i, sep = "")

    model@Theta[[theta3]] <- NULL

    for (j in 1:length(Theta)) {
      model@Theta[[theta3]] <- c(model@Theta[[theta3]], Theta[[j]][[theta3]])
    }

    model@Theta[[theta3]] <- matrix(model@Theta[[theta3]], ncol = d, byrow = TRUE)

    colnames(model@Theta[[theta3]]) <- NULL
    rownames(model@Theta[[theta3]]) <- paste(which(bsampleest@summary$c == c), sep = "")

    theta3.se <- paste("theta3.",  i, ".se", sep = "")
    theta3.cv <- paste("theta3.",  i, ".cv", sep = "")

    model@Theta.se[[theta3.se]] <- as.vector(apply(model@Theta[[theta3]], 2, sd))
    model@Theta.cv[[theta3.cv]] <- as.vector(model@Theta.se[[theta3.se]] / apply(model@Theta[[theta3]], 2, mean))
  }  

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## boot

setMethod("boot",
          signature(x = "REBMVNORM"),
function(x,
  rseed,
  pos,
  Bootstrap,
  B,
  n,
  replace,
  prob, ...)
{
  digits <- getOption("digits"); options(digits = 15)

  model <- new("REBMVNORM.boot",
    x = x,
    rseed = rseed,
    pos = pos,
    Bootstrap = Bootstrap,
    B = B,
    n = n,
    replace = replace,
    prob = prob)

  if (model@Bootstrap == .rebmix.boot$Bootstrap[1]) {
    bsample <- RNGMIX(model = "RNGMVNORM", Dataset.name = paste("bsample_", 1:model@B, sep = ""), rseed = rseed,
      n = ceiling(model@n * as.numeric(model@x@w[[model@pos]])),
      Theta = model@x@Theta[[model@pos]])

    Dataset <- bsample@Dataset
  }
  else
  if (model@Bootstrap == .rebmix.boot$Bootstrap[2]) {
    dataset <- as.matrix(model@x@Dataset[[model@pos]])

    Dataset <- list()

    for (i in 1:model@B) {
      R1 <- sample.int(n = nrow(dataset), size = model@n, replace = replace, prob = if (length(prob) == 0) NULL else prob)

      Dataset[[i]] <- as.data.frame(dataset[R1,], stringsAsFactors = FALSE)
    }
  }

  if (length(model@x@theta1) > 0) {
    theta1 <- model@x@theta1; theta1[is.na(theta1)] <- 0
  }
  else {
    theta1 <- numeric()
  }

  if (length(model@x@theta2) > 0) {
    theta2 <- model@x@theta2; theta2[is.na(theta2)] <- 0
  }
  else {
    theta2 <- numeric()
  }

  d <- length(model@x@Variables)

  # cmax.
  
  cmax <- list(...)$cmax

  if (is.null(cmax) || (length(cmax) == 0)) cmax <- as.numeric(model@x@summary[model@pos, "cmax"])

  if (!is.wholenumber(cmax)) {
    stop(sQuote("cmax"), " integer is requested!", call. = FALSE)
  }

  if (cmax < 1) {
    stop(sQuote("cmax"), " must be greater than 0!", call. = FALSE)
  }
  
  # cmin.
  
  cmin <- list(...)$cmin

  if (is.null(cmin) || (length(cmin) == 0)) cmin <- as.numeric(model@x@summary[model@pos, "cmin"])

  if (!is.wholenumber(cmin)) {
    stop(sQuote("cmin"), " integer is requested!", call. = FALSE)
  }

  if (cmin < 1) {
    stop(sQuote("cmin"), " must be greater than 0!", call. = FALSE)
  }

  if (cmin >= cmax) {
    stop(sQuote("cmax"), " must be greater than ", cmin, "!", call. = FALSE)
  }  

  bsampleest <- REBMIX(model = "REBMVNORM",
    Dataset = Dataset,
    Preprocessing = as.character(model@x@summary[model@pos, "Preprocessing"]),
    cmax = cmax,
    cmin = cmin,
    Criterion = as.character(model@x@summary[model@pos, "Criterion"]),
    Variables = model@x@Variables,
    pdf = model@x@pdf,
    theta1 = theta1,
    theta2 = theta2,
    K = eval(parse(text = as.character(model@x@summary[model@pos, "K"]))),
    y0 = model@x@y0,
    ymin = model@x@ymin,
    ymax = model@x@ymax,
    ar = as.numeric(model@x@summary[model@pos, "ar"]),
    Restraints = as.character(model@x@summary[model@pos, "Restraints"]))

  freq <- table(as.numeric(bsampleest@summary$c))

  c <- as.integer(names(freq)[which.max(freq)])

  model@c <- as.numeric(bsampleest@summary$c)
  model@c.se <- sd(model@c)
  model@c.cv <- model@c.se / mean(model@c)

  w <- bsampleest@w[model@c == c]

  Theta <- bsampleest@Theta[model@c == c]

  model@c.mode <- c
  model@c.prob <- length(w) / model@B

  model@w <- matrix(unlist(w), ncol = c, byrow = TRUE)

  colnames(model@w) <- paste(1:c, sep = "")
  rownames(model@w) <- paste(which(bsampleest@summary$c == c), sep = "")

  model@w.se <- as.vector(apply(model@w, 2, sd))
  model@w.cv <- as.vector(model@w.se / apply(model@w, 2, mean))

  for (i in 1:model@c.mode) {
    theta1 <- paste("theta1.",  i, sep = "")

    model@Theta[[theta1]] <- NULL

    for (j in 1:length(Theta)) {
      model@Theta[[theta1]] <- c(model@Theta[[theta1]], Theta[[j]][[theta1]])
    }

    model@Theta[[theta1]] <- matrix(model@Theta[[theta1]], ncol = d, byrow = TRUE)

    colnames(model@Theta[[theta1]]) <- NULL
    rownames(model@Theta[[theta1]]) <- paste(which(bsampleest@summary$c == c), sep = "")

    theta1.se <- paste("theta1.",  i, ".se", sep = "")
    theta1.cv <- paste("theta1.",  i, ".cv", sep = "")

    model@Theta.se[[theta1.se]] <- as.vector(apply(model@Theta[[theta1]], 2, sd))
    model@Theta.cv[[theta1.cv]] <- as.vector(model@Theta.se[[theta1.se]] / apply(model@Theta[[theta1]], 2, mean))
  }

  for (i in 1:model@c.mode) {
    theta2 <- paste("theta2.",  i, sep = "")

    model@Theta[[theta2]] <- NULL

    for (j in 1:length(Theta)) {
      model@Theta[[theta2]] <- c(model@Theta[[theta2]], Theta[[j]][[theta2]])
    }

    model@Theta[[theta2]] <- matrix(model@Theta[[theta2]], ncol = d * d, byrow = TRUE)

    colnames(model@Theta[[theta2]]) <- NULL
    rownames(model@Theta[[theta2]]) <- paste(which(bsampleest@summary$c == c), sep = "")

    theta2.se <- paste("theta2.",  i, ".se", sep = "")
    theta2.cv <- paste("theta2.",  i, ".cv", sep = "")

    model@Theta.se[[theta2.se]] <- as.vector(apply(model@Theta[[theta2]], 2, sd))
    model@Theta.cv[[theta2.cv]] <- as.vector(model@Theta.se[[theta2.se]] / apply(model@Theta[[theta2]], 2, mean))
  }

  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("model"))])

  return(model)
}) ## boot
