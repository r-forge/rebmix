setMethod("split",
          signature(p = "numeric"),
function(p = 0.75, Dataset, class, ...)
{
  if (missing(Dataset) || (length(Dataset) == 0)) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }

  if (!is.data.frame(Dataset)) {
    stop(sQuote("Dataset"), " data frame is requested!", call. = FALSE)
  }

  if (ncol(Dataset) < 2) {
    stop(sQuote("Dataset"), " number of columns in data frame must be greater than 1!", call. = FALSE)
  }

  if (!is.wholenumber(class)) {
    stop(sQuote("class"), " integer is requested!", call. = FALSE)
  }

  length(class) <- 1

  if ((class < 1) || (class > ncol(Dataset))) {
    stop(sQuote("class"), " must be greater than 0 and less or equal than ", ncol(Dataset), "!", call. = FALSE)
  }

  if ((p < 0.0) || (p > 1.0)) {
    stop(sQuote("p"), " must be greater or equal than 0.0 and less or equal than 1.0!", call. = FALSE)
  }

  output <- new("RCLS.chunk")

  names <- colnames(Dataset)

  f <- as.factor(Dataset[, class])

  output@levels <- levels(f)

  levels <- output@levels

  output@train <- Dataset <- base::split(Dataset, f)

  output@s <- length(Dataset)

  r <- integer()

  for (i in 1:output@s) {
    n <- nrow(Dataset[[i]])

    output@ntrain[i] <- as.integer(n * p)

    if (output@ntrain[i] > 1) {
      sample <- sample.int(n = n, size = output@ntrain[i])

      output@train[[i]] <- output@train[[i]][sample,]

      output@test <- rbind(output@test, Dataset[[i]][-sample,])
    }
    else {
      r <- c(r, i)

      output@test <- rbind(output@test, Dataset[[i]])
    }
  }

  if (length(r) > 0) {
    output@s <- output@s - length(r)

    output@levels <- output@levels[-r]

    output@train <- output@train[-r]

    output@ntrain <- output@ntrain[-r]

    levels <- c(levels[-r], levels[r])

    message("Note: Test dataset contains more classes than train datasets!")
  }

  output@Zr <- lapply(output@train, function(x) {x <- factor(x[, class], levels = output@levels); x})

  output@train <- lapply(output@train, function(x) {x <- data.frame(x[, -class]); colnames(x) <- names[-class]; x})
  
  output@Zt <- factor(output@test[, class], levels = levels)

  output@test <- data.frame(output@test[, -class])

  colnames(output@test) <- names[-class]

  output@ntest <- nrow(output@test)

  if (output@s < 1) {
    stop(sQuote("ntrain"), " number of observations in at least one train dataset must be greater than 1!", call. = FALSE)
  }

  if (output@ntest < 1) {
    stop(sQuote("ntest"), " number of observations in test dataset must be greater than 0!", call. = FALSE)
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)
}) ## split

setMethod("split",
          signature(p = "list"),
function(p = list(), Dataset, class, ...)
{
  if (missing(Dataset) || (length(Dataset) == 0)) {
    stop(sQuote("Dataset"), " must not be empty!", call. = FALSE)
  }

  if (!is.data.frame(Dataset)) {
    stop(sQuote("Dataset"), " data frame is requested!", call. = FALSE)
  }

  if (ncol(Dataset) < 3) {
    stop(sQuote("Dataset"), " number of columns in data frame must be greater than 2!", call. = FALSE)
  }

  if (!is.wholenumber(class)) {
    stop(sQuote("class"), " integer is requested!", call. = FALSE)
  }

  length(class) <- 1

  if ((class < 1) || (class > ncol(Dataset))) {
    stop(sQuote("class"), " must be greater than 0 and less or equal than ", ncol(Dataset), "!", call. = FALSE)
  }

  if (length(p) != 3) {
    stop("length of ",  sQuote("p"), " must equal 3!", call. = FALSE)
  }

  names(p) <- c("type", "train", "test")

  if (!is.wholenumber(p$type)) {
    stop(sQuote("p$type"), " integer is requested!", call. = FALSE)
  }

  length(p$type) <- 1

  if ((p$type < 1) || (p$type > ncol(Dataset))) {
    stop(sQuote("p$type"), " must be greater than 0 and less or equal than ", ncol(Dataset), "!", call. = FALSE)
  }

  type <- as.character(unique(Dataset[, p$type]))

  if (!(p$train %in% type)) {
    stop(sQuote("p$train"), " should be one of ", paste(dQuote(type), collapse = ", "), "!", call. = FALSE)
  }

  if (!(p$test %in% type)) {
    stop(sQuote("p$test"), " should be one of ", paste(dQuote(type), collapse = ", "), "!", call. = FALSE)
  }
  
  output <- new("RCLS.chunk") 

  names <- colnames(Dataset)
  
  output@test <- subset(Dataset, subset = Dataset[, p$type] == p$test)
  
  output@ntest <- nrow(output@test)
  
  f <- droplevels(as.factor(output@test[, class]))
  
  levels.test <- levels(f)
  
  output@train <- subset(Dataset, subset = Dataset[, p$type] == p$train)

  f <- droplevels(as.factor(output@train[, class]))

  r <- f %in% levels.test

  output@train <- output@train[r, ]

  f <- droplevels(f[r])

  output@train <- base::split(output@train, f)

  levels.train <- levels(f)

  output@s <- length(output@train)  

  r <- integer()  

  for (i in 1:output@s) {
    output@ntrain[i] <- nrow(output@train[[i]])

    if (output@ntrain[i] <= 1) {
      r <- c(r, i)
    }
  }

  if (length(r) > 0) {
    output@s <- output@s - length(r)

    levels.train <- levels.train[-r]

    output@train <- output@train[-r]

    output@ntrain <- output@ntrain[-r]
  }

  r <- levels.test %in% levels.train
  
  if (any(r == FALSE)) {
    message("Note: Test dataset contains more classes than train datasets!")
  }

  levels.test <- c(levels.test[r], levels.test[!r])

  output@Zr <- lapply(output@train, function(x) {x <- factor(x[, class], levels = levels.train); x})

  output@train <- lapply(output@train, function(x) {x <- data.frame(x[, c(-class, -p$type)]); colnames(x) <- names[c(-class, -p$type)]; x})

  output@Zt <- factor(output@test[, class], levels = levels.test)

  output@test <- data.frame(output@test[, c(-class, -p$type)])
  
  output@levels <- levels.test

  colnames(output@test) <- names[c(-class, -p$type)]

  if (output@s < 1) {
    stop(sQuote("ntrain"), " number of observations in at least one train dataset must be greater than 1!", call. = FALSE)
  }

  if (output@ntest < 1) {
    stop(sQuote("ntest"), " number of observations in test dataset must be greater than 0!", call. = FALSE)
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)
}) ## split
