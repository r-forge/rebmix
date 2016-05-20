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
  
  output@levels <- levels(factor(Dataset[, class]))
  
  Dataset <- base::split(Dataset, Dataset[, class])

  output@s <- length(Dataset)

  for (i in 1:output@s) {
    n <- nrow(Dataset[[i]]) 

    output@ntrain[i] <- as.integer(n * p)

    sample <- sample.int(n = n, size = output@ntrain[i], ...)

    output@train[[i]] <- Dataset[[i]][sample,]

    output@test <- rbind(output@test, Dataset[[i]][-sample,])
  }

  output@ntest <- nrow(output@test)

  output@train <- lapply(output@train, function(x) x[, -class])

  output@Zt <- factor(output@test[, class])

  output@test <- output@test[, -class]

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
  
  output <- new("RCLS.chunk") 
                 
  type <- as.character(unique(Dataset[, p$type]))
                 
  if (!(p$train %in% type)) {
    stop(sQuote("p$train"), " should be one of ", paste(dQuote(type), collapse = ", "), "!", call. = FALSE)
  }
  
  if (!(p$test %in% type)) {
    stop(sQuote("p$test"), " should be one of ", paste(dQuote(type), collapse = ", "), "!", call. = FALSE)
  }                  
                 
  output@test <- subset(Dataset, subset = Dataset[, p$type] == p$test)
  
  output@ntest <- nrow(output@test)
  
  output@Zt <- factor(output@test[, class])

  output@test <- output@test[, c(-class, -p$type)]
  
  Dataset <- subset(Dataset, subset = Dataset[, p$type] == p$train)
  
  output@levels <- levels(factor(Dataset[, class]))
  
  Dataset <- base::split(Dataset, Dataset[, class])
  
  output@s <- length(Dataset)
  
  for (i in 1:output@s) {
    output@ntrain[i] <- nrow(Dataset[[i]])
    
    output@train[[i]] <- Dataset[[i]]
  }
  
  output@train <- lapply(output@train, function(x) x[, c(-class, -p$type)])  

  rm(list = ls()[!(ls() %in% c("output"))]) 

  invisible(output)
}) ## split
