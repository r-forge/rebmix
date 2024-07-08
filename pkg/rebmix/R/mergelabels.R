setMethod("mergelabels",
          signature(A = "list"),
function(A, w, k, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  # A.
  
  if (missing(A) || (length(A) == 0)) {
    stop(sQuote("A"), " must not be empty!", call. = FALSE)
  }

  if (!is.list(A)) {
    stop(sQuote("A"), " list is requested!", call. = FALSE)
  }
  
  if (!all(unlist(lapply(A, is.matrix)))) {
    stop(sQuote("A"), " list of matrices is requested!", call. = FALSE)
  }
  
  nrow <- unique(unlist(lapply(A, nrow)))
  
  if (length(nrow) != 1) {
    stop(sQuote("A"), " numbers of rows in matrices must be equal!", call. = FALSE)
  }  
  
  ncol <- unique(unlist(lapply(A, ncol)))
  
  if (length(ncol) != 1) {
    stop(sQuote("A"), " numbers of columns in matrices must be equal!", call. = FALSE)
  }  

  if (nrow != ncol) {
    stop(sQuote("A"), " numbers of rows and columns in matrices must be equal!", call. = FALSE)
  }
    
  n <- length(A); c <- nrow
  
  set = rep(TRUE, c)

  for (i in 1:n) {
    for (j in 1:ncol) {
      tmp <- complete.cases(A[[i]][, j])
      
      if (any(tmp)) {
        set <- set & tmp; break
      }
    }  
  }
  
  nrow <- sum(set, na.rm = TRUE)
  
  for (i in 1:n) {
    A[[i]] <- A[[i]][set, ]
  }  
  
  set = rep(TRUE, c)

  for (i in 1:n) {
    for (j in 1:nrow) {
      tmp <- complete.cases(A[[i]][j, ])
      
      if (any(tmp)) {
        set <- set & tmp; break
      }
    }
  }
  
  ncol <- sum(set, na.rm = TRUE)
  
  for (i in 1:n) {
    A[[i]] <- A[[i]][, set]
  } 
  
  if (nrow != ncol) {
    stop(sQuote("A"), " numbers of rows and columns in matrices must be equal!", call. = FALSE)
  }
  
  c <- nrow
  
  if (c < 2) {
    stop(sQuote("c"), " must be greater than 1!", call. = FALSE)
  }    
  
  # w.
  
  if (length(w) != n) {
    stop(sQuote("w"), " number of weights must equal ", n, "!", call. = FALSE)
  }
  
  if (!is.number(w)) {
    stop(sQuote("w"), " numeric vector is requested!", call. = FALSE)
  }
  
  if (abs(sum(w) - 1.0) > 1.E-6){
    stop(sQuote("w"), " must sum to 1.0!", call. = FALSE) 
  }  
  
  # k.

  if (!is.wholenumber(k)) {
    stop(sQuote("k"), " integer is requested!", call. = FALSE)
  }

  if (k < 1) {
    stop(sQuote("k"), " must be greater than 0!", call. = FALSE)
  }
  
  output <- .C(C_RMergeLabels,
    n = as.integer(n),
    A = as.double(unlist(A)),
    c = as.integer(c),    
    w = as.double(w),
    L = double(c * c),
    error = character(3075),
    PACKAGE = "rebmix")

  error <- unlist(strsplit(output$error, "\n"));
      
  if (error[1] != "") {
    stop(error[1], call. = FALSE); return(NA)
  }
   
  if (error[2] != "") {
    warning(error[2], call. = FALSE, immediate. = TRUE)
  }  
    
  if (error[3] != "") {
    warning(error[3], call. = FALSE, immediate. = TRUE)
  } 
    
  output <- list(L = matrix(output$L, nrow = c, ncol = c))
  
  names <- paste(which(set), sep = "")
    
  rownames(output$L) <- names
  colnames(output$L) <- names
    
  k <- min(k, length(names))      
      
  x <- eigen(output$L, TRUE)$vectors[, 1:k]
    
  rownames(x) <- names
      
  output$cluster <- kmeans(x, k, ...)$cluster  
  
  options(digits = digits)
  
  rm(list = ls()[!(ls() %in% c("output"))])

  invisible(output)  
}) ## mergelabels
