setMethod("Zp",
          signature(x = "RCLRMIX"),
function(x, s)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class RCLRMIX is requested!", call. = FALSE)
  }

  Zp <- as.numeric(levels(x@Zp))[x@Zp]
  
  c <- x@c; s <- eval(s)
  
  if (!is.wholenumber(s)) {
    stop(sQuote("s"), " integer is requested!", call. = FALSE)
  }
  
  length(s) <- 1

  if ((s < 1) || (s > c)) {
    stop(sQuote("s"), " must be greater than 0 and less or equal than ", c, "!", call. = FALSE)
  }
 
  i <- c - 1
  
  while (s < length(unique(Zp))) {
    Zp[Zp == x@from[i]] <- x@to[i]
   
    i <- i - 1
  }  

  rm(list = ls()[!(ls() %in% c("Zp"))])
  
  as.factor(Zp)
}) ## Zp
  
setMethod("Zp",
          signature(x = "RCLRMVNORM"),
function(x, s)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class RCLRMVNORM is requested!", call. = FALSE)
  }
  
  Zp <- as.numeric(levels(x@Zp))[x@Zp]
  
  c <- x@c; s <- eval(s)
  
  if (!is.wholenumber(s)) {
    stop(sQuote("s"), " integer is requested!", call. = FALSE)
  }
  
  length(s) <- 1

  if ((s < 1) || (s > c)) {
    stop(sQuote("s"), " must be greater than 0 and less or equal than ", c, "!", call. = FALSE)
  }
 
  i <- c - 1
  
  while (s < length(unique(Zp))) {
    Zp[Zp == x@from[i]] <- x@to[i]
   
    i <- i - 1
  }  

  rm(list = ls()[!(ls() %in% c("Zp"))])
  
  as.factor(Zp)
}) ## Zp
