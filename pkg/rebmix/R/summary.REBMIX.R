summary.REBMIX <- function(object, ...)
{
  if (missing(object) || (class(object) != "REBMIX")) {
    stop(sQuote("object"), " object of class REBMIX is requested!", call. = FALSE)
  }

  p <- match(c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M"), names(object$summary), nomatch = 0)
  
  summary <- as.data.frame(rbind(object$summary[p]), stringsAsFactors = FALSE)

  rownames(summary) <- paste(1:nrow(summary), sep = "")  
  colnames(summary) <- c("Dataset", "Preprocessing", "Criterion", "c", "v/k", "IC", "logL", "M")

  print(apply(summary, c(1, 2), as.number), quote = FALSE, ...)

  cat(paste("Maximum logL = ", as.number(object$summary[object$pos, "logL"]), " at pos = ", as.number(object$pos), ".\n", sep = "", collapse = ""))
  
  rm(list = ls()) 
} ## summary.REBMIX
