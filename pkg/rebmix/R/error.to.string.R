error.to.string <- function(error)
{
  output <- rep("", 3)

  if (error[1] > 0) {
    output[1] <- paste0("File = ", .error.defaults$FileNames[error[3] + 1], "; line = ", error[2], "; code = ", .error.defaults$ErrorNames[error[1] + 1], ".") 
  }
  
  if (error[4] > 0) {
    output[2] <- paste0("File = ", .error.defaults$FileNames[error[6] + 1], "; line = ", error[5], "; code = ", .error.defaults$ErrorNames[error[4] + 1], ".") 
  }
  
  if (error[7] > 0) {
    output[3] <- paste0("File = ", .error.defaults$FileNames[error[9] + 1], "; line = ", error[8], "; code = ", .error.defaults$ErrorNames[error[7] + 1], ".") 
  }  

  output
} ## error.to.string
