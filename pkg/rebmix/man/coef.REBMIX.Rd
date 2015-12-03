\name{coef-methods}
\docType{methods}
\alias{coef}
\alias{coef-methods}
\alias{coef,REBMIX-method}
\alias{coef,REBMVNORM-method}
\title{
Prints Univariate or Multivariate REBMIX Coefficients
}
\description{
Returns the \code{w} and \code{Theta} output at \code{pos}.
}
\usage{
\S4method{coef}{REBMIX}(object = NULL, pos = 1, ...)
## ... and for other signatures
}
\arguments{
\item{object}{
see Methods section below.
}
\item{pos}{
a desired row number in \code{summary} to be printed. The default value is \code{1}.
}
\item{\dots}{
currently not used..
}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "REBMIX")}}{an object of class \code{REBMIX}.}
\item{\code{signature(x = "REBMVNORM")}}{an object of class \code{REBMVNORM}.}
}
}
\author{Marko Nagode}
\examples{
# Generate and print simulated dataset.

Theta <- list(pdf1 = "Weibull",
  theta1.1 = c(2, 10, 30),
  theta2.1 = c(2, 4, 7))

simulated <- RNGMIX(Dataset.name = paste("simulated_", 1:5, sep = ""),
  rseed = -1,
  n = c(40, 60, 50),
  Theta = Theta)
  
simulated

# Estimate number of components, component weights and component parameters.

simulatedest <- REBMIX(Dataset = simulated@Dataset, 
  Preprocessing = "histogram", 
  cmax = 6, 
  Criterion = "AIC", 
  pdf = "Weibull", 
  K = 8:25)
  
# Print coefficients and plot finite mixture.  

coef(simulatedest, pos = 3)

plot(simulatedest)
}
\keyword{coef}

