\name{coef.REBMIX}
\alias{coef.REBMIX}
\title{
Prints Univariate or Multivariate REBMIX Coefficients
}
\description{
Returns the \code{w} and \code{Theta} output at \code{pos} for class REBMIX.
}
\usage{
\method{coef}{REBMIX}(object, pos = 1, ...)
}
\arguments{
\item{object}{
an object of class \code{REBMIX}.
}
\item{pos}{
a desired row number in \code{summary} to be printed. The default value is \code{1}.
}
\item{\dots}{
further arguments to \code{\link[base]{print}}.
}
}
\examples{
## Generate and print simulated dataset.

Theta <- rbind(pdf = rep("Weibull", 3),
  theta1 = c(2.0, 10.0, 30),
  theta2 = c(2.0, 4.0, 7.0))

simulated <- RNGMIX(Dataset = "simulated",
  rseed = -1,
  n = c(40, 60, 50),
  Theta = Theta)
  
simulated

## Estimate number of components, component weights and component parameters.

simulatedest <- REBMIX(Dataset = simulated$Dataset, 
  Preprocessing = "histogram", 
  D = 0.025, 
  cmax = 6, 
  Criterion = "AIC", 
  Variables = "continuous",
  pdf = "Weibull", 
  K = 8:25, 
  b = 1.0)
  
## Print coefficients and plot finite mixture.  

coef(simulatedest)

plot(simulatedest)
}
\keyword{print}

