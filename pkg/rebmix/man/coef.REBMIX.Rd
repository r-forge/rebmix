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
# Generate and print simulated dataset.

Theta <- list(pdf1 = rep("Weibull", 3),
  theta1.1 = c(2, 10, 30),
  theta2.1 = c(2, 4, 7))

simulated <- RNGMIX(Dataset = "simulated",
  rseed = -1,
  n = c(40, 60, 50),
  Theta = Theta)
  
print(simulated)

# Estimate number of components, component weights and component parameters.

simulatedest <- REBMIX(Dataset = simulated$Dataset, 
  Preprocessing = "histogram", 
  cmax = 6, 
  Criterion = "AIC", 
  Variables = "continuous",
  pdf = "Weibull", 
  K = 8:25)
  
# Print coefficients and plot finite mixture.  

coef(simulatedest)

plot(simulatedest)
}
\keyword{print}

