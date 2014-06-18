\name{summary.REBMIX}
\alias{summary.REBMIX}
\title{
Prints Univariate or Multivariate REBMIX Summary
}
\description{
Returns the \code{summary} output for class REBMIX.
}
\usage{
\method{summary}{REBMIX}(object, ...)
}
\arguments{
\item{object}{
an object of class \code{REBMIX}.
}
\item{\dots}{
further arguments to \code{\link[base]{print}}.
}
}
\examples{
## Generate and print simulated dataset.

Theta <- rbind(pdf = rep("Dirac", 4),
  theta1 = c(1, 3, 4, 8))

simulated <- RNGMIX(Dataset = "simulated",
  rseed = -1,
  n = c(40, 20, 50, 10),
  Theta = Theta)
  
simulated

## Estimate number of components, component weights and component parameters.

simulatedest <- REBMIX(Dataset = simulated$Dataset, 
  Preprocessing = "histogram", 
  D = 0.025, 
  cmax = 8, 
  Criterion = "D", 
  Variables = "discrete",
  pdf = "Dirac", 
  K = 1)
  
## Print summary and plot finite mixture.  

summary(simulatedest)

plot(simulatedest)
}
\keyword{print}

