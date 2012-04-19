\name{print.boot.REBMIX}
\alias{print.boot.REBMIX}
\title{
Prints Univariate or Multivariate boot.REBMIX Output
}
\description{
Returns the \code{\link[base]{print}} output for class boot.REBMIX.
}
\usage{
\method{print}{boot.REBMIX}(x, ...)
}
\arguments{
\item{x}{
an object of class \code{boot.REBMIX}.
}
\item{\dots}{
further arguments to \code{\link[base]{print}}.
}
}
\value{
\code{print.boot.REBMIX} returns (invisibly) the full value of x itself.
}
\examples{
## Generate and print simulated dataset.

Theta <- rbind(pdf = rep("Weibull", 2),
  theta1 = c(10, 20),
  theta2 = c(2.0, 6.0))

simulated <- RNGMIX(Dataset = "simulated",
  rseed = -1,
  n = c(30, 70),
  Theta = Theta)

## Estimate number of components, component weights and component parameters.

simulatedest <- REBMIX(Dataset = simulated$Dataset, 
  Preprocessing = "histogram", 
  D = 0.025, 
  cmax = 4, 
  Criterion = "AIC", 
  Variables = "continuous",
  pdf = "Weibull", 
  K = 20:40, 
  b = 1.0)
  
## Bootstrap finite mixture.  

simulatedboot <- boot.REBMIX(x = simulatedest, Bootstrap = "nonparametric", B = 20)

simulatedboot
}
\keyword{print}

