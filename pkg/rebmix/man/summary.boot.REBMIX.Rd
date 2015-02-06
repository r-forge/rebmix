\name{summary.boot.REBMIX}
\alias{summary.boot.REBMIX}
\title{
Prints Univariate or Multivariate boot.REBMIX Summary
}
\description{
Returns the \code{summary} output for class boot.REBMIX.
}
\usage{
\method{summary}{boot.REBMIX}(object, ...)
}
\arguments{
\item{object}{
an object of class \code{boot.REBMIX}.
}
\item{\dots}{
further arguments to \code{\link[base]{print}}.
}
}
\examples{
## Generate and print simulated dataset.

Theta <- rbind(pdf = rep("normal", 4),
  theta1 = c(1, 3, 6, 9),
  theta2 = c(0.1, 0.2, 0.3, 0.4))

simulated <- RNGMIX(Dataset = "simulated",
  rseed = -1,
  n = c(10, 20, 30, 40),
  Theta = Theta)

## Estimate number of components, component weights and component parameters.

simulatedest <- REBMIX(Dataset = simulated$Dataset, 
  Preprocessing = "histogram", 
  cmax = 8, 
  Criterion = "AIC", 
  Variables = "continuous",
  pdf = "normal", 
  K = 30:40)

## Bootstrap finite mixture and print summary.  

simulatedboot <- boot.REBMIX(x = simulatedest, B = 20)

summary(simulatedboot)
}
\keyword{print}

