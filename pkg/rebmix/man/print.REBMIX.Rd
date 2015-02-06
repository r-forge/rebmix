\name{print.REBMIX}
\alias{print.REBMIX}
\title{
Prints Univariate or Multivariate REBMIX Output
}
\description{
Returns the \code{\link[base]{print}} output at \code{pos} for class REBMIX.
}
\usage{
\method{print}{REBMIX}(x, pos = 1, ...)
}
\arguments{
\item{x}{
an object of class \code{REBMIX}.
}
\item{pos}{
a desired row number in \code{summary} to be printed. The default value is \code{1}.
}
\item{\dots}{
further arguments to \code{\link[base]{print}}.
}
}
\value{
\code{print.REBMIX} returns (invisibly) the full value of x itself.
}
\examples{
## Generate and print simulated dataset.

Theta <- rbind(pdf = rep("binomial", 3),
  theta1 = c(10, 10, 10),
  theta2 = c(0.1, 0.4, 0.9))

simulated <- RNGMIX(Dataset = "simulated",
  rseed = -1,
  n = c(40, 60, 50),
  Theta = Theta)
  
print(simulated)

## Estimate number of components, component weights and component parameters.

simulatedest <- REBMIX(Dataset = simulated$Dataset, 
  Preprocessing = "histogram", 
  cmax = 6, 
  Criterion = "AIC", 
  Variables = "discrete",
  pdf = "binomial", 
  Theta1 = 10,
  K = 1)
  
## Print and plot finite mixture.  

print(simulatedest)

plot(simulatedest)
}
\keyword{print}

