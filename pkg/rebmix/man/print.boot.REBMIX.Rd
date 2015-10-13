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
# Generate and print simulated dataset.

Theta <- list(pdf1 = rep("Weibull", 2),
  theta1.1 = c(10, 100),
  theta2.1 = c(2, 20))
  
n <- c(400, 600)

simulated <- RNGMIX(Dataset = "simulated",
  rseed = -1,
  n = n,
  Theta = Theta)
  
# Number of classes or nearest neighbours to be processed.

K <- c(as.integer(1 + log2(sum(n))), # Minimum v follows the Sturges rule.
  as.integer(10 * log10(sum(n)))) # Maximum v follows the log10 rule.  

# Estimate number of components, component weights and component parameters.

simulatedest <- REBMIX(Dataset = simulated$Dataset, 
  Preprocessing = "histogram", 
  cmax = 4, 
  Criterion = "AIC", 
  Variables = "continuous",
  pdf = "Weibull", 
  K = K[1]:K[2])
  
# Bootstrap finite mixture.  

simulatedboot <- boot.REBMIX(x = simulatedest, Bootstrap = "nonparametric", B = 20)

print(simulatedboot)
}
\keyword{print}

