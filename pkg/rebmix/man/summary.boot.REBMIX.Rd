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
# Generate and print simulated dataset.

Theta <- list(pdf1 = rep("normal", 4),
  theta1.1 = c(1, 3, 6, 9),
  theta2.1 = c(0.1, 0.2, 0.3, 0.4))

n <- c(100, 200, 300, 400)

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
  cmax = 8, 
  Criterion = "BIC", 
  pdf = "normal", 
  K = K[1]:K[2])

plot(simulatedest )

# Bootstrap finite mixture and print summary.  

simulatedboot <- boot.REBMIX(x = simulatedest, B = 20)

summary(simulatedboot)
}
\keyword{print}

