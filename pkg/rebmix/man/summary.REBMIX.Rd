\name{summary}
\alias{summary}
\title{
Prints Univariate or Multivariate REBMIX Summary
}
\description{
Returns the \code{summary} output for class REBMIX.
}
\usage{
\S4method{summary}{REBMIX}(object, ...)
}
\arguments{
\item{object}{
an object of class \code{REBMIX}.
}
\item{\dots}{
additional arguments affecting the summary produced.
}
}
\author{Marko Nagode}
\examples{
# Generate and print simulated dataset.

Theta <- list(pdf1 = "Dirac",
  theta1.1 = c(1, 3, 4, 8),
  theta2.1 = c(NA, NA, NA, NA))

simulated <- RNGMIX(Dataset.name = "simulated",
  rseed = -1,
  n = c(40, 20, 50, 10),
  Theta = Theta)
  
show(simulated)

# Estimate number of components, component weights and component parameters.

simulatedest <- REBMIX(Dataset = simulated@Dataset, 
  Preprocessing = "histogram", 
  cmax = 8, 
  Criterion = "D", 
  pdf = "Dirac", 
  K = 1)
  
# Print summary and plot finite mixture.  

summary(simulatedest)

plot(simulatedest)
}
\keyword{summary}

