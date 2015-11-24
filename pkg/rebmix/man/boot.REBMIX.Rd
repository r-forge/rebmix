\name{boot}
\alias{boot}
\title{
Parametric or Nonparametric Bootstrap for Standard Error and Coefficient of Variation Estimation
}
\description{
Returns the boot.REBMIX output for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson or Dirac component densities.
}
\usage{
boot.REBMIX(x, pos = 1, Bootstrap = "parametric", 
     B = 100, n = NULL, replace = TRUE, prob = NULL, ...)
}
\arguments{
\item{x}{
an object of class \code{REBMIX}.
}
\item{pos}{
a desired row number in \code{x$summary} to be bootstrapped. The default value is \code{1}.
}
\item{Bootstrap}{
a character vector giving the bootstrap type. One of default \code{"parametric"} or \code{"nonparametric"}.
}
\item{B}{
number of bootstrap datasets. The default value is \code{100}.
}
\item{n}{
number of observations. The default value is \code{NULL}.
}
\item{replace}{
logical. The sampling is with replacement if \code{TRUE}, see also \code{\link[base]{sample}}. The default value is \code{TRUE}.
}
\item{prob}{
a vector of length \eqn{n} containing probability weights, see also \code{\link[base]{sample}}. The default value is \code{NULL}.
}
\item{\dots}{
further arguments to \code{\link[base]{sample}}.
}
}
\value{
Returns an object of class \code{REBMIX.boot} or \code{REBMVNORM.boot}.
}
\author{Marko Nagode}
\references{
G. McLachlan and D. Peel. Finite Mixture Models. John Wiley & Sons, New York, 2000.
}
\examples{
data("weibull")

n <- nrow(weibull)

# Number of classes or nearest neighbours to be processed.

K <- c(as.integer(1 + log2(sum(n))), # Minimum v follows the Sturges rule.
  as.integer(10 * log10(n))) # Maximum v follows the log10 rule.

# Estimate number of components, component weights and component parameters. 

weibullest <- REBMIX(Dataset = list(weibull),
  Preprocessing = "Parzen window", 
  cmax = 4,
  Criterion = "BIC", 
  pdf = "Weibull", 
  K = K[1]:K[2], 
  Restraints = "loose")
  
# Plot finite mixture.

plot(weibullest, what = c("density", "distribution", "IC", "logL", "D"), 
  nrow = 3, ncol = 2, npts = 1000)
  
# Bootstrap finite mixture.  

weibullboot <- boot(x = weibullest, Bootstrap = "nonparametric", B = 10)

weibullboot
}
\keyword{bootstrap}
