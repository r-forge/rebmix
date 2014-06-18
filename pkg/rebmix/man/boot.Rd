\name{boot}
\alias{boot}
\alias{boot.default}
\alias{boot.REBMIX}
\title{
Parametric or Nonparametric Bootstrap for Standard Error and Coefficient of Variation Estimation
}
\description{
Returns the boot.REBMIX output for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson or Dirac component densities.
}
\usage{
\method{boot}{REBMIX}(x, pos = 1, Bootstrap = "parametric", 
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
number of bootstrap samples. The default value is \code{100}.
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
\item{c}{
a data frame containing numbers of components for \eqn{B} bootstrap samples.
}
\item{c.se}{
standard error of numbers of components \code{c}.
}
\item{c.cv}{
coefficient of variation of numbers of components \code{c}.
}
\item{c.mode}{
mode of numbers of components \code{c}.
}
\item{c.prob}{
probability of mode \code{c.mode}.
}
\item{w}{
a data frame containing component weights for \eqn{\leq B} bootstrap samples.
}
\item{w.se}{
a vector containing standard errors of component weights \code{w}.
}
\item{w.cv}{
a vector containing coefficients of variation of component weights \code{w}.
}
\item{theta1.i}{
data frames containing component parameters \code{theta1.i} for \eqn{\leq B} bootstrap samples.
}
\item{theta1.i.se}{
vectors containing standard errors of component parameters \code{theta1.i}.
}
\item{theta1.i.cv}{
vectors containing coefficients of variation of component parameters \code{theta1.i}.
}
\item{theta2.i}{
data frames containing component parameters \code{theta2.i} for \eqn{\leq B} bootstrap samples.
}
\item{theta2.i.se}{
vectors containing standard errors of component parameters \code{theta2.i}.
}
\item{theta2.i.cv}{
vectors containing coefficients of variation of component parameters \code{theta2.i}.
}
}
\references{
McLachlan GJ, Peel D (2000). \emph{Finite Mixture Models}. John Wiley & Sons, New York.
}
\examples{
data("weibull")

n <- nrow(weibull)

## Number of classes or nearest neighbours to be processed.

K <- c(as.integer(1 + log2(sum(n))), ## Minimum v follows the Sturges rule.
  as.integer(10 * log10(n))) ## Maximum v follows the log10 rule.

## Estimate number of components, component weights and component parameters. 

weibullest <- REBMIX(Dataset = list(weibull),
  Preprocessing = "Parzen window", 
  D = 0.025, 
  cmax = 4,
  Criterion = "BIC", 
  Variables = "continuous",
  pdf = "Weibull", 
  K = K[1]:K[2], 
  Restraints = "loose")

## Plot finite mixture.

plot(weibullest, what = c("density", "distribution", "c", "IC", "logL", "D"), 
  nrow = 3, ncol = 2, npts = 1000)

## Bootstrap finite mixture.  

weibullboot <- boot(x = weibullest, Bootstrap = "nonparametric")

weibullboot
}
\keyword{bootstrap}
