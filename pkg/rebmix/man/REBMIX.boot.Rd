\name{boot-methods}
\docType{methods}
\alias{boot}
\alias{boot-methods}
\alias{boot,REBMIX-method}
\alias{boot,REBMVNORM-method}
\alias{show,REBMIX.boot-method}
\alias{show,REBMVNORM.boot-method}
\alias{summary,REBMIX.boot-method}
\alias{summary,REBMVNORM.boot-method}
\title{
Parametric or Nonparametric Bootstrap for Standard Error and Coefficient of Variation Estimation
}
\description{
Returns as default the boot output for mixtures of conditionally independent normal,
lognormal, Weibull, gamma, binomial, Poisson, Dirac or von Mises component densities. If
\code{x} is of class \code{RNGMVNORM} the boot output for mixtures of multivariate normal
component densities with unrestricted variance-covariance matrices is returned.
}
\usage{
\S4method{boot}{REBMIX}(x = NULL, rseed = -1, pos = 1, Bootstrap = "parametric",
     B = 100, n = numeric(), replace = TRUE, prob = numeric(), ...)
## ... and for other signatures
\S4method{summary}{REBMIX.boot}(object, ...)
## ... and for other signatures
}
\arguments{
\item{x}{
see Methods section below.
}
\item{rseed}{
set the random seed to any negative integer value to initialize the sequence. The first bootstrap dataset corresponds to it.
For each next bootstrap dataset the random seed is decremented \eqn{r_{\mathrm{seed}} = r_{\mathrm{seed}} - 1}. The default value is \code{-1}.
}
\item{pos}{
a desired row number in \code{x@summary} to be bootstrapped. The default value is \code{1}.
}
\item{Bootstrap}{
a character giving the bootstrap type. One of default \code{"parametric"} or \code{"nonparametric"}.
}
\item{B}{
number of bootstrap datasets. The default value is \code{100}.
}
\item{n}{
number of observations. The default value is \code{numeric()}.
}
\item{replace}{
logical. The sampling is with replacement if \code{TRUE}, see also \code{\link[base]{sample}}. The default value is \code{TRUE}.
}
\item{prob}{
a vector of length \eqn{n} containing probability weights, see also \code{\link[base]{sample}}. The default value is \code{numeric()}.
}
\item{\dots}{
further arguments to \code{\link[base]{sample}}; additional arguments affecting the summary produced.
}
\item{object}{
see Methods section below.
}
}
\value{
Returns an object of class \code{REBMIX.boot} or \code{REBMVNORM.boot}.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "REBMIX")}}{an object of class \code{REBMIX} for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson, Dirac or von Mises component densities.}
\item{\code{signature(x = "REBMVNORM")}}{an object of class \code{REBMVNORM} for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices.}
\item{\code{signature(object = "REBMIX")}}{an object of class \code{REBMIX}.}
\item{\code{signature(object = "REBMVNORM")}}{an object of class \code{REBMVNORM}.}
}
}
\author{Marko Nagode}
\references{
G. McLachlan and D. Peel. Finite Mixture Models. John Wiley & Sons, New York, 2000.
}
\examples{
data("weibull")

n <- nrow(weibull)

# Number of classes or nearest neighbours to be processed.

K <- c(as.integer(1 + log2(n)), # Minimum v follows Sturges rule.
  as.integer(10 * log10(n))) # Maximum v follows log10 rule.

# Estimate number of components, component weights and component parameters.

weibullest <- REBMIX(Dataset = list(weibull),
  Preprocessing = "kernel density estimation",
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
