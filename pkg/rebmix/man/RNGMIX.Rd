\name{RNGMIX-methods}
\docType{methods}
\alias{RNGMIX}
\alias{RNGMIX-methods}
\alias{RNGMIX,RNGMIX-method}
\alias{RNGMIX,RNGMVNORM-method}
\alias{show,RNGMIX-method}
\alias{show,RNGMVNORM-method}
\title{
Random Univariate or Multivariate Finite Mixture Generation
}
\description{
Returns as default the RNGMIX univariate or multivariate random datasets for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson or Dirac component densities.
If \code{model} equals \code{"RNGMVNORM"} multivariate random datasets for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices are returned.
}
\usage{
\S4method{RNGMIX}{RNGMIX}(model = "RNGMIX", Dataset.name = character(), 
       rseed = -1, n = numeric(), Theta = list(), ...)
## ... and for other signatures
}
\arguments{
\item{model}{
see Methods section below.
}
\item{Dataset.name}{
a character vector containing list names of data frames of size \eqn{n \times d} that \emph{d}-dimensional datasets are written in.
}
\item{rseed}{
set the random seed to any negative integer value to initialize the sequence. The first file in \code{Dataset.name} corresponds to it.
For each next file the random seed is decremented \eqn{r_{\mathrm{seed}} = r_{\mathrm{seed}} - 1}. The default value is \code{-1}.
}
\item{n}{
a vector containing numbers of observations in classes \eqn{n_{l}}, where number of observations \eqn{n = \sum_{l = 1}^{c} n_{l}}.
}
\item{Theta}{
a list containing \eqn{c} parametric family types \code{pdfi}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
Component parameters \code{theta1.i} follow the parametric family types. One of \eqn{\mu_{il}} for normal and lognormal distributions and \eqn{\theta_{il}} for Weibull, gamma, binomial, Poisson and Dirac distributions.
Component parameters \code{theta2.i} follow \code{theta1.i}. One of \eqn{\sigma_{il}} for normal and lognormal distributions, \eqn{\beta_{il}} for Weibull and gamma distributions and \eqn{p_{il}} for binomial distribution.
}
\item{\dots}{
currently not used.
}
}
\details{
RNGMIX is based on the "Minimal" random number generator \code{ran1} of Park and Miller with the Bays-Durham shuffle and added safeguards that returns a uniform random deviate between 0.0 and 1.0 
(exclusive of the endpoint values).
}
\value{
Returns an object of class \code{RNGMIX} or \code{RNGMVNORM}.
}
\section{Methods}{
\describe{
\item{\code{signature(model = "RNGMIX")}}{a character giving the default class name \code{"RNGMIX"} for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson or Dirac component densities.}
\item{\code{signature(model = "RNGMVNORM")}}{a character giving the class name \code{"RNGMVNORM"} for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices.}
}
}
\author{Marko Nagode}
\references{
W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery. Numerical Recipes in C: The Art
of Scientific Computing. Cambridge University Press, Cambridge, 1992.
}
\examples{
# Generate and print multivariate normal datasets with diagonal 
# variance-covariance matrices.

n <- c(75, 100, 125, 150, 175)

Theta <- list(pdf1 = rep("normal", 4),
  theta1.1 = c(10, 12, 10, 12),
  theta2.1 = c(1, 1, 1, 1),
  pdf2 = rep("normal", 4),
  theta1.2 = c(8.5, 10.5, 8.5, 10.5),
  theta2.2 = c(1, 1, 1, 1),
  pdf3 = rep("normal", 4),
  theta1.3 = c(12, 14, 12, 14),
  theta2.3 = c(1, 1, 1, 1),
  pdf4 = rep("normal", 4),
  theta1.4 = c(13, 15, 7, 9),
  theta2.4 = c(2, 2, 2, 2),
  pdf5 = rep("normal", 4),
  theta1.5 = c(7, 9, 13, 15),
  theta2.5 = c(3, 3, 3, 3))

simulated <- RNGMIX(Dataset.name = paste("simulated_", 1:25, sep = ""),
  rseed = -1,
  n = n,
  Theta = Theta)
  
simulated

# Generate and print multivariate normal datasets with unrestricted 
# variance-covariance matrices.

n <- c(200, 50, 50)

Theta <- list(pdf1 = rep("normal", 3),
  theta1.1 = c(0, 0, 0),
  theta2.1 = c(9, 0, 0, 0, 4, 0, 0, 0, 1),
  pdf2 = rep("normal", 3),
  theta1.2 = c(-6, 3, 6),
  theta2.2 = c(4, -3.2, -0.2, -3.2, 4, 0, -0.2, 0, 1),
  pdf3 = rep("normal", 3),
  theta1.3 = c(6, 6, 4),
  theta2.3 = c(4, 3.2, 2.8, 3.2, 4, 2.4, 2.8, 2.4, 2))

simulated <- RNGMIX(model = "RNGMVNORM",
  Dataset.name = paste("simulated_", 1:2, sep = ""),
  rseed = -1,
  n = n,
  Theta = Theta)
  
simulated

# Generate and print multivariate mixed continuous-discrete dataset.

n <- c(400, 100, 500)

Theta <- list(pdf1 = c("lognormal", "Poisson", "binomial", "Weibull"),
  theta1.1 = c(1, 2, 10, 2), 
  theta2.1 = c(0.3, NA, 0.9, 3),
  pdf2 = c("lognormal", "Poisson", "binomial", "Weibull"),
  theta1.2 = c(3.5, 10, 10, 10),
  theta2.2 = c(0.2, NA, 0.1, 7),
  pdf3 = c("lognormal", "Poisson", "binomial", "Weibull"),
  theta1.3 = c(2.5, 15, 10, 25), 
  theta2.3 = c(0.4, NA, 0.7, 20))
  
simulated <- RNGMIX(Dataset.name = paste("simulated_", 1:5, sep = ""),
  rseed = -1,
  n = n,
  Theta = Theta)
  
simulated  

# Generate and print univariate mixed Weibull dataset.

n <- c(75, 100, 125, 150, 175)

Theta <- list(pdf1 = "Weibull",
  theta1.1 = 12,
  theta2.1 = 2,
  pdf2 = "Weibull",
  theta1.2 = 10,
  theta2.2 = 4.1,  
  pdf3 = "Weibull",
  theta1.3 = 14,
  theta2.3 = 3.2,  
  pdf4 = "Weibull",
  theta1.4 = 15,
  theta2.4 = 7.1,  
  pdf5 = "Weibull",
  theta1.5 = 9,
  theta2.5 = 5.3)  

simulated <- RNGMIX(Dataset.name = "simulated",
  rseed = -1,
  n = n,
  Theta = Theta)
  
simulated
}
\keyword{random number generation}

