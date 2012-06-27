\name{RNGMIX}
\alias{RNGMIX}
\title{
Random Univariate or Multivariate Finite Mixture Generation
}
\description{
Returns univariate or multivariate random datasets for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson or Dirac component densities.
}
\usage{
RNGMIX(Dataset = NULL, rseed = -1, n = NULL, Theta = NULL)
}
\arguments{
\item{Dataset}{
a character vector containing list names of data frames of size \eqn{n \times d} that \emph{d}-dimensional datasets are written in.
}
\item{rseed}{
set the random seed to any negative integer value to initialize the sequence. The first file in \code{Dataset} corresponds to it.
For each next file the random seed is decremented \eqn{r_{\mathrm{seed}} = r_{\mathrm{seed}} - 1}. The default value is \code{-1}.
}
\item{n}{
a vector or a data frame containing number of observations in classes \eqn{n_{l}}, where number of observations \eqn{n = \sum_{l = 1}^{c} n_{l}}.
}
\item{Theta}{
a matrix or a data frame containing \eqn{c} parametric family types \code{pdfi}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
Component parameters \code{theta1.i} follow the parametric family types. One of \eqn{\mu_{il}} for normal and lognormal distributions and \eqn{\theta_{il}} for Weibull, gamma, binomial, Poisson and Dirac distributions.
Component parameters \code{theta2.i} follow \code{theta1.i}. One of \eqn{\sigma_{il}} for normal and lognormal distributions, \eqn{\beta_{il}} for Weibull and gamma distributions and \eqn{p_{il}} for binomial distribution.
}
}
\details{
RNGMIX is based on the "Minimal" random number generator \code{ran1} of Park and Miller with the Bays-Durham shuffle and added safeguards that returns a uniform random deviate between 0.0 and 1.0 
(exclusive of the endpoint values).
}
\value{
\item{Dataset}{
a list of data frames of size \eqn{n \times d} containing \emph{d}-dimensional datasets. Each of the \eqn{d} columns represents one random variable. Number of observations \eqn{n} equals the number of rows
in the datasets.
}
\item{w}{
a data frame containing \eqn{c} component weights \eqn{w_{l}} summing to 1.
}
\item{Theta}{
a data frame containing \eqn{c} parametric family types \code{pdfi}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
Component parameters \code{theta1.i} follow the parametric family types. One of \eqn{\mu_{il}} for normal and lognormal distributions and \eqn{\theta_{il}} for Weibull, gamma, binomial, Poisson and Dirac distributions.
Component parameters \code{theta2.i} follow \code{theta1.i}. One of \eqn{\sigma_{il}} for normal and lognormal distributions, \eqn{\beta_{il}} for Weibull and gamma distributions and \eqn{p_{il}} for binomial distribution.
}
\item{Variables}{
a character vector containing types of variables. One of \code{"continuous"} or \code{"discrete"}.
}
\item{ymin}{
a vector of length \eqn{d} containing minimum observations.
}
\item{ymax}{
a vector of length \eqn{d} containing maximum observations.
}
}
\references{
Press WH, Teukolsky SA, Vetterling VT, Flannery BP (1992). \emph{Numerical Recipes in C: The Art
of Scientific Computing}. Cambridge University Press, Cambridge.
}
\examples{
## Generate and print simulated dataset.

n <- c(75, 100, 125, 150, 175)

Theta <- rbind(pdf1 = rep("normal", 5),
  theta1.1 = c(10, 8.5, 12, 13, 7),
  theta2.1 = c(1, 1, 1, 2, 3),
  pdf2 = rep("normal", 5),
  theta1.2 = c(12, 10, 14, 15, 9),
  theta2.2 = c(1, 1, 1, 2, 3),
  pdf3 = rep("normal", 5),
  theta1.3 = c(10, 8.5, 12, 7, 13),
  theta2.3 = c(1, 1, 1, 2, 3),
  pdf4 = rep("normal", 5),
  theta1.4 = c(12, 10.5, 14, 9, 15),
  theta2.4 = c(1, 1, 1, 2, 3))

simulated <- RNGMIX(Dataset = paste("simulated_", 1:25, sep = ""),
  rseed = -1,
  n = n,
  Theta = Theta)

## Generate and print simulated dataset.

Theta <- rbind(pdf1 = rep("normal", 5),
  theta1.1 = c(10, 8.5, 12, 13, 7),
  theta2.1 = c(1, 1, 1, 2, 3),
  pdf2 = rep("Weibull", 5),
  theta1.2 = c(12, 10, 14, 15, 9),
  theta2.2 = c(2, 4.1, 3.2, 7.1, 5.3),
  pdf3 = rep("binomial", 5),
  theta1.3 = c(10, 8, 12, 7, 13),
  theta2.3 = c(0.1, 0.5, 0.9, 0.4, 0.2),
  pdf4 = rep("Dirac", 5),
  theta1.4 = c(1, 1, 1, 1, 1))

simulated <- RNGMIX(Dataset = paste("simulated_", 1:25, sep = ""),
  rseed = -1,
  n = n,
  Theta = Theta)
  
simulated

## Generate and print simulated dataset.

Theta <- rbind(pdf = rep("W", 5),
  theta1 = c(12, 10, 14, 15, 9),
  theta2 = c(2, 4.1, 3.2, 7.1, 5.3))

simulated <- RNGMIX(Dataset = paste("simulated_", 1:2, sep = ""),
  rseed = -1,
  n = n,
  Theta = Theta)
  
simulated
}
\keyword{random number generation}

