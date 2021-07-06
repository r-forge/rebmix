\name{pemix-methods}
\alias{pemix}
\alias{pemix-methods}
\alias{pemix,REBMIX-method}
\alias{pemix,REBMVNORM-method}
\title{
Empirical Distribution Function Calculation
}
\description{
Returns the data frame containing observations \eqn{\bm{x}_{1}, \ldots, \bm{x}_{n}} and empirical
distribution functions \eqn{F_{1}, \ldots, F_{n}}. Vectors \eqn{\bm{x}} are subvectors of
\eqn{\bm{y} = (y_{1}, \ldots, y_{d})^{\top}}.
}
\usage{
\S4method{pemix}{REBMIX}(x = NULL, pos = 1, variables = expression(1:d),
      lower.tail = TRUE, log.p = FALSE, ...)
## ... and for other signatures
}
\arguments{
\item{x}{
see Methods section below.
}
\item{pos}{
a desired row number in \code{x@summary} for which the empirical distribution functions are calculated. The default value is \code{1}.
}
\item{variables}{
a vector containing indices of variables in subvectors \eqn{\bm{x}}. The default value is \code{1:d}.
}
\item{lower.tail}{
logical. If \code{TRUE}, probabilities are \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}. The default value is \code{TRUE}.
}
\item{log.p}{
logical. if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}. The default value is \code{FALSE}.
}
\item{\dots}{
currently not used.
}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "REBMIX")}}{an object of class \code{REBMIX}.}
\item{\code{signature(x = "REBMVNORM")}}{an object of class \code{REBMVNORM}.}
}
}
\author{Marko Nagode}
\examples{
# Generate simulated dataset.

n <- c(15, 15)

Theta <- new("RNGMIX.Theta", c = 2, pdf = rep("normal", 3))

a.theta1(Theta, 1) <- c(10, 20, 30)
a.theta1(Theta, 2) <- c(3, 4, 5)
a.theta2(Theta, 1) <- c(3, 2, 1)
a.theta2(Theta, 2) <- c(15, 10, 5)

simulated <- RNGMIX(Dataset.name = paste("simulated_", 1:4, sep = ""),
  rseed = -1,
  n = n,
  Theta = a.Theta(Theta))

# Create object of class EM.Control.

EM <- new("EM.Control", strategy = "exhaustive", variant = "ECM",
  acceleration = "fixed", acceleration.multiplier = 1.0, tolerance = 1.0E-4,
  maximum.iterations = 1000)

# Estimate number of components, component weights and component parameters.

simulatedest <- REBMIX(Dataset = a.Dataset(simulated),
  Preprocessing = "kernel density estimation",
  cmax = 4,
  pdf = c("n", "n", "n"),
  EMControl = EM)

# Preprocess simulated dataset.

f <- pemix(simulatedest, pos = 3, variables = c(1))

f
}
\references{
M. Nagode and M. Fajdiga. The rebmix algorithm for the univariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(5):876-892, 2011a. \doi{10.1080/03610920903480890}.\cr\cr
M. Nagode and M. Fajdiga. The rebmix algorithm for the multivariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(11):2022-2034, 2011b. \doi{10.1080/03610921003725788}.\cr\cr
M. Nagode. Finite mixture modeling via REBMIX.
Journal of Algorithms and Optimization, 3(2):14-28, 2015. \doi{10.5963/JAO0302001}.
}
\keyword{distributions}
