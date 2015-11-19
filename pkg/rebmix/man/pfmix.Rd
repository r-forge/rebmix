\name{pfmix}
\alias{pfmix}
\title{
Predictive Distribution Function Calculation
}
\description{
\code{pfmix} returns (invisibly) the vector containing predictive distribution functions \eqn{F(\bm{x} | c, \bm{w}, \bm{\Theta})}. Vector \eqn{\bm{x}} ia a subvector of \eqn{\bm{y} = (y_{1}, \ldots, y_{d})^{\top}}.
}
\usage{
pfmix(x = NULL, w = NULL, Theta = NULL, lower.tail = TRUE, log.p = FALSE, ...)
}
\arguments{
\item{x}{
a vector, a matrix or a data frame containing continuous or discrete vector observations \eqn{\bm{x}}.
}
\item{w}{
a vector containing \eqn{c} component weights \eqn{w_{l}} summing to 1.
}
\item{Theta}{
a list containing \eqn{c} parametric family types \code{pdfi}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
Component parameters \code{theta1.i} follow the parametric family types. One of \eqn{\mu_{il}} for normal and lognormal distributions and \eqn{\theta_{il}} for Weibull, gamma, binomial, Poisson and Dirac distributions.
Component parameters \code{theta2.i} follow \code{theta1.i}. One of \eqn{\sigma_{il}} for normal and lognormal distributions, \eqn{\beta_{il}} for Weibull and gamma distributions and \eqn{p_{il}} for binomial distribution.
}
\item{lower.tail}{
logical. If \code{TRUE} (default), probabilities are \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.
}
\item{log.p}{
logical. if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
}
\item{\dots}{
currently not used.
}
}
\examples{
# Generate simulated dataset.

Theta <- list(pdf1 = c("normal", "Weibull"),
  theta1.1 = c(10, 3),
  theta2.1 = c(3, 2),
  pdf2 = c("normal", "Weibull"),
  theta1.2 = c(20, 2),
  theta2.2 = c(2, 10))

simulated <- RNGMIX(Dataset = "simulated",
  rseed = -1,
  n = c(15, 25),
  Theta = Theta)

# Estimate number of components, component weights and component parameters.

simulatedest <- REBMIX(Dataset = simulated$Dataset,
  Preprocessing = "Parzen window",
  cmax = 4,
  Criterion = "BIC",
  pdf = c("normal", "Weibull"),
  K = 8)

# Predictive distribution calculation.

i <- 1

x <- simulatedest$Dataset[[1]][, i]

F <- pfmix(x = x,
  w = simulatedest$w[[1]],
  Theta = lapply(simulatedest$Theta[[1]], "[", i))
  
plot(x, F)

i <- 2

x <- simulatedest$Dataset[[1]][, i]

F <- pfmix(x = x,
  w = simulatedest$w[[1]],
  Theta = lapply(simulatedest$Theta[[1]], "[", i))
  
plot(x, F)
}
\references{
M. Nagode and M. Fajdiga. The rebmix algorithm for the univariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(5):876-892, 2011a. \url{http://dx.doi.org/10.1080/03610920903480890}.\cr\cr
M. Nagode and M. Fajdiga. The rebmix algorithm for the multivariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(11):2022-2034, 2011b. \url{http://dx.doi.org/10.1080/03610921003725788}.
}
\keyword{distributions}
