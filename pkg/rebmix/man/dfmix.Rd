\name{dfmix}
\alias{dfmix}
\title{
Predictive Density Calculation
}
\description{
\code{dfmix} returns (invisibly) the vector containing predictive densities \eqn{f(\bm{x} | c, \bm{w}, \bm{\Theta})}, where \eqn{\bm{x}} stands for a subvector of 
\eqn{\bm{y} = (y_{1}, \ldots, y_{d})^{\top}}.
}
\usage{
dfmix(x = NULL, w = NULL, Theta = NULL, ...)
}
\arguments{
\item{x}{
a vector, a matrix or a data frame containing continuous or discrete vector observations \eqn{\bm{x}}.
}
\item{w}{
a vector or a data frame containing \eqn{c} component weights \eqn{w_{l}} summing to 1.
}
\item{Theta}{
a matrix or a data frame containing  \eqn{c} parametric family types \code{pdfi}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
Component parameters \code{theta1.i} follow the parametric family types. One of \eqn{\mu_{il}} for normal and lognormal distributions and \eqn{\theta_{il}} for Weibull, gamma, binomial, Poisson and Dirac distributions.
Component parameters \code{theta2.i} follow \code{theta1.i}. One of \eqn{\sigma_{il}} for normal and lognormal distributions, \eqn{\beta_{il}} for Weibull and gamma distributions and \eqn{p_{il}} for binomial distribution.
}
\item{\dots}{
currently not used.
}
}
\examples{
## Generate simulated dataset.

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

## Estimate number of components, component weights and component parameters.

simulatedest <- REBMIX(Dataset = simulated$Dataset,
  Preprocessing = "Parzen window",
  cmax = 4,
  Criterion = "BIC",
  Variables = c("continuous", "continuous"),
  pdf = c("normal", "Weibull"),
  K = 8)

## Preprocess and plot finite mixtures.

opar <- plot(simulatedest)

par(opar)

y1f <- demix(x = simulatedest$Dataset[[1]][, 1, drop = FALSE],
  Preprocessing = "Parzen window",
  Variables = "continuous",
  pdf = "normal",
  k = 8)

y1 <- seq(from = min(y1f[, 1]), to = max(y1f[, 1]), length.out = 200)

f1 <- dfmix(x = cbind(y1), w = simulatedest$w[[1]], simulatedest$Theta[[1]][1:3,])

y2f <- demix(x = simulatedest$Dataset[[1]][, 2, drop = FALSE],
  Preprocessing = "Parzen window",
  Variables = "continuous",
  pdf = "Weibull",  
  k = 8)

y2 <- seq(from = min(y2f[, 1]), to = max(y2f[, 1]), length.out = 200)

f2 <- dfmix(x = cbind(y2), w = simulatedest$w[[1]], simulatedest$Theta[[1]][4:6,])

opar <- par(mfrow = c(1, 2))

plot(y1, f1, xlab = bquote(y[1]), ylab = bquote(f(y[1])), type = "l", col = "blue")

points(y1f, pch = 3, col = "red")

plot(y2, f2, xlab = bquote(y[2]), ylab = bquote(f(y[2])), type = "l", col = "blue")

points(y2f, pch = 3, col = "red")

par(opar)
}
\references{
M. Nagode and M. Fajdiga. The rebmix algorithm for the univariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(5):876-892, 2011a. \url{http://dx.doi.org/10.1080/03610920903480890}.\cr\cr
M. Nagode and M. Fajdiga. The rebmix algorithm for the multivariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(11):2022-2034, 2011b. \url{http://dx.doi.org/10.1080/03610921003725788}. 
}
\keyword{distributions}
