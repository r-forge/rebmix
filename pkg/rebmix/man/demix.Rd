\name{demix}
\alias{demix}
\title{
Empirical Density Calculation
}
\description{
\code{demix} returns (invisibly) the data frame containing observations \eqn{\bm{x}_{1}, \ldots, \bm{x}_{n}} and empirical 
densities \eqn{f_{1}, \ldots, f_{n}} for the Parzen window or \emph{k}-nearest neighbour or bin means \eqn{\bar{\bm{x}}_{1}, \ldots, \bar{\bm{x}}_{v}} 
and empirical densities \eqn{f_{1}, \ldots, f_{v}} for the histogram preprocessing. Vectors \eqn{\bm{x}} and \eqn{\bar{\bm{x}}} are subvectors of 
\eqn{\bm{y} = (y_{1}, \ldots, y_{d})^{\top}} and \eqn{\bar{\bm{y}} = (\bar{y}_{1}, \ldots, \bar{y}_{d})^{\top}}.
}
\usage{
demix(x = NULL, Preprocessing = NULL, Variables = NULL, pdf = NULL,
      k = NULL, xmin = NULL, xmax = NULL, ...)
}
\arguments{
\item{x}{
a vector, a matrix or a data frame containing continuous or discrete vector observations \eqn{\bm{x}}.
}
\item{Preprocessing}{
a preprocessing type. One of \code{"histogram"}, \code{"Parzen window"} or \code{"k-nearest neighbour"}.
}
\item{Variables}{
a character vector of length \eqn{\leq d} containing types of variables. One of \code{"continuous"} or \code{"discrete"}.
}
\item{pdf}{
a character vector of length \eqn{d} containing continuous or discrete parametric family types. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
}
\item{k}{
a number of bins \eqn{v} for the histogram and the Parzen window or number of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
\item{xmin}{
a vector of length \eqn{\leq d} containing minimum observations. The default value is \code{NULL}.
}
\item{xmax}{
a vector of length \eqn{\leq d} containing maximum observations. The default value is \code{NULL}.
}
\item{\dots}{
currently not used.
}
}
\examples{
## Generate simulated dataset.

Theta <- rbind(pdf1 = rep("normal", 2),
  theta1.1 = c(10, 20),
  theta2.1 = c(3.0, 2.0),
  pdf1 = rep("normal", 2),
  theta1.1 = c(3, 2),
  theta2.1 = c(20, 10))

simulated <- RNGMIX(Dataset = "simulated",
  rseed = -1,
  n = c(15, 15),
  Theta = Theta)

## Preprocess simulated dataset.

y1y2f <- demix(x = simulated$Dataset[[1]],
  Preprocessing = "histogram",
  Variables = c("continuous", "continuous"),
  pdf = c("normal", "normal"),
  k = 6)
  
y1y2f
}
\references{
M. Nagode and M. Fajdiga. The rebmix algorithm for the univariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(5):876-892, 2011a. \url{http://dx.doi.org/10.1080/03610920903480890}.\cr\cr
M. Nagode and M. Fajdiga. The rebmix algorithm for the multivariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(11):2022-2034, 2011b. \url{http://dx.doi.org/10.1080/03610921003725788}. 
}
\keyword{distributions}
