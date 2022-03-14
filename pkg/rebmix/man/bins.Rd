\name{bins-methods}
\alias{bins}
\alias{bins-methods}
\alias{bins,list-method}
\title{
Binning of Data
}
\description{
Returns the list of data frames containing bin means \eqn{\bar{\bm{y}}_{1}, \ldots, \bar{\bm{y}}_{v}} and frequencies \eqn{k_{1}, \ldots, k_{v}} for the histogram preprocessing.
}
\usage{
\S4method{bins}{list}(Dataset = list(), K = matrix(),
     ymin = numeric(), ymax = numeric(), ...)
## ... and for other signatures
}
\arguments{
\item{Dataset}{
a list of length \eqn{n_{\mathrm{D}}} of data frames of size \eqn{n \times d} containing \emph{d}-dimensional datasets. Each of the \eqn{d} columns
represents one random variable. Numbers of observations \eqn{n} equal the number of rows in the datasets.
}
\item{K}{
a matrix of size \eqn{n_{\mathrm{D}} \times d} containing numbers of bins \eqn{v_{1}, \ldots, v_{d}} for the histogram. 
If, e.g., \code{K = matrix(c(10, 15, 18, 5, 7, 9), byrow = TRUE, ncol = 3)} than \eqn{d = 3} and the list \code{Dataset} contains \eqn{n_{\mathrm{D}} = 2} data frames. Hence, different numbers of bins can be assigned to \eqn{y_{1}, \ldots, y_{d}}.
The default value is \code{matrix()}.
}
\item{ymin}{
a vector of length \eqn{d} containing minimum observations. The default value is \code{numeric()}.
}
\item{ymax}{
a vector of length \eqn{d} containing maximum observations. The default value is \code{numeric()}.
}
\item{\dots}{
currently not used.
}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "list")}}{a list of data frames.}
}
}
\author{Branislav Panic, Marko Nagode}
\references{
M. Nagode. Finite mixture modeling via REBMIX.
Journal of Algorithms and Optimization, 3(2):14-28, 2015. \url{https://repozitorij.uni-lj.si/Dokument.php?id=127674&lang=eng}.
}
\examples{
# Generate multivariate normal datasets.

n <- c(7, 10)

Theta <- new("RNGMVNORM.Theta", c = 2, d = 2)

a.theta1(Theta, 1) <- c(8, 6)
a.theta1(Theta, 2) <- c(6, 8)
a.theta2(Theta, 1) <- c(8, 2, 2, 4)
a.theta2(Theta, 2) <- c(2, 1, 1, 4)

sim2d <- RNGMIX(model = "RNGMVNORM", 
  Dataset.name = paste("sim2d_", 1:2, sep = ""),
  rseed = -1,
  n = n,
  Theta = a.Theta(Theta))

# Calculate optimal numbers of bins.

opt.k <- optbins(Dataset = sim2d@Dataset,
  Rule = "Knuth equal",
  kmin = 1, 
  kmax = 20)

opt.k

Y <- bins(Dataset = sim2d@Dataset, K = opt.k)

Y

opt.k <- optbins(Dataset = sim2d@Dataset,
  Rule = "Knuth unequal",
  kmin = 1, 
  kmax = 20)

opt.k

Y <- bins(Dataset = sim2d@Dataset, K = opt.k)

Y
}
\keyword{auxiliary}
