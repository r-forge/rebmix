\name{demix-methods}
\alias{demix}
\alias{demix-methods}
\alias{demix,REBMIX-method}
\alias{demix,REBMVNORM-method}
\title{
Empirical Density Calculation
}
\description{
Returns the data frame containing observations \eqn{\bm{x}_{1}, \ldots, \bm{x}_{n}} and empirical 
densities \eqn{f_{1}, \ldots, f_{n}} for the kernel density estimation or \emph{k}-nearest neighbour or bin means \eqn{\bar{\bm{x}}_{1}, \ldots, \bar{\bm{x}}_{v}} 
and empirical densities \eqn{f_{1}, \ldots, f_{v}} for the histogram preprocessing. Vectors \eqn{\bm{x}} and \eqn{\bar{\bm{x}}} are subvectors of 
\eqn{\bm{y} = (y_{1}, \ldots, y_{d})^{\top}} and \eqn{\bar{\bm{y}} = (\bar{y}_{1}, \ldots, \bar{y}_{d})^{\top}}.
}
\usage{
\S4method{demix}{REBMIX}(x = NULL, pos = 1, variables = expression(1:d), ...)
## ... and for other signatures
}
\arguments{
\item{x}{
see Methods section below.
}
\item{pos}{
a desired row number in \code{x@summary} for which the empirical densities are calculated. The default value is \code{1}.
}
\item{variables}{
a vector containing indices of variables in subvectors \eqn{\bm{x}} or \eqn{\bar{\bm{x}}}. The default value is \code{1:d}.
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

Theta <- list(pdf1 = rep("normal", 3),
  theta1.1 = c(10, 20, 30),
  theta2.1 = c(3, 2, 1),
  pdf2 = rep("normal", 3),
  theta1.2 = c(3, 4, 5),
  theta2.2 = c(15, 10, 5))

simulated <- RNGMIX(Dataset.name = paste("simulated_", 1:4, sep = ""),
  rseed = -1,
  n = n,
  Theta = Theta)

# Number of classes or nearest neighbours to be processed.

K <- c(as.integer(1 + log2(sum(n))), # Minimum v follows Sturges rule.
  as.integer(10 * log10(sum(n)))) # Maximum v follows log10 rule.

# Estimate number of components, component weights and component parameters. 

simulatedest <- REBMIX(model = "REBMVNORM",
  Dataset = simulated@Dataset,
  Preprocessing = "h", 
  cmax = 4,
  Criterion = "BIC", 
  pdf = c("n", "n", "n"), 
  K = K[1]:K[2])

# Preprocess simulated dataset.

f <- demix(simulatedest, pos = 3, variables = c(1, 3))

f

# Plot finite mixture.

opar <- plot(simulatedest, pos = 3, nrow = 3, ncol = 1)

par(usr = opar[[2]]$usr, mfg = c(2, 1))

text(x = f[, 1], y = f[, 2], labels = format(f[, 3], digits = 3), cex = 0.8, pos = 1)
}
\references{
M. Nagode and M. Fajdiga. The rebmix algorithm for the univariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(5):876-892, 2011a. \url{http://dx.doi.org/10.1080/03610920903480890}.\cr\cr
M. Nagode and M. Fajdiga. The rebmix algorithm for the multivariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(11):2022-2034, 2011b. \url{http://dx.doi.org/10.1080/03610921003725788}.\cr\cr
M. Nagode. Finite mixture modeling via REBMIX.
Journal of Algorithms and Optimization, 3(2):14-28, 2015. \url{http://dx.doi.org/10.5963/JAO0302001}.
}
\keyword{distributions}
