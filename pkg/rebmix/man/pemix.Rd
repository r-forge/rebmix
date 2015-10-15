\name{pemix}
\alias{pemix}
\title{
Empirical Distribution Function Calculation
}
\description{
\code{pemix} returns (invisibly) the data frame containing observations \eqn{\bm{x}_{1}, \ldots, \bm{x}_{n}} and empirical 
distribution functions \eqn{F_{1}, \ldots, F_{n}}, where \eqn{\bm{x} = (x_{1}, \ldots, x_{d})^{\top}}.
}
\usage{
pemix(x = NULL, lower.tail = TRUE, log.p = FALSE, ...)
}
\arguments{
\item{x}{
a vector, a matrix or a data frame containing continuous or discrete vector observations \eqn{\bm{x}}.
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

Theta <- list(pdf1 = rep("normal", 2),
  theta1.1 = c(10, 20),
  theta2.1 = c(3, 2),
  pdf2 = rep("normal", 2),
  theta1.2 = c(3, 2),
  theta2.2 = c(20, 10))

simulated <- RNGMIX(Dataset = "simulated",
  rseed = -1,
  n = c(15, 15),
  Theta = Theta)

# Preprocess simulated dataset.

y1y2F <- pemix(x = simulated$Dataset[[1]])
  
y1y2F
}
\keyword{distributions}
