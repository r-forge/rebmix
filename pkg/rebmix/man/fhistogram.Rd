\name{fhistogram-methods}
\alias{fhistogram}
\alias{fhistogram-methods}
\alias{fhistogram,Histogram-method}
\title{
Fast Histogram Calculation
}
\description{
Returns an object of class \code{Histogram}. The method can be called recursively.
This way more than one dataset can be binned into one histogram. Set \code{shrink} 
to \code{TRUE} only when the method is called for the last time to optimize the size of the object.
The method is memory consuming.
}
\usage{
\S4method{fhistogram}{Histogram}(x = NULL, Dataset = data.frame(),
           K = numeric(), ymin = numeric(), ymax = numeric(), 
           shrink = FALSE, ...)
## ... and for other signatures
}
\arguments{
\item{x}{
an object of class \code{Histogram}.
}
\item{Dataset}{
a data frame of size \eqn{n \times d} containing \emph{d}-dimensional dataset. Each of the \eqn{d} columns
represents one random variable. Number of observations \eqn{n} equals the number of rows in the dataset.
}
\item{K}{
an integer or a vector of length \eqn{d} containing numbers of bins \eqn{v}.
}
\item{ymin}{
a vector of length \eqn{d} containing minimum observations.
}
\item{ymax}{
a vector of length \eqn{d} containing maximum observations.
}
\item{shrink}{
logical. If \code{TRUE} the output is shrank to its optimal size. The default value is \code{FALSE}.
}
\item{\dots}{
currently not used.
}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "Histogram")}}{an object of class \code{Histogram}.}
}
}
\author{Marko Nagode}
\examples{

# Create three datasets.

set.seed(1)

n <- 15

Dataset1 <- as.data.frame(cbind(rnorm(n, 157, 8), rnorm(n, 71, 10)))
Dataset2 <- as.data.frame(cbind(rnorm(n, 244, 14), rnorm(n, 61, 29)))
Dataset3 <- as.data.frame(cbind(rnorm(n, 198, 8), rnorm(n, 252, 13)))

apply(Dataset1, 2, range)
apply(Dataset2, 2, range)
apply(Dataset3, 2, range)

# Bin the first dataset.

hist <- fhistogram(Dataset = Dataset1, K = c(4, 5), ymin = c(100.0, 0.0), ymax = c(300.0, 300.0))

# Bin the second dataset.

hist <- fhistogram(x = hist, Dataset = Dataset2)

# Bin the third dataset and shrink the hist object.

hist <- fhistogram(x = hist, Dataset = Dataset3, shrink = TRUE)

hist
}
\keyword{auxiliary}
