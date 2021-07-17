\name{chunk-methods}
\alias{chunk}
\alias{chunk-methods}
\alias{chunk,RCLS.chunk-method}
\alias{show,RCLS.chunk-method}
\title{
Extracts Chunk from Train and Test Datasets
}
\description{
Returns (invisibly) the object containing train and test observations \eqn{\bm{x}_{1}, \ldots, \bm{x}_{n}} as well as true class membership \eqn{\bm{\Omega}_{g}} for the test dataset. Vectors \eqn{\bm{x}} are subvectors of
\eqn{\bm{y} = (y_{1}, \ldots, y_{d})^{\top}}.
}
\usage{
\S4method{chunk}{RCLS.chunk}(x = NULL, variables = expression(1:d))
## ... and for other signatures
}
\arguments{
\item{x}{
see Methods section below.
}
\item{variables}{
a vector containing indices of variables in subvectors \eqn{\bm{x}}. The default value is \code{1:d}.
}
}
\value{
Returns an object of class \code{RCLS.chunk}.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "RCLS.chunk")}}{an object of class \code{RCLS.chunk}.}
}
}
\author{Marko Nagode}
\examples{
data(iris)

# Split dataset into train (75%) and test (25%) subsets.

set.seed(5)

Iris <- split(p = 0.75, Dataset = iris, class = 5)

# Extract chunk from train and test datasets.

Iris14 <- chunk(x = Iris, variables = c(1,4))

Iris14
}
\keyword{classification}

