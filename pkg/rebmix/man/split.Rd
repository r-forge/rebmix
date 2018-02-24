\name{split-methods}
\alias{split}
\alias{split-methods}
\alias{split,numeric-method}
\alias{split,list-method}
\title{
Splits Dataset into Train and Test Datasets
}
\description{
Returns (invisibly) the object containing train and test observations \eqn{\bm{y}_{1}, \ldots, \bm{y}_{n}} as well as true class membership \eqn{\bm{\Omega}_{g}} for the test dataset.
}
\usage{
\S4method{split}{numeric}(p = 0.75, Dataset = data.frame(), class = numeric(), ...)
\S4method{split}{list}(p = list(), Dataset = data.frame(), class = numeric(), ...)
## ... and for other signatures
}
\arguments{
\item{p}{
see Methods section below.
}
\item{Dataset}{
a data frame containing dataset \eqn{Y} of length \eqn{n}. For the dataset the corresponding class membership \eqn{\bm{\Omega}_{g}} is known.
The default value is \code{data.frame()}.
}
\item{class}{
a column number in \code{Dataset} containing the class membership information. The default value is \code{numeric()}.
}
\item{\dots}{
further arguments to \code{\link[base]{sample}}.
}
}
\value{
Returns an object of class \code{RCLS.chunk}.
}
\section{Methods}{
\describe{
\item{\code{signature(p = "numeric")}}{a number specifying the fraction of observations for training \eqn{0.0 \leq p \leq 1.0}. The default value is \code{0.75}.}
\item{\code{signature(p = "list")}}{a list composed of column number \code{p$type} in \code{Dataset} containing the type membership information followed by the corresponding train \code{p$train} and test \code{p$test} values.
The default value is \code{list()}.}
}
}
\author{Marko Nagode}
\examples{
data("iris")

# Split dataset into train (75%) and test (25%) subsets.

set.seed(5)

Iris <- split(p = 0.75, Dataset = iris, class = 5)

Iris
}
\keyword{classification}

