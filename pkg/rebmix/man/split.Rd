\name{split-methods}
\alias{split}
\alias{split-methods}
\alias{split,numeric-method}
\alias{split,list-method}
\title{
Splits Dataset into Train and Test Datasets
}
\description{
Returns (invisibly) the \code{s}, \code{ntrain}, \code{train}, \code{ntest}, \code{test} and \code{Zt} output for \code{Dataset}.
}
\usage{
\S4method{split}{numeric}(p = 0.75, Dataset = data.frame(), class = numeric(), ...)
\S4method{split}{list}(p = list(), Dataset = data.frame(), class = numeric(), ...)
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
further arguments to \code{\link[base]{sample.int}}.
}
}
\value{
\item{s}{
finite set of size \eqn{s} of classes \eqn{\bm{\Omega} = \{\bm{\Omega}_{g}; \ g = 1, \ldots, s\}}.
}
\item{ntrain}{
a vector of length \eqn{s} containing numbers of observations in train datasets \eqn{Y_{\mathrm{train}g}}.
}
\item{train}{
a list of data frames containing train datasets \eqn{Y_{\mathrm{train}g}} of length \eqn{n_{\mathrm{train}g}}.
}
\item{ntest}{
number of observations in test dataset \eqn{Y_{\mathrm{test}}}.
}
\item{test}{
a data frame containing test dataset \eqn{Y_{\mathrm{test}}} of length \eqn{n_{\mathrm{test}}}.
}
\item{Zt}{
a factor of true class membership \eqn{\bm{\Omega}_{g}} for the test dataset.
}
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
}
\keyword{classification}

