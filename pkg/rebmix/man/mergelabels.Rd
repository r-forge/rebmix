\name{mergelabels-methods}
\alias{mergelabels}
\alias{mergelabels-methods}
\alias{mergelabels,list-method}
\title{
Merge Labels Based on Composed Adjacency Matrix
}
\description{
Returns the list with the normalised adjacency matrix \code{L} of size \eqn{c \times c}. The normalised adjacency matrix \eqn{L = \sum_{l = 1}^{n} w_{l} D_{l}^{-1/2} A_{l} D_{l}^{-1/2}} depends on the adjacency matrices \eqn{A_{l}}
and the degree matrices \eqn{D_{l}(i,j) = \sum_{j} A_{l}(i,j)}. The \eqn{A_{l}} matrices may have some \code{NA} rows and columns. The list also contains the vector of integers \code{cluster} of length \eqn{k}, which indicates the cluster to which each image moment \eqn{\bm{M}} is assigned.
}
\usage{
\S4method{mergelabels}{list}(A = list(), w = numeric(), k = 2, ...)
## ... and for other signatures
}
\arguments{
\item{A}{
a list of length \eqn{n} of adjacency matrices \eqn{A_{l}} of size \eqn{c_{\mathrm{max}} \times c_{\mathrm{max}}}, where \eqn{c_{\mathrm{max}} \geq c}.
}
\item{w}{
vector of length \eqn{n} containing weights \eqn{w_{l}} summing to 1.
}
\item{k}{
number of clusters \eqn{k} for \code{\link[stats]{kmeans}}. The default value is 2.
}
\item{\dots}{
further arguments to \code{\link[stats]{kmeans}}.
}
}
\section{Methods}{
\describe{
\item{\code{signature(A = "list")}}{a list.}
}
}
\author{Marko Nagode, Branislav Panic}
\references{
A. Ng, M. Jordan and Y. Weiss. On spectral clustering: Analysis and an algorithm. Advances in Neural Information Processing Systems 14 (NIPS 2001).
}
\examples{
i <- 1
}
\keyword{auxiliary}
