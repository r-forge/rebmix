\name{mergelabels-methods}
\alias{mergelabels}
\alias{mergelabels-methods}
\alias{mergelabels,list-method}
\title{
Merge Labels Based on Probability Adjacency Matrix
}
\description{
Returns the list with the normalised adjacency matrix \code{L} of size \eqn{c \times c}. The normalised adjacency matrix
\eqn{L = D^{-1/2} P D^{-1/2}} depends on the probability adjacency matrix \eqn{P(i,j) = \sum_{l = 1}^{n} p_{l} A_{l}(i,j)}, where \eqn{p_{l} = w_{l} / \sum_{i = 1}^{c}\sum_{j = i + 1}^{c} A_{l}(i,j)}
and the degree matrix \eqn{D(i,i) = \sum_{j = 1}^{c} P(i,j)}. The \eqn{A_{l}} matrices may contain some \code{NA} rows and columns, which are eliminated by the method. 
The list also contains the vector of integers \code{cluster} of length \eqn{k}, which indicates the cluster to which each label is assigned.
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
Zp <- array(0, dim = c(10, 10, 2))

Zp[ , ,1][10, 1:4] <- 1
Zp[ , ,1][1:4, 10] <- 2

Zp[ , ,2][9, 1:5] <- 3
Zp[ , ,2][1:6, 9] <- 4

labelmoments <- labelmoments(Zp, cmax = 4, Sigma = 1.0)

labelmoments

set.seed(3)

mergelabels <- mergelabels(list(labelmoments$A), w = 1.0, k = 2, nstart = 5)

Zp

mergelabels
}
\keyword{auxiliary}
