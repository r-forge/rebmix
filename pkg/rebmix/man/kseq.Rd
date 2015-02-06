\name{kseq}
\alias{kseq}
\title{
Sequence of Bins or Nearest Neighbours Generation
}
\description{
\code{kseq} returns (invisibly) a vector containing numbers of bins \eqn{v} for the histogram and the Parzen window or numbers of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
\usage{
kseq(from = NULL, to = NULL, f = 0.05, ...)
}
\arguments{
\item{from}{
starting value of the sequence. The default value is \code{NULL}.
}
\item{to}{
end value of the sequence. The default value is \code{NULL}.
}
\item{f}{
number specifying the fraction by which the bins or nearest neighbours should be separated \eqn{0.0 < f < 1.0}. The default value is \code{0.05}.
}
\item{\dots}{
currently not used.
}
}
\examples{
## Generate numbers of bins.

n <- 10000

Sturges <- as.integer(1 + log2(n)) ## Minimum v follows the Sturges rule.
Log10 <- as.integer(10 * log10(n)) ## Maximum v follows the Log10 rule.
RootN <- as.integer(2 * n^0.5) ## Maximum v follows the RootN rule.

K <- kseq(from = Sturges, to = Log10, f = 0.05)

K

K <- kseq(from = Sturges, to = RootN, f = 0.03)

K
}
\keyword{parameter estimation}

