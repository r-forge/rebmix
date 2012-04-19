\name{BIC}
\alias{BIC.REBMIX}
\title{
Bayesian Information Criterion
}
\description{
Returns the Bayesian information criterion at \code{pos} for class REBMIX.
}
\usage{
\method{BIC}{REBMIX}(object, pos = 1, ...)
}
\arguments{
\item{object}{
an object of class \code{REBMIX}.
}
\item{pos}{
a desired row number in \code{summary} for which the information criterion is calculated. The default value is \code{1}.
}
\item{\dots}{
currently not used.
}
}
\references{
Schwarz G (1978). Estimating the Dimension of a Model. \emph{Annals of Statistics}, 6, 461-464.
}
\keyword{information criterion}
