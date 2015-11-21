\name{HQC}
\alias{HQC}
\title{
Hannan-Quinn Information Criterion
}
\description{
Returns the Hannan-Quinn information criterion at \code{pos} for classes \code{REBMIX} and \code{REBMVNORM}.
}
\usage{
HQC(x, pos = 1, ...)
}
\arguments{
\item{x}{
an object of class \code{REBMIX}.
}
\item{pos}{
a desired row number in \code{summary} for which the information criterion is calculated. The default value is \code{1}.
}
\item{\dots}{
currently not used.
}
}
\author{Marko Nagode}
\references{
E. J. Hannan and B. G. Quinn. The determination of the order of an autoregression.
Journal of the Royal Statistical Society. Series B, 41(2):190-195, 1979. \url{http://www.jstor.org/stable/2985032}.
}
\keyword{information criterion}
