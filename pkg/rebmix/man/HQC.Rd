\name{HQC}
\alias{HQC}
\alias{HQC.default}
\alias{HQC.REBMIX}
\title{
Hannan-Quinn Information Criterion
}
\description{
Returns the Hannan-Quinn information criterion at \code{pos} for class REBMIX.
}
\usage{
\method{HQC}{REBMIX}(object, pos = 1, ...)
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
Hannan EJ, Quinn BG (1979). The Determination of the Order of an Autoregression.
\emph{Journal of the Royal Statistical Society B}, 41, 190-195.
}
\keyword{information criterion}
