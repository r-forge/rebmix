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
E. J. Hannan and B. G. Quinn. The determination of the order of an autoregression.
Journal of the Royal Statistical Society. Series B, 41(2):190-195, 1979. URL
http://www.jstor.org/stable/2985032.
}
\keyword{information criterion}
