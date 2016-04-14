\name{HQC-methods}
\alias{HQC}
\alias{HQC-methods}
\alias{HQC,REBMIX-method}
\alias{HQC,REBMVNORM-method}
\title{
Hannan-Quinn Information Criterion
}
\description{
Returns the Hannan-Quinn information criterion at \code{pos}.
}
\usage{
\S4method{HQC}{REBMIX}(x = NULL, pos = 1, ...)
## ... and for other signatures
}
\arguments{
\item{x}{
see Methods section below.
}
\item{pos}{
a desired row number in \code{x@summary} for which the information criterion is calculated. The default value is \code{1}.
}
\item{\dots}{
currently not used.
}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "REBMIX")}}{an object of class \code{REBMIX}.}
\item{\code{signature(x = "REBMVNORM")}}{an object of class \code{REBMVNORM}.}
}
}
\author{Marko Nagode}
\references{
E. J. Hannan and B. G. Quinn. The determination of the order of an autoregression.
Journal of the Royal Statistical Society. Series B, 41(2):190-195, 1979. \url{http://www.jstor.org/stable/2985032}.
}
\keyword{information criterion}
