\name{MDL}
\alias{MDL2}
\alias{MDL2.default}
\alias{MDL2.REBMIX}
\alias{MDL5}
\alias{MDL5.default}
\alias{MDL5.REBMIX}
\title{
Minimum Description Length
}
\description{
Returns the minimum desription length at \code{pos} for class REBMIX.
}
\usage{
\method{MDL2}{REBMIX}(object, pos = 1, ...)
\method{MDL5}{REBMIX}(object, pos = 1, ...)

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
Liang Z (1992). Parameter Estimation of Finite Mixture Using the EM-Algorithm and
Information Criteria with Applications to Medical Image Processing. \emph{IEEE: Nuclear
Science}, 39(4), 11-26.
}
\keyword{information criterion}
