\name{PRD}
\alias{PRD}
\alias{PRD.default}
\alias{PRD.REBMIX}
\title{
Total of Positive Relative Deviations
}
\description{
Returns the total of positive relative deviations \code{D} at \code{pos} for class REBMIX.
}
\usage{
\method{PRD}{REBMIX}(object, pos = 1, ...)
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
Nagode M, Fajdiga M (2011). The REBMIX Algorithm for the Univariate Finite Mixture
Estimation. \emph{Communications in Statistics - Theory and Methods}, 40(5), 876-892.\cr\cr
Nagode M, Fajdiga M (2011). The REBMIX Algorithm for the Multivariate Finite Mixture
Estimation. \emph{Communications in Statistics - Theory and Methods}, 40(11), 2022-2034. 
}
\keyword{information criterion}
