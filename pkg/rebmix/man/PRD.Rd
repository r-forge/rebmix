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
M. Nagode and M. Fajdiga. The rebmix algorithm for the univariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(5):876-892, 2011a. \url{http://dx.doi.org/10.1080/03610920903480890}.\cr\cr
M. Nagode and M. Fajdiga. The rebmix algorithm for the multivariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(11):2022-2034, 2011b. \url{http://dx.doi.org/10.1080/03610921003725788}.
}
\keyword{information criterion}
