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
M. H. Hansen and B. Yu. Model selection and the principle of minimum description
length. Journal of the American Statistical Association, 96(454):746-774, 2001. \url{http://www.jstor.org/stable/2670311}.
}
\keyword{information criterion}
