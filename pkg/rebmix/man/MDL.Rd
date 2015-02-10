\name{MDL}
\alias{MDL2}
\alias{MDL5}
\title{
Minimum Description Length
}
\description{
Returns the minimum desription length at \code{pos} for class REBMIX.
}
\usage{
MDL2(object, pos = 1, ...)
MDL5(object, pos = 1, ...)
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
