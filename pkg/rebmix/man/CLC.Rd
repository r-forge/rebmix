\name{CLC}
\alias{CLC}
\alias{CLC.default}
\alias{CLC.REBMIX}
\title{
Classification Likelihood Criterion
}
\description{
Returns the classification likelihood criterion at \code{pos} for class REBMIX.
}
\usage{
\method{CLC}{REBMIX}(object, pos = 1, ...)
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
C. Biernacki and G. Govaert. Using the classification likelihood to choose the number of clusters. In
E. J. Wegman and S. P. Azen, editors, Computing Science and Statistics, 1997.
}
\keyword{information criterion}
