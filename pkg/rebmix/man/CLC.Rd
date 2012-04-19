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
Biernacki C, Govaert G (1997). Using the Classification Likelihood to Choose the Number
of Clusters. \emph{Computing Science and Statistics}, 29, 451-457.
}
\keyword{information criterion}
