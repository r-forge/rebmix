\name{SSE}
\alias{SSE}
\alias{SSE.default}
\alias{SSE.REBMIX}
\title{
Sum of Squares Error
}
\description{
Returns the sum of squares error at \code{pos} for class REBMIX.
}
\usage{
\method{SSE}{REBMIX}(object, pos = 1, ...)
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
Bishop CM (1995). \emph{Neural Networks for Pattern Recognition}. Clarendon Press, Oxford.
}
\keyword{information criterion}
