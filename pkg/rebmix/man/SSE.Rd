\name{SSE-methods}
\alias{SSE}
\alias{SSE-methods}
\alias{SSE,REBMIX-method}
\alias{SSE,REBMVNORM-method}
\title{
Sum of Squares Error
}
\description{
Returns the sum of squares error at \code{pos}.
}
\usage{
\S4method{SSE}{REBMIX}(x = NULL, pos = 1, ...)
## ... and for other signatures
}
\arguments{
\item{x}{
see Methods section below.
}
\item{pos}{
a desired row number in \code{summary} for which the information criterion is calculated. The default value is \code{1}.
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
C. M. Bishop. Neural Networks for Pattern Recognition. Clarendon Press, Oxford, 1995.
}
\keyword{information criterion}
