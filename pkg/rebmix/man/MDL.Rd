\name{MDL-methods}
\alias{MDL2}
\alias{MDL5}
\alias{MDL2-methods}
\alias{MDL2,REBMIX-method}
\alias{MDL2,REBMVNORM-method}
\alias{MDL5-methods}
\alias{MDL5,REBMIX-method}
\alias{MDL5,REBMVNORM-method}
\title{
Minimum Description Length
}
\description{
Returns the minimum desription length at \code{pos}.
}
\usage{
\S4method{MDL2}{REBMIX}(x = NULL, pos = 1, ...)
\S4method{MDL5}{REBMIX}(x = NULL, pos = 1, ...)
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
M. H. Hansen and B. Yu. Model selection and the principle of minimum description
length. Journal of the American Statistical Association, 96(454):746-774, 2001. https://www.jstor.org/stable/2670311.
}
\keyword{information criterion}
