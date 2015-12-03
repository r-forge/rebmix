\name{logL}
\docType{methods}
\alias{logL}
\alias{logL-methods}
\alias{logL,REBMIX-method}
\alias{logL,REBMVNORM-method}
\title{
Log Likelihood
}
\description{
Returns the log likelihood at \code{pos}.
}
\usage{
\S4method{logL}{REBMIX}(x = NULL, pos = 1, ...)
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
G. McLachlan and D. Peel. Finite Mixture Models. John Wiley & Sons, New York, 2000.
}
\keyword{information criterion}
