\name{logL}
\alias{logL}
\title{
Log Likelihood
}
\description{
Returns the log likelihood at \code{pos} for class REBMIX.
}
\usage{
logL(object, pos = 1, ...)
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
G. McLachlan and D. Peel. Finite Mixture Models. John Wiley & Sons, New York, 2000.
}
\keyword{information criterion}
