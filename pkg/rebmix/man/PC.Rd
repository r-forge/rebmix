\name{PC}
\alias{PC}
\title{
Partition Coefficient
}
\description{
Returns the partition coefficient of Bezdek at \code{pos} for class REBMIX.
}
\usage{
PC(object, pos = 1, ...)
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
