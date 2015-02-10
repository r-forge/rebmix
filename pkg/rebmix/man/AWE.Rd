\name{AWE}
\alias{AWE}
\title{
Approximate Weight of Evidence Criterion
}
\description{
Returns the approximate weight of evidence criterion at \code{pos} for class REBMIX.
}
\usage{
AWE(object, pos = 1, ...)
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
J. D. Banfield and A. E. Raftery. Model-based gaussian and non-gaussian clustering. Biometrics, 49
(3):803-821, 1993. \url{http://dx.doi.org/10.2307/2532201}.
}
\keyword{information criterion}
