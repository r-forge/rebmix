\name{ICLBIC}
\alias{ICLBIC}
\alias{ICLBIC.default}
\alias{ICLBIC.REBMIX}
\title{
Approximate Integrated Classification Likelihood Criterion
}
\description{
Returns the approximate integrated classification likelihood criterion at \code{pos} for class REBMIX.
}
\usage{
\method{ICLBIC}{REBMIX}(object, pos = 1, ...)
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
C. Biernacki, G. Celeux, and G. Govaert. Assessing a mixture model for clustering with the integrated
classification likelihood. Technical Report 3521, INRIA, Rhône-Alpes, 1998.
}
\keyword{information criterion}
