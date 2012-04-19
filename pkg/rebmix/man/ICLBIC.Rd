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
Biernacki C, Celeux G, Govaert G (1998). Assessing a Mixture Model for Clustering with
the Integrated Classification Likelihood. \emph{Technical Report 3521}, INRIA, Rhone-Alpes.
}
\keyword{information criterion}
