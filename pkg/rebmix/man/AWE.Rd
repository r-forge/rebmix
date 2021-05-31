\name{AWE-methods}
\docType{methods}
\alias{AWE}
\alias{AWE-methods}
\alias{AWE,REBMIX-method}
\alias{AWE,REBMVNORM-method}
\title{
Approximate Weight of Evidence Criterion
}
\description{
Returns the approximate weight of evidence criterion at \code{pos}.
}
\usage{
\S4method{AWE}{REBMIX}(x = NULL, pos = 1, ...)
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
J. D. Banfield and A. E. Raftery. Model-based gaussian and non-gaussian clustering. Biometrics, 49 (3):803-821, 1993. \doi{10.2307/2532201}.
}
\keyword{information criterion}
