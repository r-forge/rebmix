\name{ICL-methods}
\alias{ICL}
\alias{ICL-methods}
\alias{ICL,REBMIX-method}
\alias{ICL,REBMVNORM-method}
\title{
Integrated Classification Likelihood Criterion
}
\description{
Returns the integrated classification likelihood criterion at \code{pos}.
}
\usage{
\S4method{ICL}{REBMIX}(x = NULL, pos = 1, ...)
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
C. Biernacki, G. Celeux, and G. Govaert. Assessing a mixture model for clustering with the integrated
classification likelihood. Technical Report 3521, INRIA, Rhone-Alpes, 1998.
}
\keyword{information criterion}
