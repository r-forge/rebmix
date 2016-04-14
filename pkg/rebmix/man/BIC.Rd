\name{BIC-methods}
\alias{BIC}
\alias{BIC-methods}
\alias{BIC,REBMIX-method}
\alias{BIC,REBMVNORM-method}
\title{
Bayesian Information Criterion
}
\description{
Returns the Bayesian information criterion at \code{pos}.
}
\usage{
\S4method{BIC}{REBMIX}(x = NULL, pos = 1, ...)
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
G. Schwarz. Estimating the dimension of the model. The Annals of Statistics, 6(2):461-464, 1978.
}
\keyword{information criterion}
