\name{BIC}
\alias{BIC}
\title{
Bayesian Information Criterion
}
\description{
Returns the Bayesian information criterion at \code{pos} for classes \code{REBMIX} and \code{REBMVNORM}.
}
\usage{
BIC(x, pos = 1, ...)
}
\arguments{
\item{x}{
an object of class \code{REBMIX}.
}
\item{pos}{
a desired row number in \code{summary} for which the information criterion is calculated. The default value is \code{1}.
}
\item{\dots}{
currently not used.
}
}
\author{Marko Nagode}
\references{
G. Schwarz. Estimating the dimension of the model. The Annals of Statistics, 6(2):461-464, 1978.
}
\keyword{information criterion}
