\name{BIC}
\alias{BIC}
\title{
Bayesian Information Criterion
}
\description{
Returns the Bayesian information criterion at \code{pos} for class REBMIX.
}
\usage{
BIC(object, pos = 1, ...)
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
G. Schwarz. Estimating the dimension of the model. The Annals of Statistics, 6(2):461-464, 1978.
}
\keyword{information criterion}
