\name{CLC-methods}
\alias{CLC}
\alias{CLC-methods}
\alias{CLC,REBMIX-method}
\alias{CLC,REBMVNORM-method}
\title{
Classification Likelihood Criterion
}
\description{
Returns the classification likelihood criterion at \code{pos}.
}
\usage{
\S4method{CLC}{REBMIX}(x = NULL, pos = 1, ...)
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
C. Biernacki and G. Govaert. Using the classification likelihood to choose the number of clusters. In
E. J. Wegman and S. P. Azen, editors, Computing Science and Statistics, 1997.
}
\keyword{information criterion}
