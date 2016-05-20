\name{REBMIX.combine-class}
\docType{class}
\alias{REBMIX.combine-class}
\alias{REBMVNORM.combine-class}
\title{Class \code{"REBMIX.combine"}}
\description{Object of class \code{REBMIX.combine}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("REBMIX.combine", ...)}. 
}
\section{Slots}{
\describe{
\item{\code{x}:}{
an object of class \code{REBMIX}.
}
\item{\code{pos}:}{
a desired row number in \code{x@summary} to be combined. The default value is \code{1}.
}
\item{\code{from}:}{
a vector of length \eqn{c} containing clusters merged to \code{to} clusters.
}
\item{\code{to}:}{
a vector of length \eqn{c} containing clusters originating from \code{from} clusters.
}
\item{\code{EN}:}{
a vector of length \eqn{c} containing entropies for combined clusters.
}
\item{\code{ED}:}{
a vector of length \eqn{c} containing decrease of entropies for combined clusters.
}
}
}
\author{Marko Nagode}
\references{
J. P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo. Combining mixture components for clustering.
Journal of Computational and Graphical Statistics, 19(2):332-353, 2010. \url{http://dx.doi.org/10.1198/jcgs.2010.08111}.
}
\keyword{classes}
