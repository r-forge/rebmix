\name{Histogram-class}
\docType{class}
\alias{Histogram-class}
\title{Class \code{"Histogram"}}
\description{Object of class \code{Histogram}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Histogram", ...)}. Accessor methods for the slots are \code{a.Y(x = NULL)} and \code{a.h(x = NULL)}.
}
\section{Slots}{
\describe{
\item{\code{Y}:}{
a data frame of size \eqn{v \times (d + 1)} containing \emph{d}-dimensional histogram.
Each of the first \eqn{d} columns represents one random variable and contains bin means 
\eqn{\bar{\bm{y}}_{1}, \ldots, \bar{\bm{y}}_{v}}. Column \eqn{d + 1} contains frequencies \eqn{k_{1}, \ldots, k_{v}}.
}
\item{\code{h}:}{
a vector of length \eqn{d} containing bin widths.
}
}
}
\author{Marko Nagode}
\examples{
}
\keyword{classes}
