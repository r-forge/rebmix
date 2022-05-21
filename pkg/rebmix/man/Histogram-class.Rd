\name{Histogram-class}
\docType{class}
\alias{Histogram-class}
\title{Class \code{"Histogram"}}
\description{Object of class \code{Histogram}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Histogram", ...)}. Accessor methods for the slots are \code{a.Y(x = NULL)}, 
\code{a.K(x = NULL)}, \code{a.ymin(x = NULL)}, \code{a.ymax(x = NULL)}, \code{a.y0(x = NULL)}, \code{a.h(x = NULL)}, \code{a.n(x = NULL)} and \code{a.ns(x = NULL)}.
}
\section{Slots}{
\describe{
\item{\code{Y}:}{
a data frame of size \eqn{v \times (d + 1)} containing \emph{d}-dimensional histogram.
Each of the first \eqn{d} columns represents one random variable and contains bin means 
\eqn{\bar{\bm{y}}_{1}, \ldots, \bar{\bm{y}}_{v}}. Column \eqn{d + 1} contains frequencies \eqn{k_{1}, \ldots, k_{v}}.
}
\item{\code{K}:}{
an integer or a vector of length \eqn{d} containing numbers of bins \eqn{v}.
}
\item{\code{ymin}:}{
a vector of length \eqn{d} containing minimum observations.
}
\item{\code{ymax}:}{
a vector of length \eqn{d} containing maximum observations.
}
\item{\code{y0}:}{
a vector of length \eqn{d} containing origins.
}
\item{\code{h}:}{
a vector of length \eqn{d} containing bin widths.
}
\item{\code{n}:}{
an integer containing total number \eqn{n} of observations.
}
\item{\code{ns}:}{
an integer containing number \eqn{n_{\mathrm{s}}} of samples.
}
}
}
\author{Marko Nagode}
\examples{
Y <- as.data.frame(matrix(0.0, nrow = 8, ncol = 3))

hist <- new("Histogram", Y = Y, K = c(4, 2), ymin = c(2, 1), ymax = c(10, 8))

a.Y(hist)
a.K(hist)
a.ymin(hist)
a.ymax(hist)
a.y0(hist)
a.h(hist)
a.n(hist)
a.ns(hist)
}
\keyword{classes}
