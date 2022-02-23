\name{Histogram-class}
\docType{class}
\alias{Histogram-class}
\title{Class \code{"Histogram"}}
\description{Object of class \code{Histogram}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Histogram", ...)}. Accessor methods for the slots are \code{a.Dataset.name(x = NULL)},
\code{a.Dataset(x = NULL)} and \code{a.h(x = NULL)}.
}
\section{Slots}{
\describe{
\item{\code{Dataset.name}:}{
a character containing the name of data frame of size \eqn{v \times (d + 1)} that \emph{d}-dimensional histogram is written in.
}
\item{\code{Dataset}:}{
a data frame of size \eqn{v \times (d + 1)} containing \emph{d}-dimensional histogram. Each of the \eqn{d}
columns represents one random variable and contains the coresponding bin means. Column \eqn{d + 1} contains frequencies.
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
