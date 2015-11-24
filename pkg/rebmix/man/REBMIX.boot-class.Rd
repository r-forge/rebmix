\name{REBMIX.boot-class}
\docType{class}
\alias{REBMIX.boot-class}
\alias{REBMVNORM.boot-class}
\title{Class "REBMIX.boot"}
\description{Object of class \code{REBMIX.boot}.}
\section{Objects from the Class}{
  Objects can be created by \code{new("REBMIX.boot", ...)}. 
}
\section{Slots}{
\describe{
\item{\code{x}:}{
an object of class \code{REBMIX}.
}
\item{\code{pos}:}{
a desired row number in \code{x@summary} to be bootstrapped. The default value is \code{1}.
}
\item{\code{Bootstrap}:}{
a character vector giving the bootstrap type. One of default \code{"parametric"} or \code{"nonparametric"}.
}
\item{\code{B}:}{
number of bootstrap datasets. The default value is \code{100}.
}
\item{\code{n}:}{
number of observations. The default value is \code{numeric(0)}.
}
\item{\code{replace}:}{
logical. The sampling is with replacement if \code{TRUE}, see also \code{\link[base]{sample}}. The default value is \code{TRUE}.
}
\item{\code{prob}:}{
a vector of length \eqn{n} containing probability weights, see also \code{\link[base]{sample}}. The default value is \code{numeric(0)}.
}
\item{\code{c}:}{
a vector containing numbers of components for \eqn{B} bootstrap datasets.
}
\item{\code{c.se}:}{
standard error of numbers of components \code{c}.
}
\item{\code{c.cv}:}{
coefficient of variation of numbers of components \code{c}.
}
\item{\code{c.mode}:}{
mode of numbers of components \code{c}.
}
\item{\code{c.prob}:}{
probability of mode \code{c.mode}.
}
\item{\code{w}:}{
a matrix containing component weights for \eqn{\leq B} bootstrap datasets.
}
\item{\code{w.se}:}{
a vector containing standard errors of component weights \code{w}.
}
\item{\code{w.cv}:}{
a vector containing coefficients of variation of component weights \code{w}.
}
\item{\code{Theta}:}{
a list of matrices containing component parameters \code{theta1.i} and \code{theta2.i} for \eqn{\leq B} bootstrap datasets.
}
\item{\code{Theta.se}:}{
a list of vectors containing standard errors of component parameters \code{theta1.i} and \code{theta2.i}.
}
\item{\code{Theta.cv}:}{
a list of vectors containing coefficients of variation of component parameters \code{theta1.i} and \code{theta2.i}.
}
}
}
\author{Marko Nagode}
\keyword{classes}
