\name{REBMIX.boot-class}
\docType{class}
\alias{REBMIX.boot-class}
\alias{REBMVNORM.boot-class}
\title{Class \code{"REBMIX.boot"}}
\description{Object of class \code{REBMIX.boot}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("REBMIX.boot", ...)}. Accessor methods for the slots are \code{a.rseed(x = NULL)},
\code{a.pos(x = NULL)}, \code{a.Bootstrap(x = NULL)}, \code{a.B(x = NULL)}, \code{a.n(x = NULL)}, \code{a.replace(x = NULL)}, \code{a.prob(x = NULL)},
\code{a.c(x = NULL)}, \code{a.c.se(x = NULL)}, \code{a.c.cv(x = NULL)}, \code{a.c.mode(x = NULL)}, \code{a.c.prob(x = NULL)}, \code{a.w(x = NULL)},
\code{a.w.se(x = NULL)}, \code{a.w.cv(x = NULL)}, \code{a.Theta(x = NULL)}, \code{a.Theta.se(x = NULL)} and \code{a.Theta.cv(x = NULL)}, where \code{x} stands for an object of class \code{REBMIX.boot}. 
}
\section{Slots}{
\describe{
\item{\code{x}:}{
an object of class \code{REBMIX}.
}
\item{\code{rseed}:}{
set the random seed to any negative integer value to initialize the sequence. The first bootstrap dataset corresponds to it.
For each next bootstrap dataset the random seed is decremented \eqn{r_{\mathrm{seed}} = r_{\mathrm{seed}} - 1}. The default value is \code{-1}.
}
\item{\code{pos}:}{
a desired row number in \code{x@summary} to be bootstrapped. The default value is \code{1}.
}
\item{\code{Bootstrap}:}{
a character giving the bootstrap type. One of default \code{"parametric"} or \code{"nonparametric"}.
}
\item{\code{B}:}{
number of bootstrap datasets. The default value is \code{100}.
}
\item{\code{n}:}{
number of observations. The default value is \code{numeric()}.
}
\item{\code{replace}:}{
logical. The sampling is with replacement if \code{TRUE}, see also \code{\link[base]{sample}}. The default value is \code{TRUE}.
}
\item{\code{prob}:}{
a vector of length \eqn{n} containing probability weights, see also \code{\link[base]{sample}}. The default value is \code{numeric()}.
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
a list of matrices containing component parameters \code{theta1.l} and \code{theta2.l} for \eqn{\leq B} bootstrap datasets.
}
\item{\code{Theta.se}:}{
a list of vectors containing standard errors of component parameters \code{theta1.l} and \code{theta2.l}.
}
\item{\code{Theta.cv}:}{
a list of vectors containing coefficients of variation of component parameters \code{theta1.l} and \code{theta2.l}.
}
}
}
\author{Marko Nagode}
\keyword{classes}
