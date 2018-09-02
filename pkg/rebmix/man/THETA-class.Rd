\name{THETA-class}
\docType{class}
\alias{THETA-class}
\title{Class \code{"THETA"}}
\description{Object of class \code{THETA}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("THETA", ...)}. Accessor functions for the slots are \code{g.c.name(x = NULL)}, ...,
\code{g.Theta(x = NULL, pos = 0)}, where \code{x} and \code{pos} stand for an object of class \code{RNGMIX} and a desired slot item, respectively.
}
\section{Slots}{
\describe{
\item{c}{
number of components \eqn{c > 0}. The default value is \code{1}.
}
\item{d}{
number of dimensions.
}
\item{\code{pdf}:}{
a character vector of length \eqn{d} containing continuous or discrete parametric family types. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"}, \code{"Dirac"} or \code{"vonMises"}.
}
\item{\code{Theta}:}{
a list containing \eqn{c} parametric family types \code{pdfl}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"}, \code{"Dirac"} or circular \code{"vonMises"} defined for \eqn{0 \leq y_{i} \leq 2 \pi}.
Component parameters \code{theta1.l} follow the parametric family types. One of \eqn{\mu_{il}} for normal, lognormal and von Mises distributions and \eqn{\theta_{il}} for Weibull, gamma, binomial, Poisson and Dirac distributions.
Component parameters \code{theta2.l} follow \code{theta1.l}. One of \eqn{\sigma_{il}} for normal and lognormal distributions, \eqn{\beta_{il}} for Weibull and gamma distributions, \eqn{p_{il}} for binomial distribution, \eqn{\kappa_{il}} for von Mises distribution and \code{"NA"} otherwise.
}
}
}
\author{Marko Nagode}
\examples{
Theta <- new("THETA", c = 5, pdf = rep("n", 10))  
}
\keyword{classes}
