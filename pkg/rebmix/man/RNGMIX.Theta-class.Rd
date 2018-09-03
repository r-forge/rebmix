\name{RNGMIX.Theta-class}
\docType{class}
\alias{RNGMIX.Theta-class}
\alias{RNGMVNORM.Theta-class}
\alias{show,RNGMIX.Theta-method}
\alias{show,RNGMVNORM.Theta-method}
\title{Class \code{"RNGMIX.Theta"}}
\description{Object of class \code{RNGMIX.Theta}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RNGMIX.Theta", ...)}. Accessor methods for the slots are \code{a.c.name(x = NULL)}, ...,
\code{a.Theta(x = NULL, pos = 0)}, where \code{x} and \code{pos} stand for an object of class \code{RNGMIX} and a desired slot item, respectively. Setter methods
\code{a.theta1(x = NULL, l = numeric())} and \code{a.theta2(x = NULL, l = numeric())} are provided to write to \code{Theta} slot, where \eqn{l = 1, \ldots, c}.
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
Theta <- new("RNGMIX.Theta", c = 2, pdf = c("normal", "Weibull"))  

a.theta1(Theta, l = 1) <- c(2, 10)
a.theta2(Theta, l = 1) <- c(0.5, 2.3)
a.theta1(Theta, l = 2) <- c(20, 50)
a.theta2(Theta, l = 2) <- c(3, 4.2)

Theta

Theta <- new("RNGMVNORM.Theta", c = 2, d = 3)

a.theta1(Theta, l = 1) <- c(2, 10, -20)
a.theta1(Theta, l = 2) <- c(-2.4, -15.1, 30)

Theta
}
\keyword{classes}
