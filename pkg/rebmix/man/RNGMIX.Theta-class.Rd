\name{RNGMIX.Theta-class}
\docType{class}
\alias{RNGMIX.Theta-class}
\alias{RNGMVNORM.Theta-class}
\alias{show,RNGMIX.Theta-method}
\alias{show,RNGMVNORM.Theta-method}
\title{Class \code{"RNGMIX.Theta"}}
\description{Object of class \code{RNGMIX.Theta}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RNGMIX.Theta", ...)}. Accessor methods for the slots are \code{a.c(x = NULL)}, \code{a.d(x = NULL)},
\code{a.pdf(x = NULL)} and \code{a.Theta(x = NULL)}, where \code{x} stands for an object of class \code{RNGMIX.Theta}. Setter methods
\code{a.theta1(x = NULL, l = numeric())}, \code{a.theta2(x = NULL, l = numeric())} and \code{a.theta3(x = NULL, l = numeric())}, 
\code{a.theta1.all(x = NULL)}, \code{a.theta2.all(x = NULL)} and \code{a.theta3.all(x = NULL)} 
are provided to write to \code{Theta} slot, where \eqn{l = 1, \ldots, c}.
}
\section{Slots}{
\describe{
\item{\code{c}:}{
number of components \eqn{c > 0}. The default value is \code{1}.
}
\item{\code{d}:}{
number of dimensions.
}
\item{\code{pdf}:}{
a character vector of length \eqn{d} containing continuous or discrete parametric family types. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"Gumbel"}, \code{"binomial"}, \code{"Poisson"}, \code{"Dirac"}, \code{"uniform"} or \code{"vonMises"}.
}
\item{\code{Theta}:}{
a list containing \eqn{c} parametric family types \code{pdfl}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"Gumbel"}, \code{"binomial"}, \code{"Poisson"}, \code{"Dirac"}, \code{"uniform"} or circular \code{"vonMises"} defined for \eqn{0 \leq y_{i} \leq 2 \pi}.
Component parameters \code{theta1.l} follow the parametric family types. One of \eqn{\mu_{il}} for normal, lognormal, Gumbel and von Mises distributions, \eqn{\theta_{il}} for Weibull, gamma, binomial, Poisson and Dirac distributions and \eqn{a} for uniform distribution.
Component parameters \code{theta2.l} follow \code{theta1.l}. One of \eqn{\sigma_{il}} for normal, lognormal and Gumbel distributions, \eqn{\beta_{il}} for Weibull and gamma distributions, \eqn{p_{il}} for binomial distribution, \eqn{\kappa_{il}} for von Mises distribution and \eqn{b} for uniform distribution.
Component parameters \code{theta3.l} follow \code{theta2.l}. One of \eqn{\xi_{il} \in \{-1, 1\}} for Gumbel distribution.
}
}
}
\author{Marko Nagode}
\examples{
Theta <- new("RNGMIX.Theta", c = 2, pdf = c("normal", "Gumbel"))  

a.theta1(Theta, l = 1) <- c(2, 10)
a.theta2(Theta, l = 1) <- c(0.5, 2.3)
a.theta3(Theta, l = 1) <- c(NA, 1.0)
a.theta1(Theta, l = 2) <- c(20, 50)
a.theta2(Theta, l = 2) <- c(3, 4.2)
a.theta3(Theta, l = 2) <- c(NA, -1.0)

Theta

Theta <- new("RNGMIX.Theta", c = 2, pdf = c("normal", "Gumbel"))  

a.theta1.all(Theta) <- c(2, 10, 20, 50)
a.theta2.all(Theta) <- c(0.5, 2.3, 3, 4.2)
a.theta3.all(Theta) <- c(NA, 1.0, NA, -1.0)

Theta

Theta <- new("RNGMVNORM.Theta", c = 2, d = 3)

a.theta1(Theta, l = 1) <- c(2, 10, -20)
a.theta1(Theta, l = 2) <- c(-2.4, -15.1, 30)

Theta
}
\keyword{classes}
