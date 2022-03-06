\name{EM.Control-class}
\docType{class}
\alias{EM.Control-class}
\alias{show,EM.Control-method}
\title{Class \code{"EM.Control"}}
\description{Object of class \code{EM.Control}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EM.Control", ...)}. Accessor methods for the slots are \code{a.strategy(x = NULL)},
\code{a.variant(x = NULL)}, \code{a.acceleration(x = NULL)}, \code{a.tolerance(x = NULL)}, \code{a.acceleration.multiplier(x = NULL)}, 
\code{a.maximum.iterations(x = NULL)}, \code{a.K(x = NULL)} and \code{a. eliminate.zero.components (x = NULL)}, where \code{x} stands for an object of 
class \code{EM.Control}. Setter methods \code{a.strategy(x = NULL)}, \code{a.variant(x = NULL)},
\code{a.acceleration(x = NULL)}, \code{a.tolerance(x = NULL)}, \code{a.acceleration.multiplier(x = NULL)}, \code{a.maximum.iterations(x = NULL)}, 
\code{a.K(x = NULL)} and \code{eliminate.zero.components} are provided to write to \code{strategy}, \code{variant}, \code{acceleration}, \code{tolerance}, 
\code{acceleration.multiplier}, \code{maximum.iterations} and \code{eliminate.zero.components} slot respectively.
}
\section{Slots}{
\describe{
\item{\code{strategy}:}{
a character containing the EM and REBMIX strategy. One of \code{"none"}, \code{"exhaustive"}, \code{"best"} and \code{"single"}. The default value is \code{"none"}.
}
\item{\code{variant}:}{
a character containing the type of the EM algorithm to be used. One of \code{"EM"} of \code{"ECM"}. The default value is \code{"EM"}.
}
\item{\code{acceleration}:}{
a character containing the type of acceleration of the EM iteration increment. One of \code{"fixed"}, \code{"line"} or \code{"golden"}. The default value is \code{"fixed"}.
}
\item{\code{tolerance}:}{
tolerance value for the EM convergence criteria. The default value is 1.0E-4.
}
\item{\code{acceleration.multiplier}:}{
acceleration.multiplier \eqn{a_{\mathrm{EM}}}, \eqn{1.0 \leq a_{\mathrm{EM}} \leq 2.0}. acceleration.multiplier for the EM step increment. The default value is 1.0.
}
\item{\code{maximum.iterations}:}{
a positive integer containing the maximum allowed number of iterations of the EM algorithm. The default value is 1000.
}
\item{\code{K}:}{
an integer containing the number of bins for the histogram based EM algorithm. This option can reduce computational time drastically if the datasets contain a large number of observations \eqn{n} and \code{K} is set to the value \eqn{\ll n}. The default value of 0 means that the EM algorithm runs over all \eqn{n}.
}
\item{\code{eliminate.zero.components}:}{
a logical indicating if the componenets with \eqn{w_{l} = 0} should be eliminated from output. Only used with \code{EMMIX-methods}.
}
}
}
\author{Branislav Panic}
\references{
B. Panic, J. Klemenc, M. Nagode. Improved initialization of the EM algorithm for mixture model parameter estimation.
Mathematics, 8(3):373, 2020.
\doi{10.3390/math8030373}.\cr\cr
A. P. Dempster et al. Maximum likelihood from incomplete data via the EM algorithm. Journal of the Royal Statistical Society. Series B, 39(1):1-38, 1977.
\url{https://www.jstor.org/stable/2984875}.\cr\cr
G. Celeux and G. Govaert. A classification EM algorithm for clustering and two stochastic versions, Computational Statistics & Data Analysis, 14(3):315:332, 1992.
\doi{10.1016/0167-9473(92)90042-E}.
}
\examples{
# Inline creation by function new call.

EM <- new("EM.Control", strategy = "exhaustive", 
  variant = "EM", acceleration = "fixed", 
  tolerance = 1e-4, acceleration.multiplier = 1.0, 
  maximum.iterations = 1000, K = 0)

EM

# Creation of EM object with setter functions.

EM <- new("EM.Control")

a.strategy(EM) <- "exhaustive"
a.variant(EM) <- "EM"
a.acceleration(EM) <- "fixed"
a.tolerance(EM) <- 1e-4
a.acceleration.multiplier(EM) <- 1.0
a.maximum.iterations(EM) <- 1000
a.K(EM) <- 256

EM
}
\keyword{classes}
