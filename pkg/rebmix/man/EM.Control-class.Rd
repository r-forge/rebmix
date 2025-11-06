\name{EM.Control-class}
\docType{class}
\alias{EM.Control-class}
\alias{show,EM.Control-method}
\title{Class \code{"EM.Control"}}
\description{Object of class \code{EM.Control}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EM.Control", ...)}. Accessor methods for the slots are \code{a.strategy(x = NULL)},
\code{a.variant(x = NULL)}, \code{a.acceleration(x = NULL)}, \code{a.tolerance(x = NULL)}, \code{a.acceleration.multiplier(x = NULL)}, 
\code{a.maximum.iterations(x = NULL)}, \code{a.K(x = NULL)}, \code{a.eliminate.zero.components(x = NULL)}, \code{a.likelihood.tolerance.check(x = NULL)} and \code{a.likelihood.estimation.rule(x = NULL)} where \code{x} stands for an object of 
class \code{EM.Control}. Setter methods \code{a.strategy(x = NULL)}, \code{a.variant(x = NULL)},
\code{a.acceleration(x = NULL)}, \code{a.tolerance(x = NULL)}, \code{a.acceleration.multiplier(x = NULL)}, \code{a.maximum.iterations(x = NULL)}, 
\code{a.K(x = NULL)}, \code{a.eliminate.zero.components(x = NULL)}, \code{a.likelihood.tolerance.check{x = NULL}} and \code{a.likelihood.estimation.rule{x = NULL}} are provided to write to \code{strategy}, \code{variant}, \code{acceleration}, \code{tolerance}, 
\code{acceleration.multiplier}, \code{maximum.iterations}, \code{eliminate.zero.components}, \code{likelihood.tolerance.check} and \code{likelihood.estimation.rule} slot respectively.
}
\section{Slots}{
\describe{
\item{\code{strategy}:}{
a character containing the EM and REBMIX strategy. One of \code{"none"}, \code{"exhaustive"}, \code{"best"} and \code{"single"}. The default value is \code{"none"}.
}
\item{\code{variant}:}{
a character containing the type of the EM algorithm to be used. One of \code{"EM"}, \code{"ECM"}, \code{"SEM"}, \code{"SEM-EM"}, \code{"ECM-EM"}. The default value is \code{"EM"}.
}
\item{\code{acceleration}:}{ 
a character containing the type of acceleration of the EM iteration increment. One of \code{"fixed"}, \code{"line"}, \code{"golden"}, \code{"stem1"}, \code{"stem2"}, \code{"stem3"}, \code{"square1"}, \code{"square2"} or \code{"square3"}. The default value is \code{"fixed"}.
}
\item{\code{tolerance}:}{
tolerance value for the EM convergence criteria. The default value is 1.0E-4.
}
\item{\code{acceleration.multiplier}:}{
acceleration.multiplier \eqn{a_{\mathrm{EM}}}, \eqn{1.0 \leq a_{\mathrm{EM}} \leq 2.0}. acceleration.multiplier for the EM step increment. The default value is 1.0. Only used when \code{acceleration == "fixed"}.
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
\item{\code{likelihood.tolerance.check}:}{
a character containing the type of log likelihood convergence check. One of \code{"absolute"}, \code{"normalised"} and \code{"percentage"}. The default value is \code{"normalised"}.
}
\item{\code{likelihood.estimation.rule}:}{
a character containing the type of log likelihood estimate. One of \code{"standard"}, \code{"aitken-lindsay"}, \code{"aitken-bohning"} and \code{"aitken-nicholas"}. The default value is \code{"standard"}.
}
}
}
\author{Branislav Panic}
\references{
B. Panic, J. Klemenc, M. Nagode. Improved initialization of the EM algorithm for mixture model parameter estimation.
Mathematics, 8(3):373, 2020.
\doi{10.3390/math8030373}.\cr\cr
A. P. Dempster et al. Maximum likelihood from incomplete data via the EM algorithm. Journal of the Royal Statistical Society. Series B, 39(1):1-38, 1977.
https://www.jstor.org/stable/2984875.\cr\cr
G. Celeux and G. Govaert. A classification EM algorithm for clustering and two stochastic versions, Computational Statistics & Data Analysis, 14(3):315:332, 1992.
\doi{10.1016/0167-9473(92)90042-E}.\cr\cr
R. Varadhan and C. Roland. Simple and Globally Convergent Methods for Accelerating the Convergence of Any EM Algorithm. Scandinavian Journal of Statistics, 35(2):335:353, 2008. 
\doi{https://doi.org/10.1111/j.1467-9469.2007.00585.x}.\cr\cr
P. D. McNicholas et al. Serial and parallel implementations of model-based clustering via parsimonious Gaussian mixture models. Computational Statistics & Data Analysis, 54(3):711:723. 2010.
\doi{https://doi.org/10.1016/j.csda.2009.02.011}
}
\examples{
# Inline creation by new call.

EM <- new("EM.Control", strategy = "exhaustive", 
  variant = "EM", acceleration = "fixed", 
  tolerance = 1e-4, acceleration.multiplier = 1.0, 
  maximum.iterations = 1000, K = 0, likelihood.tolerance.check = "absolute",
  likelihood.estimation.rule = "standard")

EM

# Creation of EM object with setter method.

EM <- new("EM.Control")

a.strategy(EM) <- "exhaustive"
a.variant(EM) <- "EM"
a.acceleration(EM) <- "fixed"
a.tolerance(EM) <- 1e-4
a.acceleration.multiplier(EM) <- 1.0
a.maximum.iterations(EM) <- 1000
a.K(EM) <- 256
a.likelihood.tolerance.check(EM) <- "normalised"
a.likelihood.estimation.rule(EM) <- "standard"

EM
}
\keyword{classes}
