\name{EM.Control-class}
\docType{class}
\alias{EM.Control-class}
\alias{show,EM.Control-method}
\title{Class \code{"EM.Control"}}
\description{Object of class \code{EM.Control}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EM.Control", ...)}. Accessor methods for the slots are \code{a.strategy(x = NULL)}, \code{a.variant(x = NULL)},
\code{a.acceleration(x = NULL)}, \code{a.tolerance(x = NULL)}, \code{a.acceleration.multiplier(x = NULL)} and \code{a.maximum.iterations(x = NULL)} where \code{x} stands for an object of class \code{EM.Control}. Setter methods \code{a.strategy(x = NULL)}, \code{a.variant(x = NULL)},
\code{a.acceleration(x = NULL)}, \code{a.tolerance(x = NULL)}, \code{a.acceleration.multiplier(x = NULL)} and \code{a.maximum.iterations(x = NULL)} are provided to write to \code{strategy}, \code{variant}, \code{acceleration}, \code{tolerance}, \code{acceleration.multiplier} and \code{maximum.iterations} slot respectively.
}
\section{Slots}{
\describe{
\item{\code{strategy}:}{
a character containing the EM and REBMIX strategy. One of \code{"None"}, \code{"Exhaustive"}, \code{"Best"} and \code{"Single"}. The default value is \code{"None"}.  
}
\item{\code{variant}:}{
a character containing the type of the EM algorithm to be used. One of \code{"EM"} of \code{"ECM"}. The default value is \code{"EM"}.
}
\item{\code{acceleration}:}{
a character containing the type of acceleration of the EM iteration increment. One of \code{"Fixed"}, \code{"Line"} or \code{"Golden"}. The default value is \code{"Fixed"}.
}
\item{\code{tolerence}:}{
tolerance value for the EM convergence criteria. The default value is 1e-4
}
\item{\code{acceleration.multiplier}:}{
acceleration.multiplier \eqn{a_{\mathrm{EM}}}, \eqn{1.0 \leq a_{\mathrm{EM}} \leq 2.0}. acceleration.multiplier for the EM step increment. The default value is 1.0.
}
\item{\code{maximum.iterations}:}{
maximum.iterations is positive integer containing the maximum allowed number of iterations of the EM algorithm. The default value is 1000.
}
}
}
\author{Panic Branislav}
\references{
A. P. Dempster et al. Maximum likelihood from incomplete data via the EM algorithm. Journal of the Royal Statistical Society. Series B, 39(1):1-38, 1977. 
\url{https://www.jstor.org/stable/2984875}.\cr\cr
G. Celeux and G. Govaert. A classification EM algorithm for clustering and two stochastic versions, Computational Statistics & Data Analysis, 14(3):315:332, 1992.
\url{https://doi.org/10.1016/0167-9473(92)90042-E}.
}
\examples{
# Inline creation by function new call.

EM <- new("EM.Control", strategy = "Exhaustive", 
  variant = "EM", acceleration = "Fixed", 
  tolerance = 1e-4, acceleration.multiplier = 1.0, 
  maximum.iterations = 1000)

EM

# Creation of EM object with setter functions.

EM <- new("EM.Control")

a.strategy(EM) <- "Exhaustive"
a.variant(EM) <- "EM"
a.acceleration(EM) <- "Fixed"
a.tolerance(EM) <- 1e-4
a.acceleration.multiplier(EM) <- 1.0
a.maximum.iterations(EM) <- 1000

EM
}
\keyword{classes}
