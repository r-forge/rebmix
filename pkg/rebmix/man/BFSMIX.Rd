\name{BFSMIX-methods}
\docType{methods}
\alias{BFSMIX}
\alias{BFSMIX-methods}
\alias{BFSMIX,RCLSMIX-method}
\alias{BFSMIX,RCLSMVNORM-method}
\title{
Predicts Class Membership Based Upon the Best First Search Algorithm
}
\description{
Returns as default the optimized RCLSMIX algorithm output for mixtures of conditionally independent normal, lognormal, Weibull, gamma, Gumbel, binomial, Poisson, Dirac or von Mises component densities. If \code{model} equals \code{"RCLSMVNORM"} optimized output for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices is returned.
}
\usage{
\S4method{BFSMIX}{RCLSMIX}(model = "RCLSMIX", x = list(), Dataset = data.frame(),
       Zt = factor(), ...)
## ... and for other signatures
}
\arguments{
\item{model}{
see Methods section below.
}
\item{x}{
a list of objects of class \code{REBMIX} of length \eqn{o} obtained by running \code{\link[rebmix]{REBMIX}} on \eqn{g = 1, \ldots, s} train datasets \eqn{Y_{\mathrm{train}g}} all of length \eqn{n_{\mathrm{train}g}}.
For the train datasets the corresponding class membership \eqn{\bm{\Omega}_{g}} is known. This yields
\eqn{n_{\mathrm{train}} = \sum_{g = 1}^{s} n_{\mathrm{train}g}}, while \eqn{Y_{\mathrm{train}q} \cap Y_{\mathrm{train}g} = \emptyset} for all \eqn{q \neq g}.
Each object in the list corresponds to one chunk, e.g., \eqn{(y_{1j}, y_{3j})^{\top}}. The default value is \code{list()}.
}
\item{Dataset}{
a data frame containing test dataset \eqn{Y_{\mathrm{test}}} of length \eqn{n_{\mathrm{test}}}. For the test dataset the corresponding class membership \eqn{\bm{\Omega}_{g}} is not known.
The default value is \code{data.frame()}.
}
\item{Zt}{
a factor of true class membership \eqn{\bm{\Omega}_{g}} for the test dataset. The default value is \code{factor()}.
}
\item{\dots}{
currently not used.
}
}
\value{
Returns an optimized object of class \code{RCLSMIX} or \code{RCLSMVNORM}.
}
\section{Methods}{
\describe{
\item{\code{signature(model = "RCLSMIX")}}{a character giving the default class name \code{"RCLSMIX"} for mixtures of conditionally independent normal, lognormal, Weibull, gamma, Gumbel, binomial, Poisson, Dirac or von Mises component densities.}
\item{\code{signature(model = "RCLSMVNORM")}}{a character giving the class name \code{"RCLSMVNORM"} for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices.}
}
}
\author{Marko Nagode}
\references{
R. Kohavi and G. H. John. Wrappers for feature subset selection, Artificial Intelligence, 97(1-2):273-324, 1997.
\url{https://doi.org/10.1016/S0004-3702(97)00043-X}.
}
\keyword{classification}
