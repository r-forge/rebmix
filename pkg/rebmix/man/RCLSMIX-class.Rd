\name{RCLSMIX-class}
\docType{class}
\alias{RCLSMIX-class}
\alias{RCLSMVNORM-class}
\title{Class \code{"RCLSMIX"}}
\description{Object of class \code{RCLSMIX}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RCLSMIX", ...)}.
}
\section{Slots}{
\describe{
\item{\code{x}:}{
a list of objects of class \code{REBMIX} of length \eqn{o} obtained by running \code{\link[rebmix]{REBMIX}} on \eqn{g = 1, \ldots, s} train datasets \eqn{Y_{\mathrm{train}g}} all of length \eqn{n_{\mathrm{train}g}}.
For the train datasets the corresponding class membership \eqn{\bm{\Omega}_{g}} is known. This yields
\eqn{n_{\mathrm{train}} = \sum_{g = 1}^{s} n_{\mathrm{train}g}}, while \eqn{Y_{\mathrm{train}q} \cap Y_{\mathrm{train}g} = \emptyset} for all \eqn{q \neq g}.
}
\item{\code{o}:}{
number of chunks \eqn{o}. \eqn{Y = \{\bm{y}_{j}; \ j = 1, \ldots, n\}} is an observed \eqn{d}-dimensional dataset of size \eqn{n} of vector observations \eqn{\bm{y}_{j} = (y_{1j}, \ldots, y_{dj})^{\top}} and
is partitioned into train and test datasets. Vector observations \eqn{\bm{y}_{j}} may further be split into \eqn{o} chunks when running \code{\link[rebmix]{REBMIX}}, e.g.,
for \eqn{d = 6} and \eqn{o = 3} the set of chunks substituting \eqn{\bm{y}_{j}} may be as follows \eqn{(y_{1j}, y_{3j})^{\top}}, \eqn{(y_{2j}, y_{4j}, y_{6j})^{\top}} and \eqn{y_{5j}}.
}
\item{\code{Dataset}:}{
a data frame containing test dataset \eqn{Y_{\mathrm{test}}} of length \eqn{n_{\mathrm{test}}}. For the test dataset the corresponding class membership \eqn{\bm{\Omega}_{g}} is not known.
}
\item{\code{s}:}{
finite set of size \eqn{s} of classes \eqn{\bm{\Omega} = \{\bm{\Omega}_{g}; \ g = 1, \ldots, s\}}.
}
\item{\code{ntrain}:}{
a vector of length \eqn{s} containing number of observations in train datasets \eqn{Y_{\mathrm{train}g}}.
}
\item{\code{P}:}{
a vector of length \eqn{s} containing prior probabilities \eqn{P(\bm{\Omega}_{g}) = \frac{n_{\mathrm{train}g}}{n_{\mathrm{train}}}}.
}
\item{\code{ntest}:}{
number of observations in test dataset \eqn{Y_{\mathrm{test}}}.
}
\item{\code{Z}:}{
predictive classes.
}
}
}
\author{Marko Nagode}
\keyword{classes}
