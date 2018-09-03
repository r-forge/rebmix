\name{RCLSMIX-class}
\docType{class}
\alias{RCLSMIX-class}
\alias{RCLSMVNORM-class}
\title{Class \code{"RCLSMIX"}}
\description{Object of class \code{RCLSMIX}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RCLSMIX", ...)}. Accessor methods for the slots are \code{a.o(x = NULL)}, ..., where \code{x} stands for an object of class \code{RCLSMIX}. 
}
\section{Slots}{
\describe{
\item{\code{x}:}{
a list of objects of class \code{REBMIX} of length \eqn{o} obtained by running \code{\link[rebmix]{REBMIX}} on \eqn{g = 1, \ldots, s} train datasets \eqn{Y_{\mathrm{train}g}} all of length \eqn{n_{\mathrm{train}g}}.
For the train datasets the corresponding class membership \eqn{\bm{\Omega}_{g}} is known. This yields
\eqn{n_{\mathrm{train}} = \sum_{g = 1}^{s} n_{\mathrm{train}g}}, while \eqn{Y_{\mathrm{train}q} \cap Y_{\mathrm{train}g} = \emptyset} for all \eqn{q \neq g}.
Each object in the list corresponds to one chunk, e.g., \eqn{(y_{1j}, y_{3j})^{\top}}.
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
a vector of length \eqn{s} containing numbers of observations in train datasets \eqn{Y_{\mathrm{train}g}}.
}
\item{\code{P}:}{
a vector of length \eqn{s} containing prior probabilities \eqn{P(\bm{\Omega}_{g}) = \frac{n_{\mathrm{train}g}}{n_{\mathrm{train}}}}.
}
\item{\code{ntest}:}{
number of observations in test dataset \eqn{Y_{\mathrm{test}}}.
}
\item{\code{Zt}:}{
a factor of true class membership \eqn{\bm{\Omega}_{g}} for the test dataset.
}
\item{\code{Zp}:}{
a factor of predictive class membership \eqn{\bm{\Omega}_{g}} for the test dataset.
}
\item{\code{CM}:}{
a table containing confusion matrix for multiclass classifier. It contains 
number \eqn{x_{qg}} of test observations with the true class \eqn{q} that are classified into the class \eqn{g}, where \eqn{q, g = 1, \ldots, s}.
}
\item{\code{Accuracy}:}{
proportion of all test observations that are classified correctly. \eqn{\mathrm{Accuracy} = \frac{\sum_{g = 1}^{s} x_{gg}}{n_{\mathrm{test}}}}.
}
\item{\code{Error}:}{
proportion of all test observations that are classified wrongly. \eqn{\mathrm{Error} = 1 - \mathrm{Accuracy}}.
}
\item{\code{Precision}:}{
a vector containing proportions of predictive observations in class \eqn{g} that are
classified correctly into class \eqn{g}. \eqn{\mathrm{Precision}(g) = \frac{x_{gg}}{\sum_{q = 1}^{s} x_{qg}}}.
}
\item{\code{Sensitivity}:}{
a vector containing proportions of test observations in class \eqn{g} that are classified
correctly into class \eqn{g}. \eqn{\mathrm{Sensitivity}(g) = \frac{x_{gg}}{\sum_{q = 1}^{s} x_{gq}}}.
}
\item{\code{Specificity}:}{
a vector containing proportions of test observations that are not in class \eqn{g} and
are classified into the non \eqn{g} class. \eqn{\mathrm{Specificity}(g) = \frac{n_{\mathrm{test}} - \sum_{q = 1}^{s} x_{qg}}{n_{\mathrm{test}} - \sum_{q = 1}^{s} x_{gq}}}.
}
\item{\code{Chunks}:}{
a vector containing selected chunks.
}
}
}
\author{Marko Nagode}
\references{
D. M. Dziuda. Data Mining for Genomics and Proteomics: Analysis of Gene and Protein
Expression Data. John Wiley & Sons, New York, 2010.
}
\keyword{classes}
