\name{RCLS.chunk-class}
\docType{class}
\alias{RCLS.chunk-class}
\title{Class \code{"RCLS.chunk"}}
\description{Object of class \code{RCLS.chunk}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RCLS.chunk", ...)}. Accessor methods for the slots are \code{a.s(x = NULL)},
\code{a.levels(x = NULL)}, \code{a.ntrain(x = NULL)}, \code{a.train(x = NULL)}, \code{a.Zr(x = NULL)}, \code{a.ntest(x = NULL)}, \code{a.test(x = NULL)} and \code{a.Zt(x = NULL)},
where \code{x} stands for an object of class \code{RCLS.chunk}.
}
\section{Slots}{
\describe{
\item{\code{s}:}{
finite set of size \eqn{s} of classes \eqn{\bm{\Omega} = \{\bm{\Omega}_{g}; \ g = 1, \ldots, s\}}.
}
\item{\code{levels}:}{
a character vector of length \eqn{s} containing class names \eqn{\bm{\Omega}_{g}}.
}
\item{\code{ntrain}:}{
a vector of length \eqn{s} containing numbers of observations in train datasets \eqn{Y_{\mathrm{train}g}}.
}
\item{\code{train}:}{
a list of length \eqn{n_{\mathrm{D}}} of data frames containing train datasets \eqn{Y_{\mathrm{train}g}} of length \eqn{n_{\mathrm{train}g}}.
}
\item{\code{Zr}:}{
a list of factors of true class membership \eqn{\bm{\Omega}_{g}} for the train datasets.
}
\item{\code{ntest}:}{
number of observations in test dataset \eqn{Y_{\mathrm{test}}}.
}
\item{\code{test}:}{
a data frame containing test dataset \eqn{Y_{\mathrm{test}}} of length \eqn{n_{\mathrm{test}}}.
}
\item{\code{Zt}:}{
a factor of true class membership \eqn{\bm{\Omega}_{g}} for the test dataset.
}
}
}
\author{Marko Nagode}
\references{
D. M. Dziuda. Data Mining for Genomics and Proteomics: Analysis of Gene and Protein
Expression Data. John Wiley & Sons, New York, 2010.
}
\keyword{classes}
