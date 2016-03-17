\name{REBMIX-class}
\docType{class}
\alias{REBMIX-class}
\alias{REBMVNORM-class}
\title{Class \code{"REBMIX"}}
\description{Object of class \code{REBMIX}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("REBMIX", ...)}.
}
\section{Slots}{
\describe{
\item{\code{Dataset}:}{
a list of data frames of size \eqn{n \times d} containing \emph{d}-dimensional datasets. Each of the \eqn{d} columns
represents one random variable. Number of observations \eqn{n} equals the number of rows in the datasets.
}
\item{\code{Preprocessing}:}{
a character vector, giving the preprocessing types. One of \code{"histogram"}, \code{"Parzen window"} or \code{"k-nearest neighbour"}.
}
\item{\code{cmax}:}{
maximum number of components \eqn{c_{\mathrm{max}} > 0}. The default value is \code{15}.
}
\item{\code{Criterion}:}{
a character vector giving the infromation criterion types. One of default Akaike \code{"AIC"}, \code{"AIC3"}, \code{"AIC4"} or \code{"AICc"}, 
Bayesian \code{"BIC"}, consistent Akaike \code{"CAIC"}, Hannan-Quinn \code{"HQC"}, minimum description length \code{"MDL2"} or \code{"MDL5"}, 
approximate weight of evidence \code{"AWE"}, classification likelihood \code{"CLC"}, 
integrated classification likelihood \code{"ICL"} or \code{"ICL-BIC"}, partition coefficient \code{"PC"}, 
total of positive relative deviations \code{"D"} or sum of squares error \code{"SSE"}.
}
\item{\code{Variables}:}{
a character vector of length \eqn{d} containing types of variables. One of \code{"continuous"} or \code{"discrete"}.
}
\item{\code{pdf}:}{
a character vector of length \eqn{d} containing continuous or discrete parametric family types. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
}
\item{\code{theta1}:}{
a vector of length \eqn{d} containing initial component parameters. One of \eqn{n_{il} = \textrm{Number of categories} - 1} for \code{"binomial"} distribution or \code{"NA"} otherwise. 
}
\item{\code{theta2}:}{
a vector of length \eqn{d} containing initial component parameters. Currently not used.
}
\item{\code{K}:}{
a vector or a list of vectors containing numbers of bins \eqn{v} for the histogram and the Parzen window or numbers of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour. There is no genuine rule to identify \eqn{v} or \eqn{k}. Consequently,
the REBMIX algorithm identifies them from the set \code{K} of input values by
minimizing the information criterion. The Sturges rule \eqn{v = 1 + \mathrm{log_{2}}(n)}, \eqn{\mathrm{Log}_{10}} rule \eqn{v = 10 \mathrm{log_{10}}(n)} or RootN 
rule \eqn{v = 2 \sqrt{n}} can be applied to estimate the limiting numbers of bins
or the rule of thumb \eqn{k = \sqrt{n}} to guess the intermediate number of nearest neighbours. If, e.g., \code{K = c(10, 20, 40, 60)} and minimum \code{IC} coincides, e.g., \code{40}, brackets are set to \code{20} and \code{60} and the golden section is applied to refine the minimum search. See also \code{\link[rebmix]{kseq}} for sequence of bins or nearest neighbours generation. 
}
\item{\code{y0}:}{
a vector of length \eqn{d} containing origins. The default value is \code{numeric()}.
}
\item{\code{ymin}:}{
a vector of length \eqn{d} containing minimum observations. The default value is \code{numeric()}.
}
\item{\code{ymax}:}{
a vector of length \eqn{d} containing maximum observations. The default value is \code{numeric()}.
}
\item{\code{ar}:}{
acceleration rate \eqn{0 < a_{\mathrm{r}} \leq 1}. The default value is \code{0.1} and in most cases does not have to be altered.
}
\item{\code{Restraints}:}{
a character giving the restraints type. One of \code{"rigid"} or default \code{"loose"}.
The rigid restraints are obsolete and applicable for well separated components only. 
}
\item{\code{w}:}{
a list of vectors of length \eqn{c} containing component weights \eqn{w_{l}} summing to 1.
}
\item{\code{Theta}:}{
a list of lists each containing  \eqn{c} parametric family types \code{pdfi}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
Component parameters \code{theta1.i} follow the parametric family types. One of \eqn{\mu_{il}} for normal and lognormal distributions and \eqn{\theta_{il}} for Weibull, gamma, binomial, Poisson and Dirac distributions.
Component parameters \code{theta2.i} follow \code{theta1.i}. One of \eqn{\sigma_{il}} for normal and lognormal distributions, \eqn{\beta_{il}} for Weibull and gamma distributions and \eqn{p_{il}} for binomial distribution.
}
\item{\code{summary}:}{
a data frame with additional information about dataset, preprocessing, \eqn{c_{\mathrm{max}}}, information criterion type, 
\eqn{a_{\mathrm{r}}}, restraints type, optimal \eqn{c}, optimal \eqn{v} or \eqn{k}, \eqn{K}, \eqn{y_{i0}}, optimal \eqn{h_{i}},
information criterion \eqn{\mathrm{IC}}, log likelihood \eqn{\mathrm{log}\, L} and degrees of freedom \eqn{M}.
}
\item{\code{pos}:}{
position in the \code{summary} data frame at which log likelihood \eqn{\mathrm{log}\, L} attains its maximum.
}
\item{\code{opt.c}:}{
a list of vectors containing numbers of components for optimal \eqn{v} for the histogram and the Parzen window or for optimal number of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
\item{\code{opt.IC}:}{
a list of vectors containing information criteria for optimal \eqn{v} for the histogram and the Parzen window or for optimal number of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
\item{\code{opt.logL}:}{
a list of vectors containing log likelihoods for optimal \eqn{v} for the histogram and the Parzen window or for optimal number of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
\item{\code{opt.D}:}{
a list of vectors containing totals of positive relative deviations for optimal \eqn{v} for the histogram and the Parzen window or for optimal number of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
\item{\code{all.K}:}{
a list of vectors containing all processed numbers of bins \eqn{v} for the histogram and the Parzen window or all processed numbers of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
\item{\code{all.IC}:}{
a list of vectors containing information criteria for all processed numbers of bins \eqn{v} for the histogram and the Parzen window or for all processed numbers of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
}
}
\author{Marko Nagode}
\keyword{classes}
