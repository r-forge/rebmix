\name{REBMIX-class}
\docType{class}
\alias{REBMIX-class}
\alias{REBMVNORM-class}
\title{Class \code{"REBMIX"}}
\description{Object of class \code{REBMIX}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("REBMIX", ...)}. Accessor methods for the slots are \code{a.Dataset(x = NULL, pos = 0)},
\code{a.Preprocessing(x = NULL)}, \code{a.cmax(x = NULL)}, \code{a.cmin(x = NULL)}, \code{a.Criterion(x = NULL)}, \code{a.Variables(x = NULL)},
\code{a.pdf(x = NULL)}, \code{a.theta1(x = NULL)}, \code{a.theta2(x = NULL)}, \code{a.theta3(x = NULL)}, \code{a.K(x = NULL)}, \code{a.y0(x = NULL)}, \code{a.ymin(x = NULL)},
\code{a.ymax(x = NULL)}, \code{a.ar(x = NULL)}, \code{a.Restraints(x = NULL)}, \code{a.w(x = NULL, pos = 0)}, \code{a.Theta(x = NULL, pos = 0)}, \code{a.summary(x = NULL, col.name = character(), pos = 0)},
\code{a.summary.EM(x = NULL, col.name = character(), pos = 0)}, \code{a.pos(x = NULL)},
\code{a.opt.c(x = NULL)}, \code{a.opt.IC(x = NULL)}, \code{a.opt.logL(x = NULL)}, \code{a.opt.D(x = NULL)}, \code{a.all.K(x = NULL)}, \code{a.all.IC(x = NULL)},
\code{a.theta1.all(x = NULL, pos = 1)}, \code{a.theta2.all(x = NULL, pos = 1)} and \code{a.theta3.all(x = NULL, pos = 1)}, where \code{x}, \code{pos} and \code{col.name} stand for an object of class \code{REBMIX},
a desired slot item and a desired column name, respectively.
}
\section{Slots}{
\describe{
\item{\code{Dataset}:}{
a list of length \eqn{n_{\mathrm{D}}} of data frames or objects of class \code{Histogram}. 
Data frames should have size \eqn{n \times d} containing \emph{d}-dimensional datasets. Each of the \eqn{d}
columns represents one random variable. Numbers of observations \eqn{n} equal the number of rows in the datasets.
}
\item{\code{Preprocessing}:}{
a character vector giving the preprocessing types. One of \code{"histogram"}, \cr
\code{"kernel density estimation"} or \code{"k-nearest neighbour"}.
}
\item{\code{cmax}:}{
maximum number of components \eqn{c_{\mathrm{max}} > 0}. The default value is \code{15}.
}
\item{\code{cmin}:}{
minimum number of components \eqn{c_{\mathrm{min}} > 0}. The default value is \code{1}.
}
\item{\code{Criterion}:}{
a character giving the information criterion type. One of default Akaike \code{"AIC"}, \code{"AIC3"}, \code{"AIC4"} or \code{"AICc"},
Bayesian \code{"BIC"}, consistent Akaike \code{"CAIC"}, Hannan-Quinn \code{"HQC"}, minimum description length \code{"MDL2"} or \code{"MDL5"},
approximate weight of evidence \code{"AWE"}, classification likelihood \code{"CLC"},
integrated classification likelihood \code{"ICL"} or \code{"ICL-BIC"}, partition coefficient \code{"PC"},
total of positive relative deviations \code{"D"} or sum of squares error \code{"SSE"}.
}
\item{\code{Variables}:}{
a character vector of length \eqn{d} containing types of variables. One of \code{"continuous"} or \code{"discrete"}.
}
\item{\code{pdf}:}{
a character vector of length \eqn{d} containing continuous or discrete parametric family types. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"Gumbel"}, \code{"binomial"}, \code{"Poisson"}, \code{"Dirac"}, \code{"uniform"} or \code{"vonMises"}.
}
\item{\code{theta1}:}{
a vector of length \eqn{d} containing initial component parameters. One of \eqn{n_{il} = \textrm{number of categories} - 1} for \code{"binomial"} distribution.
}
\item{\code{theta2}:}{
a vector of length \eqn{d} containing initial component parameters. Currently not used.
}
\item{\code{theta3}:}{
a vector of length \eqn{d} containing initial component parameters. One of \eqn{\xi_{il} \in \{-1, \textrm{NA}, 1\}} for \code{"Gumbel"} distribution.
}
\item{\code{K}:}{
a character or a vector or a list of vectors containing numbers of bins \eqn{v} for the histogram and the kernel density estimation or numbers of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour. There is no genuine rule to identify \eqn{v} or \eqn{k}. Consequently,
the REBMIX algorithm identifies them from the set \code{K} of input values by
minimizing the information criterion. The Sturges rule \eqn{v = 1 + \mathrm{log_{2}}(n)}, \eqn{\mathrm{Log}_{10}} rule \eqn{v = 10 \mathrm{log_{10}}(n)} or RootN
rule \eqn{v = 2 \sqrt{n}} can be applied to estimate the limiting numbers of bins
or the rule of thumb \eqn{k = \sqrt{n}} to guess the intermediate number of nearest neighbours. If, e.g., \code{K = c(10, 20, 40, 60)} and minimum \code{IC} coincides, e.g., \code{40}, brackets are set to \code{20} and \code{60} and the golden section is applied to refine the minimum search. See also \code{\link[rebmix]{kseq}} for sequence of bins or nearest neighbours generation. The default value is \code{"auto"}.
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
a list of lists each containing \eqn{c} parametric family types \code{pdfl}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"Gumbel"}, \code{"binomial"}, \code{"Poisson"}, \code{"Dirac"}, \code{"uniform"} or circular \code{"vonMises"} defined for \eqn{0 \leq y_{i} \leq 2 \pi}.
Component parameters \code{theta1.l} follow the parametric family types. One of \eqn{\mu_{il}} for normal, lognormal, Gumbel and von Mises distributions, \eqn{\theta_{il}} for Weibull, gamma, binomial, Poisson and Dirac distributions and \eqn{a} for uniform distribution.
Component parameters \code{theta2.l} follow \code{theta1.l}. One of \eqn{\sigma_{il}} for normal, lognormal and Gumbel distributions, \eqn{\beta_{il}} for Weibull and gamma distributions, \eqn{p_{il}} for binomial distribution, \eqn{\kappa_{il}} for von Mises distribution and \eqn{b} for uniform distribution.
Component parameters \code{theta3.l} follow \code{theta2.l}. One of \eqn{\xi_{il}} for Gumbel distribution.
}
\item{\code{summary}:}{
a data frame with additional information about dataset, preprocessing, \eqn{c_{\mathrm{max}}}, \eqn{c_{\mathrm{min}}}, information criterion type,
\eqn{a_{\mathrm{r}}}, restraints type, optimal \eqn{c}, optimal \eqn{v} or \eqn{k}, \eqn{K}, \eqn{y_{i0}}, \eqn{y_{i\mathrm{min}}}, \eqn{y_{i\mathrm{max}}}, optimal \eqn{h_{i}},
information criterion \eqn{\mathrm{IC}}, log likelihood \eqn{\mathrm{log}\, L} and degrees of freedom \eqn{M}.
}
\item{\code{summary.EM}:}{
a data frame with additional information about dataset, strategy for the EM algorithm \code{strategy},
variant of the EM algorithm \code{variant}, acceleration type \code{acceleration}, tolerance \code{tolerance}, acceleration multilplier \code{acceleration.multiplier},
maximum allowed number of iterations \code{maximum.iterations}, number of iterations used for obtaining optimal solution \code{opt.iterations.nbr} and total number of iterations of the EM algorithm \code{total.iterations.nbr}.
}
\item{\code{pos}:}{
position in the \code{summary} data frame at which log likelihood \eqn{\mathrm{log}\, L} attains its maximum.
}
\item{\code{opt.c}:}{
a list of vectors containing numbers of components for optimal \eqn{v} for the histogram and the kernel density estimation or for optimal number of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
\item{\code{opt.IC}:}{
a list of vectors containing information criteria for optimal \eqn{v} for the histogram and the kernel density estimation or for optimal number of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
\item{\code{opt.logL}:}{
a list of vectors containing log likelihoods for optimal \eqn{v} for the histogram and the kernel density estimation or for optimal number of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
\item{\code{opt.D}:}{
a list of vectors containing totals of positive relative deviations for optimal \eqn{v} for the histogram and the kernel density estimation or for optimal number of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
\item{\code{all.K}:}{
a list of vectors containing all processed numbers of bins \eqn{v} for the histogram and the kernel density estimation or all processed numbers of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
\item{\code{all.IC}:}{
a list of vectors containing information criteria for all processed numbers of bins \eqn{v} for the histogram and the kernel density estimation or for all processed numbers of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour.
}
}
}
\author{Marko Nagode}
\keyword{classes}
