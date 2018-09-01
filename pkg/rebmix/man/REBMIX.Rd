\name{REBMIX-methods}
\docType{methods}
\alias{REBMIX}
\alias{REBMIX-methods}
\alias{REBMIX,REBMIX-method}
\alias{REBMIX,REBMVNORM-method}
\alias{show,REBMIX-method}
\alias{show,REBMVNORM-method}
\alias{summary,REBMIX-method}
\alias{summary,REBMVNORM-method}
\alias{g.Preprocessing}
\alias{g.cmax}
\alias{g.Criterion}
\alias{g.pdf}
\alias{g.theta1}
\alias{g.theta2}
\alias{g.K}
\alias{g.y0}
\alias{g.ar}
\alias{g.Restraints}
\alias{g.summary}
\alias{g.pos}
\alias{g.opt.c}
\alias{g.opt.IC}
\alias{g.opt.logL}
\alias{g.opt.D}
\alias{g.all.K}
\alias{g.all.IC}
\alias{g.Dataset,REBMIX-method}
\alias{g.Preprocessing,REBMIX-method}
\alias{g.cmax,REBMIX-method}
\alias{g.Criterion,REBMIX-method}
\alias{g.Variables,REBMIX-method}
\alias{g.pdf,REBMIX-method}
\alias{g.theta1,REBMIX-method}
\alias{g.theta2,REBMIX-method}
\alias{g.K,REBMIX-method}
\alias{g.y0,REBMIX-method}
\alias{g.ymin,REBMIX-method}
\alias{g.ymax,REBMIX-method}
\alias{g.ar,REBMIX-method}
\alias{g.Restraints,REBMIX-method}
\alias{g.w,REBMIX-method}
\alias{g.Theta,REBMIX-method}
\alias{g.summary,REBMIX-method}
\alias{g.pos,REBMIX-method}
\alias{g.opt.c,REBMIX-method}
\alias{g.opt.IC,REBMIX-method}
\alias{g.opt.logL,REBMIX-method}
\alias{g.opt.D,REBMIX-method}
\alias{g.all.K,REBMIX-method}
\alias{g.all.IC,REBMIX-method}
\alias{g.Dataset,REBMVNORM-method}
\alias{g.Preprocessing,REBMVNORM-method}
\alias{g.cmax,REBMVNORM-method}
\alias{g.Criterion,REBMVNORM-method}
\alias{g.Variables,REBMVNORM-method}
\alias{g.pdf,REBMVNORM-method}
\alias{g.theta1,REBMVNORM-method}
\alias{g.theta2,REBMVNORM-method}
\alias{g.K,REBMVNORM-method}
\alias{g.y0,REBMVNORM-method}
\alias{g.ymin,REBMVNORM-method}
\alias{g.ymax,REBMVNORM-method}
\alias{g.ar,REBMVNORM-method}
\alias{g.Restraints,REBMVNORM-method}
\alias{g.w,REBMVNORM-method}
\alias{g.Theta,REBMVNORM-method}
\alias{g.summary,REBMVNORM-method}
\alias{g.pos,REBMVNORM-method}
\alias{g.opt.c,REBMVNORM-method}
\alias{g.opt.IC,REBMVNORM-method}
\alias{g.opt.logL,REBMVNORM-method}
\alias{g.opt.D,REBMVNORM-method}
\alias{g.all.K,REBMVNORM-method}
\alias{g.all.IC,REBMVNORM-method}
\title{
REBMIX Algorithm for Univariate or Multivariate Finite Mixture Estimation
}
\description{
Returns as default the REBMIX algorithm output for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson, Dirac or von Mises component densities. If \code{model} equals \code{"REBMVNORM"} output for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices is returned.
}
\usage{
\S4method{REBMIX}{REBMIX}(model = "REBMIX", Dataset = list(), Preprocessing = character(), 
       cmax = 15, Criterion = "AIC", pdf = character(), theta1 = numeric(), 
       theta2 = numeric(), K = "auto", y0 = numeric(), ymin = numeric(), 
       ymax = numeric(), ar = 0.1, Restraints = "loose", ...)
## ... and for other signatures
\S4method{summary}{REBMIX}(object, ...)   
## ... and for other signatures 
}
\arguments{
\item{model}{
see Methods section below.
}
\item{Dataset}{
a list of data frames of size \eqn{n \times d} containing \emph{d}-dimensional datasets. Each of the \eqn{d} columns
represents one random variable. Numbers of observations \eqn{n} equal the number of rows in the datasets.
}
\item{Preprocessing}{
a character vector giving the preprocessing types. One of \code{"histogram"}, \code{"kernel density estimation"} or \code{"k-nearest neighbour"}.
}
\item{cmax}{
maximum number of components \eqn{c_{\mathrm{max}} > 0}. The default value is \code{15}.
}
\item{Criterion}{
a character vector giving the information criterion types. One of default Akaike \code{"AIC"}, \code{"AIC3"}, \code{"AIC4"} or \code{"AICc"}, 
Bayesian \code{"BIC"}, consistent Akaike \code{"CAIC"}, Hannan-Quinn \code{"HQC"}, minimum description length \code{"MDL2"} or \code{"MDL5"}, 
approximate weight of evidence \code{"AWE"}, classification likelihood \code{"CLC"}, 
integrated classification likelihood \code{"ICL"} or \code{"ICL-BIC"}, partition coefficient \code{"PC"}, 
total of positive relative deviations \code{"D"} or sum of squares error \code{"SSE"}.
}
\item{pdf}{
a character vector of length \eqn{d} containing continuous or discrete parametric family types. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"}, \code{"Dirac"} or \code{"vonMises"}.
}
\item{theta1}{
a vector of length \eqn{d} containing initial component parameters. One of \eqn{n_{il} = \textrm{number of categories} - 1} for \code{"binomial"} distribution or \code{"NA"} otherwise. 
}
\item{theta2}{
a vector of length \eqn{d} containing initial component parameters. Currently not used.
}
\item{K}{
a character or a vector or a list of vectors containing numbers of bins \eqn{v} for the histogram and the kernel density estimation or numbers of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour. There is no genuine rule to identify \eqn{v} or \eqn{k}. Consequently,
the REBMIX algorithm identifies them from the set \code{K} of input values by
minimizing the information criterion. The Sturges rule \eqn{v = 1 + \mathrm{log_{2}}(n)}, \eqn{\mathrm{Log}_{10}} rule \eqn{v = 10 \mathrm{log_{10}}(n)} or RootN 
rule \eqn{v = 2 \sqrt{n}} can be applied to estimate the limiting numbers of bins
or the rule of thumb \eqn{k = \sqrt{n}} to guess the intermediate number of nearest neighbours. If, e.g., \code{K = c(10, 20, 40, 60)} and minimum \code{IC} coincides, e.g., \code{40}, brackets are set to \code{20} and \code{60} and the golden section is applied to refine the minimum search. See also \code{\link[rebmix]{kseq}} for sequence of bins or nearest neighbours generation. The default value is \code{"auto"}.
}
\item{y0}{
a vector of length \eqn{d} containing origins. The default value is \code{numeric()}.
}
\item{ymin}{
a vector of length \eqn{d} containing minimum observations. The default value is \code{numeric()}.
}
\item{ymax}{
a vector of length \eqn{d} containing maximum observations. The default value is \code{numeric()}.
}
\item{ar}{
acceleration rate \eqn{0 < a_{\mathrm{r}} \leq 1}. The default value is \code{0.1} and in most cases does not have to be altered.
}
\item{Restraints}{
a character giving the restraints type. One of \code{"rigid"} or default \code{"loose"}.
The rigid restraints are obsolete and applicable for well separated components only. 
}
\item{object}{
see Methods section below.
}
\item{\dots}{
currently not used; additional arguments affecting the summary produced.
}
}
\value{
Returns an object of class \code{REBMIX} or \code{REBMVNORM}.
}
\section{Methods}{
\describe{
\item{\code{signature(model = "REBMIX")}}{a character giving the default class name \code{"REBMIX"} for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson, Dirac or von Mises component densities.}
\item{\code{signature(model = "REBMVNORM")}}{a character giving the class name \code{"REBMVNORM"} for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices.}
\item{\code{signature(object = "REBMIX")}}{an object of class \code{REBMIX}.}
\item{\code{signature(object = "REBMVNORM")}}{an object of class \code{REBMVNORM}.}
}
}
\author{Marko Nagode}
\references{
H. A. Sturges. The choice of a class interval. Journal of American Statistical Association, 21(153):
65-66, 1926. \url{http://www.jstor.org/stable/2965501}.\cr\cr
P. F. Velleman. Interactive computing for exploratory data analysis I: display algorithms. Proceedings of the Statistical Computing Section,
American Statistical Association, 1976.\cr\cr 
W. J. Dixon and R. A. Kronmal. The Choice of origin and scale for graphs. Journal of the ACM, 12(2):
259-261, 1965. \url{http://dx.doi.org/10.1145/321264.321277}.\cr\cr
M. Nagode and M. Fajdiga. A general multi-modal probability density function suitable for the
rainflow ranges of stationary random processes. International Journal of Fatigue, 20(3):211-223,
1998. \url{http://dx.doi.org/10.1016/S0142-1123(97)00106-0}.\cr\cr
M. Nagode and M. Fajdiga. An improved algorithm for parameter estimation suitable for mixed
weibull distributions. International Journal of Fatigue, 22(1):75-80, 2000. \url{http://dx.doi.org/10.1016/S0142-1123(99)00112-7}.\cr\cr
M. Nagode, J. Klemenc, and M. Fajdiga. Parametric modelling and scatter prediction of rainflow
matrices. International Journal of Fatigue, 23(6):525-532, 2001. \url{http://dx.doi.org/10.1016/S0142-1123(01)00007-X}.\cr\cr
M. Nagode and M. Fajdiga. An alternative perspective on the mixture estimation problem. Reliability
Engineering & System Safety, 91(4):388-397, 2006. \url{http://dx.doi.org/10.1016/j.ress.2005.02.005}.\cr\cr
M. Nagode and M. Fajdiga. The rebmix algorithm for the univariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(5):876-892, 2011a. \url{http://dx.doi.org/10.1080/03610920903480890}.\cr\cr
M. Nagode and M. Fajdiga. The rebmix algorithm for the multivariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(11):2022-2034, 2011b. \url{http://dx.doi.org/10.1080/03610921003725788}.\cr\cr
M. Nagode. Finite mixture modeling via REBMIX.
Journal of Algorithms and Optimization, 3(2):14-28, 2015. \url{http://dx.doi.org/10.5963/JAO0302001}.
}
\examples{
# Generate and plot univariate normal dataset.

n <- c(998, 263, 1086, 487)

Theta <- list(pdf1 = "n",
  theta1.1 = c(688, 265, 30, 934),
  theta2.1 = c(72, 54, 34, 28))

normal <- RNGMIX(Dataset.name = "complex1",
  rseed = -1,
  n = n,
  Theta = Theta)
  
normal 

g.Dataset(normal, 1)[1:20,]  

# Estimate number of components, component weights and component parameters. 

normalest <- REBMIX(Dataset = g.Dataset(normal), 
  Preprocessing = "h", 
  cmax = 8, 
  Criterion = "BIC", 
  pdf = "n")
                 
normalest

BIC(normalest)

logL(normalest)

# Plot finite mixture.

plot(normalest, nrow = 2, what = c("den", "dis"), npts = 1000)
}
\keyword{parameter estimation}
