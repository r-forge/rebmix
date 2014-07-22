\name{REBMIX}
\alias{REBMIX}
\title{
REBMIX Algorithm for Univariate or Multivariate Finite Mixture Estimation
}
\description{
Returns the REBMIX algorithm output for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson or Dirac component densities.
}
\usage{
REBMIX(Dataset = NULL, Preprocessing = NULL, D = 0.025, cmax = 15, 
       Criterion = "AIC", Variables = NULL, pdf = NULL,
       Theta1 = NULL, Theta2 = NULL, K = NULL, ymin = NULL, 
       ymax = NULL, ar = 0.1, Restraints = "loose", ...)
}
\arguments{
\item{Dataset}{
a list of data frames of size \eqn{n \times d} containing \emph{d}-dimensional datasets. Each of the \eqn{d} columns
represents one random variable. Number of observations \eqn{n} equals the number of rows in the datasets.
}
\item{Preprocessing}{
a character vector, giving the preprocessing types. One of \code{"histogram"}, \code{"Parzen window"} or \code{"k-nearest neighbour"}.
}
\item{D}{
a total of positive relative deviations standing for the maximum acceptable measure of
distance between predictive and empirical densities. It satisfies the relation \eqn{0 \leq D \leq 1}.
The default value is \code{0.025}. However, if components with a low probability of occurrence are expected, it has to
decrease. If \eqn{D = 0} the mixture tends to \eqn{c = c_{\mathrm{max}}} components.
}
\item{cmax}{
maximum number of components \eqn{c_{\mathrm{max}} > 0}. The default value is \code{15}.
}
\item{Criterion}{
a character vector giving the infromation criterion types. One of default Akaike \code{"AIC"}, \code{"AIC3"}, \code{"AIC4"} or \code{"AICc"}, 
Bayesian \code{"BIC"}, consistent Akaike \code{"CAIC"}, Hannan-Quinn \code{"HQC"}, minimum description length \code{"MDL2"} or \code{"MDL5"}, 
approximate weight of evidence \code{"AWE"}, classification likelihood \code{"CLC"}, 
integrated classification likelihood \code{"ICL"} or \code{"ICL-BIC"}, partition coefficient \code{"PC"}, 
total of positive relative deviations \code{"D"} or sum of squares error \code{"SSE"}.
}
\item{Variables}{
a character vector of length \eqn{d} containing types of variables. One of \code{"continuous"} or \code{"discrete"}.
}
\item{pdf}{
a character vector of length \eqn{d} containing continuous or discrete parametric family types. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
}
\item{Theta1}{
a vector of length \eqn{d} containing initial component parameters. One of \eqn{n_{il} = \textrm{Number of categories} - 1} for \code{"binomial"} distribution or \code{"NA"} otherwise. 
}
\item{Theta2}{
a vector of length \eqn{d} containing initial component parameters. The value is \code{NULL}.
}
\item{K}{
a vector or a list of vectors containing numbers of bins \eqn{v} for the histogram and the Parzen window or numbers of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour. There is no genuine rule to identify \eqn{v} or \eqn{k}. Consequently,
the REBMIX algorithm identifies them from the set \code{K} of input values by
minimizing the information criterion. The Sturges rule \eqn{v = 1 + \mathrm{log_{2}}(n)}, \eqn{\mathrm{Log}_{10}} rule \eqn{v = 10 \mathrm{log_{10}}(n)} or RootN 
rule \eqn{v = 2 \sqrt{n}} can be applied to estimate the limiting numbers of bins
or the rule of thumb \eqn{k = \sqrt{n}} to guess the intermediate number of nearest neighbours.
}
\item{ymin}{
a vector of length \eqn{d} containing minimum observations. The default value is \code{NULL}.
}
\item{ymax}{
a vector of length \eqn{d} containing maximum observations. The default value is \code{NULL}.
}
\item{ar}{
acceleration rate \eqn{0 < a_{\mathrm{r}} \leq 1}. The default value is \code{0.1} and in most cases does not have to be altered.
}
\item{Restraints}{
a character string giving the restraints type. One of \code{"rigid"} or default \code{"loose"}.
The rigid restraints are obsolete and applicable for well separated components only. 
}
\item{...}{
potential further arguments of the method. 
}
}
\value{
\item{Dataset}{
a list of data frames of size \eqn{n \times d} containing \emph{d}-dimensional datasets. Each of the \eqn{d} columns
represents one random variable. Number of observations \eqn{n} equals the number of rows in the datasets.
}
\item{w}{
a list of data frames each containing \eqn{c} component weights \eqn{w_{l}} summing to 1.
}
\item{Theta}{
a list of data frames each containing  \eqn{c} parametric family types \code{pdfi}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
Component parameters \code{theta1.i} follow the parametric family types. One of \eqn{\mu_{il}} for normal and lognormal distributions and \eqn{\theta_{il}} for Weibull, gamma, binomial, Poisson and Dirac distributions.
Component parameters \code{theta2.i} follow \code{theta1.i}. One of \eqn{\sigma_{il}} for normal and lognormal distributions, \eqn{\beta_{il}} for Weibull and gamma distributions and \eqn{p_{il}} for binomial distribution.
}
\item{Variables}{
a character vector containing types of variables. One of \code{"continuous"} or \code{"discrete"}.
}
\item{pdf}{
a character vector containing continuous or discrete parametric family types. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
}
\item{Theta1}{
a vector containing initial component parameters. One of \eqn{n_{il} = \textrm{Number of categories} - 1} for \code{"binomial"} distribution or \code{"NA"} otherwise. 
}
\item{Theta2}{
a vector containing initial component parameters. The value is \code{NULL}.
}
\item{summary}{
a data frame with additional information about dataset, preprocessing, \eqn{D}, \eqn{c_{\mathrm{max}}}, information criterion type, 
\eqn{a_{\mathrm{r}}}, restraints type, optimal \eqn{c}, optimal \eqn{v} or \eqn{k}, \eqn{y_{i0}}, optimal \eqn{h_{i}},
information criterion \eqn{\mathrm{IC}}, log likelihood \eqn{\mathrm{log}\, L} and degrees of freedom \eqn{M}.
}
\item{pos}{
position in the \code{summary} data frame at which log likelihood \eqn{\mathrm{log}\, L} attains its maximum.
}

\item{all.Imax}{
a list of all numbers of iterations.
}
\item{all.c}{
a list of all numbers of components.
}
\item{all.IC}{
a list of all information criteria.
}
\item{all.logL}{
a list of all log lekelihoods.
}
\item{all.D}{
a list of all totals of positive relative deviations.
}
}
\references{
H. A. Sturges. The choice of a class interval. Journal of American Statistical Association, 21(153):
65-66, 1926. URL http://www.jstor.org/stable/2965501.\cr\cr
M. Nagode and M. Fajdiga. A general multi-modal probability density function suitable for the
rainflow ranges of stationary random processes. International Journal of Fatigue, 20(3):211-223,
1998. doi: 10.1016/S0142-1123(97)00106-0.\cr\cr
M. Nagode and M. Fajdiga. An improved algorithm for parameter estimation suitable for mixed
weibull distributions. International Journal of Fatigue, 22(1):75-80, 2000. doi: 10.1016/S0142-
1123(99)00112-7.\cr\cr
M. Nagode, J. Klemenc, and M. Fajdiga. Parametric modelling and scatter prediction of rainflow
matrices. International Journal of Fatigue, 23(6):525-532, 2001. doi: 10.1016/S0142-1123(01)00007-
X.\cr\cr
M. Nagode and M. Fajdiga. An alternative perspective on the mixture estimation problem. Reliability
Engineering & System Safety, 91(4):388-397, 2006. doi: 10.1016/j.ress.2005.02.005.\cr\cr
M. Nagode and M. Fajdiga. The rebmix algorithm for the univariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(5):876-892, 2011a. doi:
10.1080/03610920903480890.\cr\cr
M. Nagode and M. Fajdiga. The rebmix algorithm for the multivariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(11):2022-2034, 2011b. doi:
10.1080/03610921003725788.
}
\examples{
\dontrun{
## Generate the complex 1 dataset.

n <- c(998, 263, 1086, 487, 213, 1076, 232, 
  784, 840, 461, 773, 24, 811, 1091, 861)

Theta <- rbind(pdf = "normal",
  theta1 = c(688.4, 265.1, 30.8, 934, 561.6, 854.9, 883.7, 
  758.3, 189.3, 919.3, 98, 143, 202.5, 628, 977),
  theta2 = c(12.4, 14.6, 14.8, 8.4, 11.7, 9.2, 6.3, 10.2,
  9.5, 8.1, 14.7, 11.7, 7.4, 10.1, 14.6))

complex1 <- RNGMIX(Dataset = "complex1",
  rseed = -1,
  n = n,
  Theta = Theta)
  
complex1

complex1$Dataset[[1]][1:20, ]  

## Estimate number of components, component weights and component parameters. 

v <- c(as.integer(1 + log2(sum(n))), ## Minimum v follows the Sturges rule.
  as.integer(2 * sum(n)^0.5)) ## Maximum v follows the RootN rule.

## Number of classes or nearest neighbours to be processed.

N <- as.integer(log(v[2] / (v[1] + 1)) / log(1 + 1 / v[1]))

K <- c(v[1], as.integer((v[1] + 1) * (1 + 1 / v[1])^(0:N)))

complex1est <- REBMIX(Dataset = complex1$Dataset, 
  Preprocessing = "histogram", 
  D = 0.0025, 
  cmax = 30, 
  Criterion = "BIC", 
  Variables = "continuous",
  pdf = "normal", 
  K = K)
                 
complex1est

## Plot the finite mixture.

plot(complex1est, npts = 1000)
}
}
\keyword{parameter estimation}
