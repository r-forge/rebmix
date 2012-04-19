\name{REBMIX}
\alias{REBMIX}
\title{
REBMIX Algorithm for Univariate or Multivariate Finite Mixture Estimation
}
\description{
Returns the REBMIX algorithm output for mixtures of conditionally independent normal, lognormal, Weibull, binomial, Poisson or Dirac component densities.
}
\usage{
REBMIX(Dataset = NULL, Preprocessing = NULL, D = 0.025, cmax = 15, 
       Criterion = "AIC", Variables = NULL, pdf = NULL,
       Theta1 = NULL, Theta2 = NULL, K = NULL, ymin = NULL, 
       ymax = NULL, b = 1.0, ar = 0.1, Restraints = "loose")
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
decrease. If \eqn{D = 0} and \eqn{b = 0} the mixture tends to \eqn{c = c_{\mathrm{max}}} components.
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
a character vector of length \eqn{d} containing continuous or discrete parametric family types. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
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
\item{b}{
minimum weight multiplier \eqn{0 \leq b \leq 1} influences the number of components \eqn{w_{\mathrm{min}} = 2 D_{\mathrm{min}} ((l - 1)b + 1)}. The
default value is \code{1.0}.
}
\item{ar}{
acceleration rate \eqn{0 < a_{\mathrm{r}} \leq 1}. The default value is \code{0.1} and in most cases does not have to be altered.
}
\item{Restraints}{
a character string giving the restraints type. One of \code{"rigid"} or default \code{"loose"}.
The rigid restraints are obsolete and applicable for well separated components only. 
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
a list of data frames each containing  \eqn{c} parametric family types \code{pdfi}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
Component parameters \code{theta1.i} follow the parametric family types. One of \eqn{\mu_{il}} for normal and lognormal distributions and \eqn{\theta_{il}} for Weibull, binomial, Poisson and Dirac distributions.
Component parameters \code{theta2.i} follow \code{theta1.i}. One of \eqn{\sigma_{il}} for normal and lognormal distributions, \eqn{\beta_{il}} for Weibull distribution and \eqn{p_{il}} for binomial distribution.
}
\item{Variables}{
a character vector containing types of variables. One of \code{"continuous"} or \code{"discrete"}.
}
\item{pdf}{
a character vector containing continuous or discrete parametric family types. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"binomial"}, \code{"Poisson"} or \code{"Dirac"}.
}
\item{Theta1}{
a vector containing initial component parameters. One of \eqn{n_{il} = \textrm{Number of categories} - 1} for \code{"binomial"} distribution or \code{"NA"} otherwise. 
}
\item{Theta2}{
a vector containing initial component parameters. The value is \code{NULL}.
}
\item{summary}{
a data frame with additional information about dataset, preprocessing, \eqn{D}, \eqn{c_{\mathrm{max}}}, information criterion type, 
\eqn{b}, \eqn{a_{\mathrm{r}}}, restraints type, optimal \eqn{c}, optimal \eqn{v} or \eqn{k}, \eqn{y_{i0}}, optimal \eqn{h_{i}},
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
Sturges HA (1926). The choice of a class interval. \emph{Journal of American Statistical Association},
21, 65-66.\cr\cr
Nagode M, Fajdiga M (1998). A General Multi-Modal Probability Density Function Suitable for 
the Rainflow Ranges of Stationary Random Processes. \emph{International Journal of Fatigue}, 20, 211-223.\cr\cr
Nagode M, Fajdiga M (2000). An Improved Algorithm for Parameter Estimation Suitable
for Mixed Weibull Distributions. \emph{International Journal of Fatigue}, 22, 75-80.\cr\cr
Nagode M, Klemenc J, Fajdiga M (2001). Parametric Modelling and Scatter Prediction of
Rainflow Matrices. \emph{International Journal of Fatigue}, 23, 525-532.\cr\cr
Nagode M, Fajdiga M (2006). An Alternative Perspective on the Mixture Estimation Problem. 
\emph{Reliability Engineering & System Safety}, 91, 388-397.\cr\cr
Nagode M, Fajdiga M (2011). The REBMIX Algorithm for the Univariate Finite Mixture
Estimation. \emph{Communications in Statistics - Theory and Methods}, 40(5), 876-892.\cr\cr
Nagode M, Fajdiga M (2011). The REBMIX Algorithm for the Multivariate Finite Mixture
Estimation. \emph{Communications in Statistics - Theory and Methods}. 40(11), 2022-2034.
}
\examples{
## Generate complex dataset.

n <- c(99, 26, 109, 49, 21, 108, 23, 
  78, 84, 46, 77, 2, 81, 109, 86)

Theta <- rbind(pdf = "n",
  theta1 = c(688.4, 265.1, 30.8, 934, 561.6, 854.9, 883.7, 
  758.3, 189.3, 919.3, 98, 143, 202.5, 628, 977),
  theta2 = c(12.4, 14.6, 14.8, 8.4, 11.7, 9.2, 6.3, 10.2,
  9.5, 8.1, 14.7, 11.7, 7.4, 10.1, 14.6))

complex <- RNGMIX(Dataset = "complex",
  rseed = -1,
  n = n,
  Theta = Theta)

complex 

complex$Dataset$complex

## Estimate number of components, component weights and component parameters. 

K <- c(as.integer(1 + log2(sum(n))), ## Minimum v follows the Sturges rule.
  as.integer(2 * sum(n)^0.5)) ## Maximum v follows the RootN rule.

complexest <- REBMIX(Dataset = complex$Dataset, 
  Preprocessing = "histogram", 
  D = 0.0025, 
  cmax = 30  , 
  Criterion = "BIC", 
  Variables = "continuous",
  pdf = "normal", 
  K = K, 
  b = 0.0)

complexest
}
\keyword{parameter estimation}
