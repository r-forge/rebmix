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
\title{
REBMIX Algorithm for Univariate or Multivariate Finite Mixture Estimation
}
\description{
Returns as default the REBMIX algorithm output for mixtures of conditionally independent normal, lognormal, Weibull, gamma, Gumbel, binomial, Poisson, Dirac, uniform or von Mises component densities. If \code{model} equals \code{"REBMVNORM"} output for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices is returned.
}
\usage{
\S4method{REBMIX}{REBMIX}(model = "REBMIX", Dataset = list(), Preprocessing = character(),
       cmax = 15, cmin = 1, Criterion = "AIC", pdf = character(), 
       theta1 = numeric(), theta2 = numeric(), theta3 = numeric(), K = "auto", 
       ymin = numeric(), ymax = numeric(), ar = 0.1, 
       Restraints = "loose", EMcontrol = NULL, ...)
## ... and for other signatures
\S4method{summary}{REBMIX}(object, ...)
## ... and for other signatures
}
\arguments{
\item{model}{
see Methods section below.
}
\item{Dataset}{
a list of length \eqn{n_{\mathrm{D}}} of data frames or objects of class \code{Histogram}. 
Data frames should have size \eqn{n \times d} containing \emph{d}-dimensional datasets. Each of the \eqn{d}
columns represents one random variable. Numbers of observations \eqn{n} equal the number of rows in the datasets.
}
\item{Preprocessing}{
a character giving the preprocessing type. One of \code{"histogram"}, \cr
\code{"kernel density estimation"} or \code{"k-nearest neighbour"}.
}
\item{cmax}{
maximum number of components \eqn{c_{\mathrm{max}} > 0}. The default value is \code{15}.
}
\item{cmin}{
minimum number of components \eqn{c_{\mathrm{min}} > 0}. The default value is \code{1}.
}
\item{Criterion}{
a character giving the information criterion type. One of default Akaike \code{"AIC"}, \code{"AIC3"}, \code{"AIC4"} or \code{"AICc"},
Bayesian \code{"BIC"}, consistent Akaike \code{"CAIC"}, Hannan-Quinn \code{"HQC"}, minimum description length \code{"MDL2"} or \code{"MDL5"},
approximate weight of evidence \code{"AWE"}, classification likelihood \code{"CLC"},
integrated classification likelihood \code{"ICL"} or \code{"ICL-BIC"}, partition coefficient \code{"PC"},
total of positive relative deviations \code{"D"} or sum of squares error \code{"SSE"}.
}
\item{pdf}{
a character vector of length \eqn{d} containing continuous or discrete parametric family types. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"Gumbel"}, \code{"binomial"}, \code{"Poisson"}, \code{"Dirac"}, \code{"uniform"} or \code{"vonMises"}.
}
\item{theta1}{
a vector of length \eqn{d} containing initial component parameters. One of \eqn{n_{il} = \textrm{number of categories} - 1} for \code{"binomial"} distribution.
}
\item{theta2}{
a vector of length \eqn{d} containing initial component parameters. Currently not used.
}
\item{theta3}{
a vector of length \eqn{d} containing initial component parameters. One of \eqn{\xi_{il} \in \{-1, \textrm{NA}, 1\}} for \code{"Gumbel"} distribution.
}
\item{K}{
a character or a vector or a matrix of size \eqn{n_{\mathrm{D}} \times d} containing numbers of bins \eqn{v} or \eqn{v_{1}, \ldots, v_{d}} for the histogram and the kernel density estimation or numbers of nearest
neighbours \eqn{k} for the \emph{k}-nearest neighbour. There is no genuine rule to identify \eqn{v} or \eqn{k}. Consequently,
the REBMIX algorithm identifies them from the set \code{K} of input values by
minimizing the information criterion. The Sturges rule \eqn{v = 1 + \mathrm{log_{2}}(n)}, \eqn{\mathrm{Log}_{10}} rule \eqn{v = 10 \mathrm{log_{10}}(n)} or RootN
rule \eqn{v = 2 \sqrt{n}} can be applied to estimate the limiting numbers of bins
or the rule of thumb \eqn{k = \sqrt{n}} to guess the intermediate number of nearest neighbours. If, e.g., \code{K = c(10, 20, 40, 60)} and minimum \code{IC} coincides, e.g., \code{40}, brackets are set to \code{20} and \code{60} and the golden section is applied to refine the minimum search. 
If, e.g., \code{K = matrix(c(10, 15, 18, 5, 7, 9), byrow = TRUE, ncol = 3)} than \eqn{d = 3} and the list \code{Dataset} contains \eqn{n_{\mathrm{D}} = 2} frames. Hence, different numbers of bins can be assigned to \eqn{y_{1}, \ldots, y_{d}}.
See also \code{\link[rebmix]{kseq}} for sequence of bins or nearest neighbours generation. The default value is \code{"auto"}.
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
\item{EMcontrol}{
an object of class \code{EM.Control}.
}
\item{object}{
see Methods section below.
}
\item{\dots}{
currently not used.
}
}
\value{
Returns an object of class \code{REBMIX} or \code{REBMVNORM}.
}
\section{Methods}{
\describe{
\item{\code{signature(model = "REBMIX")}}{a character giving the default class name \code{"REBMIX"} for mixtures of conditionally independent normal, lognormal, Weibull, gamma, Gumbel, binomial, Poisson, Dirac, uniform or von Mises component densities.}
\item{\code{signature(model = "REBMVNORM")}}{a character giving the class name \code{"REBMVNORM"} for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices.}
\item{\code{signature(object = "REBMIX")}}{an object of class \code{REBMIX}.}
\item{\code{signature(object = "REBMVNORM")}}{an object of class \code{REBMVNORM}.}
}
}
\author{Marko Nagode}
\references{
H. A. Sturges. The choice of a class interval. Journal of American Statistical Association, 21(153):
65-66, 1926. \url{https://www.jstor.org/stable/2965501}.\cr\cr
P. F. Velleman. Interactive computing for exploratory data analysis I: display algorithms. Proceedings of the Statistical Computing Section,
American Statistical Association, 1976.\cr\cr
W. J. Dixon and R. A. Kronmal. The Choice of origin and scale for graphs. Journal of the ACM, 12(2):
259-261, 1965. \doi{10.1145/321264.321277}.\cr\cr
M. Nagode and M. Fajdiga. A general multi-modal probability density function suitable for the
rainflow ranges of stationary random processes. International Journal of Fatigue, 20(3):211-223,
1998. \doi{10.1016/S0142-1123(97)00106-0}.\cr\cr
M. Nagode and M. Fajdiga. An improved algorithm for parameter estimation suitable for mixed
weibull distributions. International Journal of Fatigue, 22(1):75-80, 2000. \doi{10.1016/S0142-1123(99)00112-7}.\cr\cr
M. Nagode, J. Klemenc and M. Fajdiga. Parametric modelling and scatter prediction of rainflow
matrices. International Journal of Fatigue, 23(6):525-532, 2001. \doi{10.1016/S0142-1123(01)00007-X}.\cr\cr
M. Nagode and M. Fajdiga. An alternative perspective on the mixture estimation problem. Reliability
Engineering & System Safety, 91(4):388-397, 2006. \doi{10.1016/j.ress.2005.02.005}.\cr\cr
M. Nagode and M. Fajdiga. The rebmix algorithm for the univariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(5):876-892, 2011a. \doi{10.1080/03610920903480890}.\cr\cr
M. Nagode and M. Fajdiga. The rebmix algorithm for the multivariate finite mixture estimation.
Communications in Statistics - Theory and Methods, 40(11):2022-2034, 2011b. \doi{10.1080/03610921003725788}.\cr\cr
M. Nagode. Finite mixture modeling via REBMIX.
Journal of Algorithms and Optimization, 3(2):14-28, 2015. \url{https://repozitorij.uni-lj.si/Dokument.php?id=127674&lang=eng}.\cr\cr
B. Panic, J. Klemenc, M. Nagode. Improved initialization of the EM algorithm for mixture model parameter estimation.
Mathematics, 8(3):373, 2020.
\doi{10.3390/math8030373}.
}
\examples{
# Generate and plot univariate normal dataset.

n <- c(998, 263, 1086, 487)

Theta <- new("RNGMIX.Theta", c = 4, pdf = "normal")

a.theta1(Theta) <- c(688, 265, 30, 934)
a.theta2(Theta) <- c(72, 54, 34, 28)

normal <- RNGMIX(Dataset.name = "complex1",
  rseed = -1,
  n = n,
  Theta = a.Theta(Theta))

normal

a.Dataset(normal, 1)[1:20,]

# Estimate number of components, component weights and component parameters.

normalest <- REBMIX(Dataset = a.Dataset(normal),
  Preprocessing = "h",
  cmax = 8,
  Criterion = "BIC",
  pdf = "n")

normalest

BIC(normalest)

logL(normalest)

# Plot finite mixture.

plot(normalest, nrow = 2, what = c("pdf", "marginal cdf"), npts = 1000)

# EM algorithm utilization

# Load iris data.

data(iris)

Dataset <- list(data.frame(iris[, c(1:4)]))

# Create EM.Control object.

EM <- new("EM.Control", 
  strategy = "exhaustive", 
  variant = "EM", 
  acceleration = "fixed", 
  tolerance = 1e-4, 
  acceleration.multiplier = 1.0, 
  maximum.iterations = 1000)

# Mixture parameter estimation using REBMIX and EM algorithm.

irisest <- REBMIX(model = "REBMVNORM", 
  Dataset = Dataset, 
  Preprocessing = "histogram",
  cmax = 10,
  Criterion = "BIC",
  EMcontrol = EM)

irisest

# Print total number of EM iterations used in Ehxaustive strategy from summary.EM slot.

a.summary.EM(irisest, col.name = "total.iterations.nbr", pos = 1)
}
\keyword{parameter estimation}
