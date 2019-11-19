\name{RNGMIX-methods}
\docType{methods}
\alias{RNGMIX}
\alias{RNGMIX-methods}
\alias{RNGMIX,RNGMIX-method}
\alias{RNGMIX,RNGMVNORM-method}
\alias{show,RNGMIX-method}
\alias{show,RNGMVNORM-method}
\title{
Random Univariate or Multivariate Finite Mixture Generation
}
\description{
Returns as default the RNGMIX univariate or multivariate random datasets for mixtures of conditionally independent normal, lognormal, Weibull, gamma, Gumbel, binomial, Poisson, Dirac, uniform or von Mises component densities.
If \code{model} equals \code{"RNGMVNORM"} multivariate random datasets for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices are returned.
}
\usage{
\S4method{RNGMIX}{RNGMIX}(model = "RNGMIX", Dataset.name = character(),
       rseed = -1, n = numeric(), Theta = list(), ...)
## ... and for other signatures
}
\arguments{
\item{model}{
see Methods section below.
}
\item{Dataset.name}{
a character vector containing list names of data frames of size \eqn{n \times d} that \emph{d}-dimensional datasets are written in.
}
\item{rseed}{
set the random seed to any negative integer value to initialize the sequence. The first file in \code{Dataset.name} corresponds to it.
For each next file the random seed is decremented \eqn{r_{\mathrm{seed}} = r_{\mathrm{seed}} - 1}. The default value is \code{-1}.
}
\item{n}{
a vector containing numbers of observations in classes \eqn{n_{l}}, where number of observations \eqn{n = \sum_{l = 1}^{c} n_{l}}.
}
\item{Theta}{
a list containing \eqn{c} parametric family types \code{pdfl}. One of \code{"normal"}, \code{"lognormal"}, \code{"Weibull"}, \code{"gamma"}, \code{"Gumbel"}, \code{"binomial"}, \code{"Poisson"}, \code{"Dirac"}, \code{"uniform"} or circular \code{"vonMises"} defined for \eqn{0 \leq y_{i} \leq 2 \pi}.
Component parameters \code{theta1.l} follow the parametric family types. One of \eqn{\mu_{il}} for normal, lognormal, Gumbel and von Mises distributions and \eqn{\theta_{il}} for Weibull, gamma, binomial, Poisson and Dirac distributions.
Component parameters \code{theta2.l} follow \code{theta1.l}. One of \eqn{\sigma_{il}} for normal and lognormal distributions, \eqn{\beta_{il}} for Weibull, gamma and Gumbel distributions, \eqn{p_{il}} for binomial distribution, \eqn{\kappa_{il}} for von Mises distribution and \code{"NA"} otherwise.
}
\item{\dots}{
currently not used.
}
}
\details{
RNGMIX is based on the "Minimal" random number generator \code{ran1} of Park and Miller with the Bays-Durham shuffle and added safeguards that returns a uniform random deviate between 0.0 and 1.0
(exclusive of the endpoint values).
}
\value{
Returns an object of class \code{RNGMIX} or \code{RNGMVNORM}.
}
\section{Methods}{
\describe{
\item{\code{signature(model = "RNGMIX")}}{a character giving the default class name \code{"RNGMIX"} for mixtures of conditionally independent normal, lognormal, Weibull, gamma, Gumbel, binomial, Poisson, Dirac or von Mises component densities.}
\item{\code{signature(model = "RNGMVNORM")}}{a character giving the class name \code{"RNGMVNORM"} for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices.}
}
}
\author{Marko Nagode}
\references{
W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery. Numerical Recipes in C: The Art
of Scientific Computing. Cambridge University Press, Cambridge, 1992.
}
\examples{
devAskNewPage(ask = TRUE)

# Generate and print multivariate normal datasets with diagonal
# variance-covariance matrices.

n <- c(75, 100, 125, 150, 175)

Theta <- new("RNGMIX.Theta", c = 5, pdf = rep("normal", 4))

a.theta1(Theta, 1) <- c(10, 12, 10, 12)
a.theta1(Theta, 2) <- c(8.5, 10.5, 8.5, 10.5)
a.theta1(Theta, 3) <- c(12, 14, 12, 14)
a.theta1(Theta, 4) <- c(13, 15, 7, 9)
a.theta1(Theta, 5) <- c(7, 9, 13, 15)
a.theta2(Theta, 1) <- c(1, 1, 1, 1)
a.theta2(Theta, 2) <- c(1, 1, 1, 1)
a.theta2(Theta, 3) <- c(1, 1, 1, 1)
a.theta2(Theta, 4) <- c(2, 2, 2, 2)
a.theta2(Theta, 5) <- c(3, 3, 3, 3)

simulated <- RNGMIX(Dataset.name = paste("simulated_", 1:25, sep = ""),
  rseed = -1,
  n = n,
  Theta = a.Theta(Theta))

simulated

plot(simulated, pos = 22, nrow = 2, ncol = 3)

# Generate and print multivariate normal datasets with unrestricted
# variance-covariance matrices.

n <- c(200, 50, 50)

Theta <- new("RNGMVNORM.Theta", c = 3, d = 3)

a.theta1(Theta, 1) <- c(0, 0, 0)
a.theta1(Theta, 2) <- c(-6, 3, 6)
a.theta1(Theta, 3) <- c(6, 6, 4)
a.theta2(Theta, 1) <- c(9, 0, 0, 0, 4, 0, 0, 0, 1)
a.theta2(Theta, 2) <- c(4, -3.2, -0.2, -3.2, 4, 0, -0.2, 0, 1)
a.theta2(Theta, 3) <- c(4, 3.2, 2.8, 3.2, 4, 2.4, 2.8, 2.4, 2)

simulated <- RNGMIX(model = "RNGMVNORM",
  Dataset.name = paste("simulated_", 1:2, sep = ""),
  rseed = -1,
  n = n,
  Theta = a.Theta(Theta))

simulated

plot(simulated, pos = 2, nrow = 3, ncol = 1)

# Generate and print multivariate mixed continuous-discrete datasets.

n <- c(400, 100, 500)

Theta <- new("RNGMIX.Theta", c = 3, pdf = c("lognormal", "Poisson", "binomial", "Weibull"))

a.theta1(Theta, 1) <- c(1, 2, 10, 2)
a.theta1(Theta, 2) <- c(3.5, 10, 10, 10)
a.theta1(Theta, 3) <- c(2.5, 15, 10, 25)
a.theta2(Theta, 1) <- c(0.3, NA, 0.9, 3)
a.theta2(Theta, 2) <- c(0.2, NA, 0.1, 7)
a.theta2(Theta, 3) <- c(0.4, NA, 0.7, 20)

simulated <- RNGMIX(Dataset.name = paste("simulated_", 1:5, sep = ""),
  rseed = -1,
  n = n,
  Theta = a.Theta(Theta))

simulated

plot(simulated, pos = 4, nrow = 2, ncol = 3)

# Generate and print univariate mixed Weibull dataset.

n <- c(75, 100, 125, 150, 175)

Theta <- new("RNGMIX.Theta", c = 5, pdf = "Weibull")

a.theta1(Theta) <- c(12, 10, 14, 15, 9)
a.theta2(Theta) <- c(2, 4.1, 3.2, 7.1, 5.3)

simulated <- RNGMIX(Dataset.name = "simulated",
  rseed = -1,
  n = n,
  Theta = a.Theta(Theta))

simulated

plot(simulated, pos = 1)

# Generate and print multivariate normal datasets with unrestricted
# variance-covariance matrices.

# Set dimension, dataset size, number of components and seed.

d <- 2; n <- 1000; c <- 10; set.seed(123)

# Component weights are generated.

w <- runif(c, 0.1, 0.9); w <- w / sum(w)

# Set range of means and rang of eigenvalues.

mu <- c(-100, 100); lambda <- c(1, 100)

# Component means and variance-covariance matrices are calculated.

Mu <- list(); Sigma <- list()
for (l in 1:c) {
  Mu[[l]] <- runif(d, mu[1], mu[2])
  Lambda <- diag(runif(d, lambda[1], lambda[2]), nrow = d, ncol = d)
  P <- svd(matrix(runif(d * d, -1, 1), nc = d))$u
  Sigma[[l]] <- P %*% Lambda %*% t(P)
}

# Numbers of observations are calculated and component means and 
# variance-covariance matrices are stored.

n <- round(w * n); Theta <- list()
for (l in 1:c) {
  Theta[[paste0("pdf", l)]] <- rep("normal", d)
  Theta[[paste0("theta1.", l)]] <- Mu[[l]]
  Theta[[paste0("theta2.", l)]] <- as.vector(Sigma[[l]])
}

# Dataset is generated.

simulated <- RNGMIX(model = "RNGMVNORM", Dataset.name = "mvnorm_1",
  rseed = -1, n = n, Theta = Theta)

plot(simulated)
}
\keyword{random number generation}

