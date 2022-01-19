\name{EMMIX-methods}
\docType{methods}
\alias{EMMIX}
\alias{EMMIX-methods}
\alias{EMMIX,REBMIX-method}
\alias{EMMIX,REBMVNORM-method}
\title{
EM Algorithm for Univariate or Multivariate Finite Mixture Estimation
}
\description{
Returns as default the EM algorithm output for mixtures of conditionally independent normal, lognormal, Weibull, gamma, 
Gumbel, binomial, Poisson, Dirac or von Mises component densities. If \code{model} equals \code{"REBMVNORM"} output 
for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices is returned.
}
\usage{
\S4method{EMMIX}{REBMIX}(model = "REBMIX", Dataset = list(), Theta = list(), ...)
## ... and for other signatures
}
\arguments{
\item{model}{
see Methods section below.
}
\item{Dataset}{
a list of length \eqn{n_{\mathrm{D}}} of data frames of size \eqn{n \times d} containing \emph{d}-dimensional datasets. Each of the \eqn{d} columns
represents one random variable. Numbers of observations \eqn{n} equal the number of rows in the datasets.
}
\item{Theta}{
a list of objects of \code{"EMMIX.Theta"} or \code{"EMMVNORM.Theta"}.
}
\item{\dots}{
Currently not used.  
}
}
\value{
Returns an object of class \code{REBMIX} or \code{REBMVNORM}.
}
\section{Methods}{
\describe{
\item{\code{signature(model = "REBMIX")}}{a character giving the default class name \code{"REBMIX"} for mixtures of conditionally 
independent normal, lognormal, Weibull, gamma, Gumbel, binomial, Poisson, Dirac or von Mises component densities.}
\item{\code{signature(model = "REBMVNORM")}}{a character giving the class name \code{"REBMVNORM"} for mixtures 
of multivariate normal component densities with unrestricted variance-covariance matrices.}
}
}
\author{Branislav Panic}
\references{
B. Panic, J. Klemenc, M. Nagode. Improved initialization of the EM algorithm for mixture model parameter estimation.
Mathematics, 8(3):373, 2020.
\doi{10.3390/math8030373}.
}
\examples{
\dontrun{
devAskNewPage(ask = TRUE)

# Load faithful dataset.
data(faithful)
# Plot faithfull dataset.
plot(faithful)

# Number of dimensions.

d <- ncol(faithful)

# Obtain 2 component solution with Gaussian mixtures.

c <- 2

# Set some initial parameter guesses.
# Weights.

w <- c(0.5, 0.5)

# Means.

theta1 <- c(2.0, 55.0, 4.5, 80.0)

# Covariances.

theta2 <- c(1, 0, 0, 1, 1, 0, 0, 1)

# Create EMMVNORM.Theta object with new call.

Theta <- new("EMMVNORM.Theta", d = d, c = c)

# Set parameters of Theta.
# Weights.

a.w(Theta) <- w

# Means.

a.theta1.all(Theta) <- theta1

# Covariances

a.theta2.all(Theta) <- theta2

# Run EMMIX method.

model <- EMMIX(model = "REBMVNORM", Dataset = list(faithful), Theta = Theta)

# show

model

# summary

summary(model)

# plot

plot(model)

# customize EM.control object
# strategy argument should be always set to single
# changing other arguments can have effect
# for example remove components with zero weight

EM <- new("EM.Control", strategy = "single", eliminate.zero.components = TRUE)

c <- 4

w <- c(0.25, 0.25, 0.25, 0.25)

theta1 <- c(0, 0, 2, 50, 3.5, 65, 5, 100)

theta2 <- c(1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1)

# Create EMMVNORM.Theta object with new call.

Theta <- new("EMMVNORM.Theta", d = d, c = c)

# Set parameters of Theta.
# Weights.

a.w(Theta) <- w

# Means.

a.theta1.all(Theta) <- theta1

# Covariances

a.theta2.all(Theta) <- theta2


# Run EMMIX method.

model <- EMMIX(model = "REBMVNORM", Dataset = list(faithful), Theta = Theta, EMcontrol = EM)

summary(model)
}
}
\keyword{parameter estimation}
