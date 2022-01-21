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
\S4method{EMMIX}{REBMIX}(model = "REBMIX", Dataset = list(), 
       Theta = NULL, EMcontrol = NULL, ...)
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
an object of class \code{"EMMIX.Theta"} or \code{"EMMVNORM.Theta"}.
}
\item{EMcontrol}{
an object of class \code{EM.Control}.
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

# Create EMMVNORM.Theta object with new call.

Theta <- new("EMMVNORM.Theta", d = d, c = c)

# Set parameters of Theta.
# Weights.

a.w(Theta) <- c(0.5, 0.5)

# Means.

a.theta1.all(Theta) <- c(2.0, 55.0, 4.5, 80.0)

# Covariances.

a.theta2.all(Theta) <- c(1, 0, 0, 1, 1, 0, 0, 1)

# Run EMMIX method.

model <- EMMIX(model = "REBMVNORM", Dataset = list(faithful), Theta = Theta)

# show.

model

# summary.

summary(model)

# plot.

plot(model, nrow = 3, ncol = 2, what = c("pdf", "marginal pdf", "marginal cdf"))

# Create EMMIX.Theta object with new call.

Theta <- new("EMMIX.Theta", c = c, pdf = c("normal", "normal"))

# Set parameters of Theta.
# Weights.

a.w(Theta) <- c(0.5, 0.5)

# Means.

a.theta1.all(Theta) <- c(2.0, 55.0, 4.5, 80.0)

# Covariances.

a.theta2.all(Theta) <- c(1, 1, 1, 1)

# Run EMMIX method.

model <- EMMIX(Dataset = list(faithful), Theta = Theta)

# show.

model

# summary.

summary(model)

# plot.

plot(model, nrow = 3, ncol = 2, what = c("pdf", "marginal pdf", "marginal cdf"))
}
}
\keyword{parameter estimation}
