\name{histogram-methods}
\alias{histogram}
\alias{histogram-methods}
\alias{histogram,REBMIX-method}
\alias{histogram,REBMVNORM-method}
\title{
Histogram Calculation
}
\description{
Returns an object of class \code{Histogram}. The method can be called recursively. This way more than one dataset can be
binned into one histogram.  
}
\usage{
\S4method{histogram}{REBMIX}(x = NULL, Dataset = data.frame(),
          K = numeric(), ymin = numeric(), ymax = numeric(), ...)
## ... and for other signatures
}
\arguments{
\item{x}{
an object of class \code{Histogram}.
}
\item{Dataset}{
a data frame of size \eqn{n \times d} containing \emph{d}-dimensional dataset. Each of the \eqn{d} columns
represents one random variable. Number of observations \eqn{n} equals the number of rows in the dataset.
}
\item{K}{
an integer or a vector of length \eqn{d} containing numbers of bins \eqn{v}.
}
\item{ymin}{
a vector of length \eqn{d} containing minimum observations.
}
\item{ymax}{
a vector of length \eqn{d} containing maximum observations.
}
\item{\dots}{
currently not used.
}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "REBMIX")}}{an object of class \code{REBMIX}.}
\item{\code{signature(x = "REBMVNORM")}}{an object of class \code{REBMVNORM}.}
}
}
\author{Marko Nagode}
\examples{
# Generate multivariate normal datasets.

n <- c(750, 1000)

Theta <- new("RNGMVNORM.Theta", c = 2, d = 2)

a.theta1(Theta, 1) <- c(8, 6)
a.theta1(Theta, 2) <- c(6, 8)
a.theta2(Theta, 1) <- c(8, 2, 2, 4)
a.theta2(Theta, 2) <- c(2, 1, 1, 4)

sim2d <- RNGMIX(model = "RNGMVNORM", 
  Dataset.name = paste("sim2d_", 1:5, sep = ""),
  rseed = -1,
  n = n,
  Theta = a.Theta(Theta))

# Calculate optimal numbers of bins.

opt.k <- optbins(Dataset = sim2d@Dataset,
  Rule = "Knuth equal",
  ymin = sim2d@ymin,
  ymax = sim2d@ymax,
  kmin = 2, 
  kmax = 20)

opt.k

# Create object of class EM.Control.

EM <- new("EM.Control", strategy = "exhaustive", variant = "EM",
  acceleration = "fixed", acceleration.multiplier = 1.0, tolerance = 1.0E-4,
  maximum.iterations = 1000)

# Estimate number of components, component weights and component parameters.

sim2dest <- REBMIX(model = "REBMVNORM", 
  Dataset = a.Dataset(sim2d),
  Preprocessing = "h",
  cmax = 10,
  ymin = a.ymin(sim2d),
  ymax = a.ymax(sim2d),
  K = opt.k,
  Criterion = "BIC",
  EMcontrol = EM)

# Plot finite mixture.

plot(sim2dest, pos = 3, nrow = 4, what = c("pdf", "marginal pdf", "IC"))

# Estimate number of components, component weights and component 
# parameters for well known Iris dataset.

Dataset <- list(iris[, c(1:4)])

# Calculate optimal numbers of bins using non-equal number of bins in each dimension.

opt.k <- optbins(Dataset = Dataset,
  Rule = "Knuth unequal",
  kmin = 2, 
  kmax = 20)

opt.k

# Estimate number of components, component weights and component parameters.

irisest <- REBMIX(model = "REBMVNORM", 
  Dataset = Dataset,
  Preprocessing = "h",
  cmax = 10,
  K = opt.k,
  Criterion = "BIC",
  EMcontrol = EM)
  
irisest
}
\keyword{auxiliary}
