\name{optbins-methods}
\alias{optbins}
\alias{optbins-methods}
\alias{optbins,list-method}
\title{
Optimal Numbers of Bins Calculation
}
\description{
Returns the matrix of size \eqn{n_{\mathrm{D}} \times d} containing optimal numbers of bins \eqn{v_{1}, \ldots, v_{d}} for all processed datasets.
}
\usage{
\S4method{optbins}{list}(Dataset = list(), Rule = "Knuth equal", y0 = numeric(), 
        ymin = numeric(), ymax = numeric(), kmin = numeric(),
        kmax = numeric(), ...)
## ... and for other signatures
}
\arguments{
\item{Dataset}{
a list of length \eqn{n_{\mathrm{D}}} of data frames of size \eqn{n \times d} containing \emph{d}-dimensional datasets. Each of the \eqn{d} columns
represents one random variable. Numbers of observations \eqn{n} equal the number of rows in the datasets.
}
\item{Rule}{
a character giving the histogram binning rule. One of \code{"Sturges"}, \code{"Log10"}, \code{"RootN"}, default \code{"Knuth equal"} or \code{"Knuth unequal"}.
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
\item{kmin}{
lower limit of the number of bins. The default value is \code{numeric()}.
}
\item{kmax}{
upper limit of the number of bins. The default value is \code{numeric()}.
}
\item{\dots}{
currently not used.
}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "list")}}{a list of data frames.}
}
}
\author{Branislav Panic, Marko Nagode}
\references{
K. K. Knuth. Optimal data-based binning for histograms and histogram-based probability density models. 
Digital Signal Processing, 95:102581, 2019.
\url{https://doi.org/10.1016/j.dsp.2019.102581}.
}
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
  Criterion = "BIC")

# Plot finite mixture.

plot(sim2dest, pos = 3, nrow = 4, what = c("pdf", "marginal pdf", "IC"))
}
\keyword{auxiliary}
