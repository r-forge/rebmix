\name{RCLRMIX-methods}
\docType{methods}
\alias{RCLRMIX}
\alias{RCLRMIX-methods}
\alias{RCLRMIX,RCLRMIX-method}
\alias{RCLRMIX,RCLRMVNORM-method}
\alias{show,RCLRMIX-method}
\alias{show,RCLRMVNORM-method}
\alias{summary,RCLRMIX-method}
\alias{summary,RCLRMVNORM-method}
\title{
Predicts Cluster Membership Based Upon a Model Trained by REBMIX
}
\description{
Returns as default the RCLRMIX algorithm output for mixtures of conditionally independent normal, lognormal, Weibull, gamma, Gumbel, binomial, Poisson, Dirac, uniform or von Mises component densities, following the methodology proposed in the article cited in the references. If \code{model} equals \code{"RCLRMVNORM"} output for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices is returned.
}
\usage{
\S4method{RCLRMIX}{RCLRMIX}(model = "RCLRMIX", x = NULL, Dataset = NULL, 
        pos = 1, Zt = factor(), Rule = character(), ...)
## ... and for other signatures
\S4method{summary}{RCLRMIX}(object, ...)
## ... and for other signatures
}
\arguments{
\item{model}{
see Methods section below.
}
\item{x}{
an object of class \code{REBMIX}.
}
\item{Dataset}{
a data frame or an object of class \code{Histogram} to be clustered.
}
\item{pos}{
a desired row number in \code{x@summary} for which the clustering is performed. The default value is \code{1}.
}
\item{Zt}{
a factor of true cluster membership. The default value is \code{factor()}.
}
\item{Rule}{
a character containing the merging rule. One of \code{"Entropy"} or \code{"Demp"}. The default value is \code{"Entropy"}.
}
\item{object}{
see Methods section below.
}
\item{\dots}{
currently not used.
}
}
\value{
Returns an object of class \code{RCLRMIX} or \code{RCLRMVNORM}.
}
\section{Methods}{
\describe{
\item{\code{signature(model = "RCLRMIX")}}{a character giving the default class name \code{"RCLRMIX"} for mixtures of conditionally independent normal, lognormal, Weibull, gamma, Gumbel, binomial, Poisson, Dirac, uniform or von Mises component densities.}
\item{\code{signature(model = "RCLRMVNORM")}}{a character giving the class name \code{"RCLRMVNORM"} for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices.}
\item{\code{signature(object = "RCLRMIX")}}{an object of class \code{RCLRMIX}.}
\item{\code{signature(object = "RCLRMVNORM")}}{an object of class \code{RCLRMVNORM}.}
}
}
\author{Marko Nagode}
\references{
J. P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo. Combining mixture components for clustering.
Journal of Computational and Graphical Statistics, 19(2):332-353, 2010. \doi{10.1198/jcgs.2010.08111}
}
\examples{
devAskNewPage(ask = TRUE)

# Generate Poisson dataset.

n <- c(500, 200, 400)

Theta <- new("RNGMIX.Theta", c = 3, pdf = "Poisson")

a.theta1(Theta) <- c(3, 12, 36)

poisson <- RNGMIX(Dataset.name = "Poisson_1", n = n, Theta = a.Theta(Theta))

# Estimate number of components, component weights and component parameters.

EM <- new("EM.Control", strategy = "exhaustive")

poissonest <- REBMIX(Dataset = a.Dataset(poisson),
  Preprocessing = "histogram",
  cmax = 6,
  Criterion = "BIC",
  pdf = rep("Poisson", 1),
  EMcontrol = EM)

summary(poissonest)

# Plot finite mixture.

plot(poissonest)

# Cluster dataset.

poissonclu <- RCLRMIX(x = poissonest, Zt = a.Zt(poisson))

summary(poissonclu)

# Plot clusters.

plot(poissonclu)

# Create new dataset.

Dataset <- sample.int(n = 50, size = 10, replace = TRUE)

Dataset <- as.data.frame(Dataset)

# Cluster the dataset.

poissonclu <- RCLRMIX(x = poissonest, Dataset = Dataset, Rule = "Demp")

a.Dataset(poissonclu)
}
\keyword{clustering}
