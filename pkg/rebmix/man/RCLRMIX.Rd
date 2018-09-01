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
\alias{g.Zp}
\alias{g.p}
\alias{g.pi}
\alias{g.P}
\alias{g.tau}
\alias{g.from}
\alias{g.to}
\alias{g.EN}
\alias{g.ED}
\alias{g.pos,RCLRMIX-method}
\alias{g.Zt,RCLRMIX-method}
\alias{g.Zp,RCLRMIX-method}
\alias{g.c,RCLRMIX-method}
\alias{g.p,RCLRMIX-method}
\alias{g.pi,RCLRMIX-method}
\alias{g.P,RCLRMIX-method}
\alias{g.tau,RCLRMIX-method}
\alias{g.prob,RCLRMIX-method}
\alias{g.from,RCLRMIX-method}
\alias{g.to,RCLRMIX-method}
\alias{g.EN,RCLRMIX-method}
\alias{g.ED,RCLRMIX-method}
\alias{g.pos,RCLRMVNORM-method}
\alias{g.Zt,RCLRMVNORM-method}
\alias{g.Zp,RCLRMVNORM-method}
\alias{g.c,RCLRMVNORM-method}
\alias{g.p,RCLRMVNORM-method}
\alias{g.pi,RCLRMVNORM-method}
\alias{g.P,RCLRMVNORM-method}
\alias{g.tau,RCLRMVNORM-method}
\alias{g.prob,RCLRMVNORM-method}
\alias{g.from,RCLRMVNORM-method}
\alias{g.to,RCLRMVNORM-method}
\alias{g.EN,RCLRMVNORM-method}
\alias{g.ED,RCLRMVNORM-method}
\title{
Predicts Cluster Membership Based Upon a Model Trained by REBMIX
}
\description{
Returns as default the RCLRMIX algorithm output for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson, Dirac or von Mises component densities, following the methodology proposed in the article cited in the references. If \code{model} equals \code{"RCLRMVNORM"} output for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices is returned.
}
\usage{
\S4method{RCLRMIX}{RCLRMIX}(model = "RCLRMIX", x = NULL, pos = 1, Zt = factor(), ...)
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
\item{pos}{
a desired row number in \code{x@summary} to be clustered. The default value is \code{1}.
}
\item{Zt}{
a factor of true cluster membership. The default value is \code{factor()}.
}
\item{object}{
see Methods section below.
}
\item{\dots}{
currently not used; additional arguments affecting the summary produced.
}
}
\value{
Returns an object of class \code{RCLRMIX} or \code{RCLRMVNORM}.
}
\section{Methods}{
\describe{
\item{\code{signature(model = "RCLRMIX")}}{a character giving the default class name \code{"RCLRMIX"} for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson, Dirac or von Mises component densities.}
\item{\code{signature(model = "RCLRMVNORM")}}{a character giving the class name \code{"RCLRMVNORM"} for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices.}
\item{\code{signature(object = "RCLRMIX")}}{an object of class \code{RCLRMIX}.}
\item{\code{signature(object = "RCLRMVNORM")}}{an object of class \code{RCLRMVNORM}.}
}
}
\author{Marko Nagode}
\references{
J. P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo. Combining mixture components for clustering.
Journal of Computational and Graphical Statistics, 19(2):332-353, 2010. \url{http://dx.doi.org/10.1198/jcgs.2010.08111}
}
\examples{
devAskNewPage(ask = TRUE)

# Generate Poisson dataset.

n <- c(500, 200, 400)

Theta <- list(pdf1 = "Poisson",
  theta1.1 = 3,
  theta2.1 = NA,
  pdf2 = "Poisson",
  theta1.2 = 12,
  theta2.2 = NA,
  pdf3 = "Poisson",
  theta1.3 = 36,
  theta2.3 = NA)

poisson <- RNGMIX(Dataset.name = "Poisson_1", n = n, Theta = Theta)

# Estimate number of components, component weights and component parameters.

poissonest <- REBMIX(Dataset = poisson@Dataset,
  Preprocessing = "histogram",
  cmax = 6,
  Criterion = "BIC",
  pdf = rep("Poisson", 1),
  K = 1)

summary(poissonest)

# Plot finite mixture.

plot(poissonest)

# Cluster dataset.

poissonclu <- RCLRMIX(x = poissonest, Zt = poisson@Zt)

summary(poissonclu)

# Plot clusters.

plot(poissonclu)

summary(poissonclu)
}
\keyword{clustering}
