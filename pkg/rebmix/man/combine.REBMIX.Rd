\name{combine-methods}
\docType{methods}
\alias{combine}
\alias{combine-methods}
\alias{combine,REBMIX-method}
\alias{combine,REBMVNORM-method}
\alias{show,REBMIX.combine-method}
\alias{show,REBMVNORM.combine-method}
\title{
Combining Mixture Components for Clustering
}
\description{
Returns as default the combined mixture component output for mixtures of conditionally independent normal, 
lognormal, Weibull, gamma, binomial, Poisson or Dirac component densities. If 
\code{x} is of class \code{RNGMVNORM} the combined mixture component output for mixtures of multivariate normal 
component densities with unrestricted variance-covariance matrices is returned.
} 
\usage{ 
\S4method{combine}{REBMIX}(x = NULL, pos = 1, ...)
## ... and for other signatures
}
\arguments{
\item{x}{
see Methods section below.
}
\item{pos}{
a desired row number in \code{x@summary} to be combined. The default value is \code{1}.
}
\item{\dots}{
currently not used.
}
}
\value{
Returns an object of class \code{REBMIX.combine} or \code{REBMVNORM.combine}.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "REBMIX")}}{an object of class \code{REBMIX} for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson or Dirac component densities.}
\item{\code{signature(x = "REBMVNORM")}}{an object of class \code{REBMVNORM} for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices.}
}
}
\author{Marko Nagode}
\references{
J. P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo. Combining mixture components for clustering.
Journal of Computational and Graphical Statistics, 19(2):332-353, 2010. \url{http://dx.doi.org/10.1198/jcgs.2010.08111}
}
\examples{
a <- 1
}
\keyword{clustering}
