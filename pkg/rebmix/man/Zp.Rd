\name{Zp-methods}
\alias{Zp}
\alias{Zp-methods}
\alias{Zp,RCLRMIX-method}
\alias{Zp,RCLRMVNORM-method}
\title{
Predictive Cluster Membership
}
\description{
Returns the factor of predictive cluster membership at \code{s}.
}
\usage{
\S4method{Zp}{RCLRMIX}(x = NULL, s = expression(c))
## ... and for other signatures
}
\arguments{
\item{x}{
see Methods section below.
}
\item{s}{
a desired number of clusters for which the predictive cluster membership is calculated. The default value is \code{expression(c)}.
}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "RCLRMIX")}}{an object of class \code{RCLRMIX}.}
\item{\code{signature(x = "RCLRMVNORM")}}{an object of class \code{RCLRMVNORM}.}
}
}
\author{Marko Nagode}
\keyword{clustering}
