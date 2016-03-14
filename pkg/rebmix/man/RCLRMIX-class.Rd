\name{RCLRMIX-class}
\docType{class}
\alias{RCLRMIX-class}
\alias{RCLRMVNORM-class}
\title{Class \code{"RCLRMIX"}}
\description{Object of class \code{RCLRMIX}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RCLRMIX", ...)}.
}
\section{Slots}{
\describe{
\item{\code{x}:}{
an object of class \code{REBMIX}.
}
\item{\code{pos}:}{
a desired row number in \code{x@summary} to be clustered. The default value is \code{1}.
}
\item{\code{Zt}:}{
a factor of true cluster membership.
}
\item{\code{Zp}:}{
a factor of predictive cluster membership.
}
}
}
\author{Marko Nagode}
\keyword{clusters}
