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
\examples{
devAskNewPage(ask = TRUE)

# Generate normal dataset.

n <- c(500, 200, 400)

Theta <- list(pdf1 = rep("normal", 2),
  theta1.1 = c(3, 10),
  theta2.1 = c(3, 0.3, 0.3, 2),
  pdf2 = rep("normal", 2),
  theta1.2 = c(8, 6),
  theta2.2 = c(5.7, -2.3, -2.3, 3.5),
  pdf3 = rep("normal", 2),
  theta1.3 = c(12, 11),
  theta2.3 = c(2, 1, 1, 2))

normal <- RNGMIX(model = "RNGMVNORM", Dataset.name = "normal_1", n = n, Theta = Theta)

# Number of classes or nearest neighbours to be processed.

K <- c(as.integer(1 + log2(sum(n))), # Minimum v follows Sturges rule.
  as.integer(10 * log10(sum(n)))) # Maximum v follows log10 rule.

# Estimate number of components, component weights and component parameters.

normalest <- REBMIX(model = "REBMVNORM",
  Dataset = normal@Dataset,
  Preprocessing = "histogram",
  cmax = 6,
  Criterion = "BIC",
  pdf = rep("normal", 2),
  K = K[1]:K[2])

summary(normalest)

# Plot finite mixture.

plot(normalest)

# Cluster dataset.

normalclu <- RCLRMIX(model = "RCLRMVNORM", x = normalest, Zt = normal@Zt)

# Plot clusters.

plot(normalclu)
}
\keyword{classes}
