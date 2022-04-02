\name{mapclusters-methods}
\alias{mapclusters}
\alias{mapclusters-methods}
\alias{mapclusters,RCLRMIX-method}
\alias{mapclusters,RCLRMVNORM-method}
\title{
Map Clusters
}
\description{
Returns a factor of predictive cluster membership for dataset.
}
\usage{
\S4method{mapclusters}{RCLRMIX}(x = NULL, Dataset = data.frame(),
            s = expression(c), ...)
## ... and for other signatures
}
\arguments{
\item{x}{
see Methods section below.
}
\item{Dataset}{
a data frame of size \eqn{n \times d} containing \emph{d}-dimensional dataset. Each of the \eqn{d} columns
represents one random variable. Number of observations \eqn{n} equal the number of rows in the dataset.
}
\item{s}{
a desired number of clusters to be created. The default value is \code{expression(c)}.
}
\item{\dots}{
currently not used.
}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "RCLRMIX")}}{an object of class \code{RCLRMIX}.}
\item{\code{signature(x = "RCLRMVNORM")}}{an object of class \code{RCLRMVNORM}.}
}
}
\author{Marko Nagode, Branislav Panic}
\examples{
devAskNewPage(ask = TRUE)

# Generate normal dataset.

n <- c(50, 20, 40)

Theta <- new("RNGMVNORM.Theta", c = 3, d = 2)

a.theta1(Theta, 1) <- c(3, 10)
a.theta1(Theta, 2) <- c(8, 6)
a.theta1(Theta, 3) <- c(12, 11)
a.theta2(Theta, 1) <- c(3, 0.3, 0.3, 2)
a.theta2(Theta, 2) <- c(5.7, -2.3, -2.3, 3.5)
a.theta2(Theta, 3) <- c(2, 1, 1, 2)

normal <- RNGMIX(model = "RNGMVNORM", Dataset.name = paste("normal_", 1:10, sep = ""), n = n, Theta = a.Theta(Theta))

# Convert all datasets to single histogram.

hist <- NULL

n <- length(normal@Dataset)

hist <- fhistogram(Dataset = normal@Dataset[[1]], K = c(10, 10), ymin = a.ymin(normal), ymax = a.ymax(normal))

for (i in 2:n) {
  hist <- fhistogram(x = hist, Dataset = normal@Dataset[[i]], shrink = i == n)
}

# Estimate number of components, component weights and component parameters.

normalest <- REBMIX(model = "REBMVNORM",
  Dataset = list(hist),
  Preprocessing = "histogram",
  cmax = 6,
  Criterion = "BIC")

summary(normalest)

# Plot finite mixture.

plot(normalest)

# Cluster dataset.

normalclu <- RCLRMIX(model = "RCLRMVNORM", x = normalest)

# Plot clusters.

plot(normalclu)

summary(normalclu)

# Map clusters.

Zp <- mapclusters(x = normalclu, Dataset = a.Dataset(normal, 4))

Zt <- a.Zt(normal)

Zp

Zt
}
\keyword{clustering}
