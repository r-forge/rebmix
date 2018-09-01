\name{iris}
\docType{data}
\alias{iris}
\title{Iris Data Set}
\description{
This is perhaps the best known database to be found in the pattern recognition literature. 
Fisher's paper is a classic in the field and is referenced frequently to this day. The data set contains 3 classes of 50 instances each, 
where each class refers to a type of iris plant. One class is linearly separable from the other 2; the latter are NOT linearly separable from each other.
}
\usage{
data("iris")
}
\format{
\code{iris} is a data frame with 150 cases (rows) and 5 variables (columns) named:
\enumerate{
\item{
\code{Sepal.Length} continuous.
}
\item{
\code{Sepal.Width} continuous.
}
\item{
\code{Petal.Length} continuous.
}
\item{
\code{Petal.Width} continuous.
}
\item{
\code{Class} discrete \code{iris-setosa}, \code{iris-versicolour} or \code{iris-virginica}.
}
}
}
\source{
A. Asuncion and D. J. Newman. Uci machine learning repository, 2007. \url{http://archive.ics.uci.edu/ml}.
}
\references{
R. A. Fisher. The use of multiple measurements in taxonomic problems. Annals of Eugenics,
7(2):179-188, 1936.
}
\examples{
\dontrun{
devAskNewPage(ask = TRUE)

data("iris")

# Show level attributes.

levels(iris[["Class"]])

# Split dataset into train (75%) and test (25%) subsets.

set.seed(5)

Iris <- split(p = 0.75, Dataset = iris, class = 5)

# Estimate number of components, component weights and component 
# parameters for train subsets.

n <- range(Iris@ntrain)

K <- c(as.integer(1 + log2(n[1])), # Minimum v follows Sturges rule.
  as.integer(10 * log10(n[2]))) # Maximum v follows log10 rule.

K <- c(floor(K[1]^(1/4)), ceiling(K[2]^(1/4)))

irisest <- REBMIX(model = "REBMVNORM",
  Dataset = Iris@train,
  Preprocessing = "kernel density estimation",
  cmax = 10,
  Criterion = "ICL-BIC",
  pdf = rep("normal", 4),
  K = K[1]:K[2])

plot(irisest, pos = 1, nrow = 3, ncol = 2, what = c("den"))
plot(irisest, pos = 2, nrow = 3, ncol = 2, what = c("den"))
plot(irisest, pos = 3, nrow = 3, ncol = 2, what = c("den"))

# Selected chunks.

iriscla <- RCLSMIX(model = "RCLSMVNORM", 
  x = list(irisest),
  Dataset = Iris@test,
  Zt = Iris@Zt)

iriscla

summary(iriscla)

# Plot selected chunks.

plot(iriscla, nrow = 3, ncol = 2)
}
}
\keyword{datasets}
