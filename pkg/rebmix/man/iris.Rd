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

# Show level attributes discrete variables.

levels(iris[["Class"]])

# Split iris dataset into three subsets for three Classes
# and remove Class column.

iris_set <- subset(iris, subset = Class == "iris-setosa", select = c(-Class))
iris_ver <- subset(iris, subset = Class == "iris-versicolor", select = c(-Class))
iris_vir <- subset(iris, subset = Class == "iris-virginica", select = c(-Class))

# Split datasets into train (75%) and test (25%) subsets.

set.seed(5)

Prob <- 0.75

n_set <- nrow(iris_set); s_set <- sample.int(n = n_set, size = as.integer(n_set * Prob))

iris_set_train <- iris_set[s_set,]; iris_set_test <- iris_set[-s_set,]

n_ver <- nrow(iris_ver); s_ver <- sample.int(n = n_ver, size = as.integer(n_ver * Prob))

iris_ver_train <- iris_ver[s_ver,]; iris_ver_test <- iris_ver[-s_ver,]

n_vir <- nrow(iris_vir); s_vir <- sample.int(n = n_vir, size = as.integer(n_vir * Prob))

iris_vir_train <- iris_vir[s_vir,]; iris_vir_test <- iris_vir[-s_vir,]

iris_test = rbind(iris_set_test, iris_ver_test, iris_vir_test)

Zt <- factor(c(rep(0, nrow(iris_set_test)), 
  rep(1, nrow(iris_ver_test)), 
  rep(2, nrow(iris_vir_test))))

# Estimate number of components, component weights and component 
# parameters for train subsets.

n <- range(nrow(iris_set_train), nrow(iris_ver_train), nrow(iris_vir_train))

K <- c(as.integer(1 + log2(n[1])), # Minimum v follows Sturges rule.
  as.integer(10 * log10(n[2]))) # Maximum v follows log10 rule.

K <- c(floor(K[1]^(1/4)), ceiling(K[2]^(1/4)))

irisest <- REBMIX(model = "REBMVNORM",
  Dataset = list(iris_set_train = iris_set_train,
                 iris_ver_train = iris_ver_train,
                 iris_vir_train = iris_vir_train),
  Preprocessing = "Parzen window",
  cmax = 10,
  Criterion = "ICL-BIC",
  pdf = rep("normal", 4),
  K = K[1]:K[2])

plot(irisest, pos = 1, nrow = 3, ncol = 2, what = c("den"))
plot(irisest, pos = 2, nrow = 3, ncol = 2, what = c("den"))
plot(irisest, pos = 3, nrow = 3, ncol = 2, what = c("den"))

# Selected features.

predictive <- RCLSMIX(model = "RCLSMVNORM", 
  x = list(irisest),
  Dataset = iris_test,
  Zt = Zt)

predictive

summary(predictive)

# Plot selected features.

plot(predictive, nrow = 3, ncol = 2)
}
}
\keyword{datasets}
