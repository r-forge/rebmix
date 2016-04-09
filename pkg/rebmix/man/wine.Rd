\name{wine}
\docType{data}
\alias{wine}
\title{Wine Recognition Data}
\description{
These data are the results of a chemical analysis of wines grown in the same region in Italy but derived from three
different cultivars (1-3). The analysis determined the quantities of 13 constituents: alcohol, malic acid, ash, alcalinity of ash,
magnesium, total phenols, flavanoids, nonflavanoid phenols, proanthocyanins, colour intensity, hue, OD280/OD315 of diluted wines,
and proline found in each of the three types of the wines. The number of instances in classes 1 to 3 is 59, 71 and 48, respectively.
}
\usage{
data("wine")
}
\format{
\code{wine} is a data frame with 178 cases (rows) and 14 variables (columns) named:
\enumerate{
\item{
\code{Alcohol} continuous.
}
\item{
\code{Malic.Acid} continuous.
}
\item{
\code{Ash} continuous.
}
\item{
\code{Alcalinity.of.Ash} continuous.
}
\item{
\code{Magnesium} continuous.
}
\item{
\code{Total.Phenols} continuous.
}
\item{
\code{Flavanoids} continuous.
}
\item{
\code{Nonflavanoid.Phenols} continuous.
}
\item{
\code{Proanthocyanins} continuous. 
}
\item{
\code{Color.Intensity} continuous.
}
\item{
\code{Hue} continuous.
}
\item{
\code{OD280.OD315.of.Diluted.Wines} continuous.
}
\item{
\code{Proline} continuous.
}
\item{
\code{Cultivar} discrete \code{1}, \code{2} or \code{3}.
}
}
}
\source{
A. Asuncion and D. J. Newman. Uci machine learning repository, 2007. \url{http://archive.ics.uci.edu/ml}.
}
\references{
S. J. Roberts, R. Everson, and I. Rezek. Maximum certainty data partitioning. Pattern Recognition,
33(5):833-839, 2000. \url{http://dx.doi.org/10.1016/S0031-3203(99)00086-2}.
}
\examples{
\dontrun{
devAskNewPage(ask = TRUE)

data("wine")

# Show level attributes discrete variables.

levels(factor(wine[["Cultivar"]]))

# Split wine dataset into three subsets for three Cultivars
# and remove Cultivar column.

wine_1 <- subset(wine, subset = Cultivar == 1, select = c(-Cultivar))
wine_2 <- subset(wine, subset = Cultivar == 2, select = c(-Cultivar))
wine_3 <- subset(wine, subset = Cultivar == 3, select = c(-Cultivar))

# Split datasets into train (75%) and test (25%) subsets.

set.seed(3)

Prob <- 0.75

n_1 <- nrow(wine_1); s_1 <- sample.int(n = n_1, size = as.integer(n_1 * Prob))

wine_1_train <- wine_1[s_1,]; wine_1_test <- wine_1[-s_1,]

n_2 <- nrow(wine_2); s_2 <- sample.int(n = n_2, size = as.integer(n_2 * Prob))

wine_2_train <- wine_2[s_2,]; wine_2_test <- wine_2[-s_2,]

n_3 <- nrow(wine_3); s_3 <- sample.int(n = n_3, size = as.integer(n_3 * Prob))

wine_3_train <- wine_3[s_3,]; wine_3_test <- wine_3[-s_3,]

wine_test = rbind(wine_1_test, wine_2_test, wine_3_test)

Zt <- factor(c(rep(0, nrow(wine_1_test)), rep(1, nrow(wine_2_test)), rep(2, nrow(wine_3_test))))

# Estimate number of components, component weights and component 
# parameters for train subsets.

n <- range(nrow(wine_1_train), nrow(wine_2_train), nrow(wine_3_train))

K <- c(as.integer(1 + log2(n[1])), # Minimum v follows Sturges rule.
  as.integer(10 * log10(n[2]))) # Maximum v follows log10 rule.

K <- c(floor(K[1]^(1/13)), ceiling(K[2]^(1/13)))

wineest <- REBMIX(model = "REBMVNORM",
  Dataset = list(wine_1_train = wine_1_train,
                 wine_2_train = wine_2_train,
                 wine_3_train = wine_3_train),
  Preprocessing = "Parzen window",
  cmax = 10,
  Criterion = "ICL-BIC",
  pdf = rep("normal", 13),
  K = K[1]:K[2],
  Restraints = "loose")

plot(wineest, pos = 1, nrow = 7, ncol = 6, what = c("den"))
plot(wineest, pos = 2, nrow = 7, ncol = 6, what = c("den"))
plot(wineest, pos = 3, nrow = 7, ncol = 6, what = c("den"))

# Selected features.

predictive <- RCLSMIX(model = "RCLSMVNORM", 
  x = list(wineest),
  Dataset = wine_test,
  Zt = Zt)

predictive

summary(predictive)

# Plot selected features.

plot(predictive, nrow = 7, ncol = 6)
}
}
\keyword{datasets}
