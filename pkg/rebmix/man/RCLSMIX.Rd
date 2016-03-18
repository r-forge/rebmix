\name{RCLSMIX-methods}
\docType{methods}
\alias{RCLSMIX}
\alias{RCLSMIX-methods}
\alias{RCLSMIX,RCLSMIX-method}
\alias{RCLSMIX,RCLSMVNORM-method}
\alias{show,RCLSMIX-method}
\alias{show,RCLSMVNORM-method}
\alias{summary,RCLSMIX-method}
\alias{summary,RCLSMVNORM-method}
\title{
Predicts Class Membership Based Upon a Model Trained by REBMIX
}
\description{
Returns as default the RCLSMIX algorithm output for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson or Dirac component densities. If \code{model} equals \code{"RCLSMVNORM"} output for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices is returned.
}
\usage{
\S4method{RCLSMIX}{RCLSMIX}(model = "RCLSMIX", x = list(), Dataset = data.frame(),
        Zt = factor(), ...)
## ... and for other signatures
\S4method{summary}{RCLSMIX}(object, ...)
## ... and for other signatures     
}
\arguments{
\item{model}{
see Methods section below.
}
\item{x}{
a list of objects of class \code{REBMIX} of length \eqn{o} obtained by running \code{\link[rebmix]{REBMIX}} on \eqn{g = 1, \ldots, s} train datasets \eqn{Y_{\mathrm{train}g}} all of length \eqn{n_{\mathrm{train}g}}.
For the train datasets the corresponding class membership \eqn{\bm{\Omega}_{g}} is known. This yields
\eqn{n_{\mathrm{train}} = \sum_{g = 1}^{s} n_{\mathrm{train}g}}, while \eqn{Y_{\mathrm{train}q} \cap Y_{\mathrm{train}g} = \emptyset} for all \eqn{q \neq g}.
Each object in the list corresponds to one chunk, e.g., \eqn{(y_{1j}, y_{3j})^{\top}}. The default value is \code{list()}.
}
\item{Dataset}{
a data frame containing test dataset \eqn{Y_{\mathrm{test}}} of length \eqn{n_{\mathrm{test}}}. For the test dataset the corresponding class membership \eqn{\bm{\Omega}_{g}} is not known.
The default value is \code{data.frame()}.
}
\item{Zt}{
a factor of true class membership \eqn{\bm{\Omega}_{g}} for the test dataset. The default value is \code{factor()}.
}
\item{object}{
see Methods section below.
}
\item{\dots}{
currently not used; additional arguments affecting the summary produced.
}
}
\value{
Returns an object of class \code{RCLSMIX} or \code{RCLSMVNORM}.
}
\section{Methods}{
\describe{
\item{\code{signature(model = "RCLSMIX")}}{a character giving the default class name \code{"RCLSMIX"} for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson or Dirac component densities.}
\item{\code{signature(model = "RCLSMVNORM")}}{a character giving the class name \code{"RCLSMVNORM"} for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices.}
\item{\code{signature(object = "RCLSMIX")}}{an object of class \code{RCLSMIX}.}
\item{\code{signature(object = "RCLSMVNORM")}}{an object of class \code{RCLSMVNORM}.}
}
}
\author{Marko Nagode}
\references{
R. O. Duda and P. E. Hart. Pattern Classification and Scene Analysis. John Wiley & Sons, New
York, 1973.
}
\examples{
\dontrun{
devAskNewPage(ask = TRUE)

data("adult")

# Find complete cases.

adult <- adult[complete.cases(adult), ]

# Map metric attributes.

adult[["Age"]] <- ordered(cut(adult[["Age"]], 2000))
adult[["Hours.Per.Week"]] <- ordered(cut(adult[["Hours.Per.Week"]], 2000))
adult[["Capital.Gain"]] <- ordered(cut(adult[["Capital.Gain"]], 2000))
adult[["Capital.Loss"]] <- ordered(cut(adult[["Capital.Loss"]], 2000))

# Show level attributes for binary and discrete variables.

levels(adult[["Type"]])
levels(adult[["Workclass"]])
levels(adult[["Education"]])
levels(adult[["Marital.Status"]])
levels(adult[["Occupation"]])
levels(adult[["Relationship"]])
levels(adult[["Race"]])
levels(adult[["Sex"]])
levels(adult[["Native.Country"]])
levels(adult[["Income"]])

# Replace levels with numbers.

adult <- as.data.frame(data.matrix(adult))

# Levels should start with 0 for discrete distributions except for 
# Dirac distribution.

f <- c("Type", "Age", "Workclass", "Education", "Marital.Status", 
  "Occupation", "Relationship", "Race", "Sex", "Capital.Gain", 
  "Capital.Loss", "Hours.Per.Week", "Native.Country", "Income")

adult[, f] <- adult[, f] - 1

# Split adult dataset into two train subsets for two Incomes
# and remove Type and Income columns.

trainle50k <- subset(adult, subset = (Type == 1) & (Income == 0), 
  select = c(-Type, -Income))
traingt50k <- subset(adult, subset = (Type == 1) & (Income == 1), 
  select = c(-Type, -Income))

trainall <- subset(adult, subset = Type == 1, select = c(-Type, -Income))

train <- as.factor(subset(adult, subset = Type == 1, select = c(Income))[, 1])

# Extract test dataset form adult dataset and remove Type 
# and Income columns.

testle50k <- subset(adult, subset = (Type == 0) & (Income == 0), 
  select = c(-Type, -Income))
testgt50k <- subset(adult, subset = (Type == 0) & (Income == 1), 
  select = c(-Type, -Income))

testall <- subset(adult, subset = Type == 0, select = c(-Type, -Income))

test <- as.factor(subset(adult, subset = Type == 0, select = c(Income))[, 1])

# Estimate number of components, component weights and component 
# parameters for Naive Bayes.

Variables <- c("discrete", "discrete", "continuous", "discrete",
 "discrete", "discrete", "discrete", "discrete", "discrete",
 "discrete", "discrete", "discrete", "discrete", "discrete")

pdf <- c("Dirac", "Dirac", "gamma", "Dirac",
  "Dirac", "Dirac", "Dirac", "Dirac", "Dirac",
  "Dirac", "Dirac", "Dirac", "Dirac", "Dirac")

n <- range(nrow(trainle50k), nrow(traingt50k))

K <- c(as.integer(1 + log2(sum(n[1]))), # Minimum v follows Sturges rule.
  as.integer(10 * log10(n[2]))) # Maximum v follows log10 rule.

K <- list(1, 1, K[1]:K[2], 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

adultest <- list(0)

for (i in c(1:14)) {
  adultest[[i]] <- REBMIX(Dataset = list(as.data.frame(trainle50k[, i]), 
    as.data.frame(traingt50k[, i])),
    Preprocessing = "histogram",
    cmax = if (Variables[i] == "continuous") 8 else 100,
    Criterion = if (Variables[i] == "continuous") "BIC" else "D",
    pdf = pdf[i],
    K = K[[i]])

    plot(adultest[[i]], pos = 1)
    plot(adultest[[i]], pos = 2)
}

# Best-first feature subset selection.

c <- NULL; rvs <- 1:14; Error <- 1.0

for (i in 1:14) {
  k <- NA

  for (j in rvs) {
    predictive <- RCLSMIX(x = adultest[c(c, j)],
      Dataset = as.data.frame(trainall[, c(c, j)]),
      Zt = train) 

    if (predictive@Error < Error) {
      Error <- predictive@Error; k <- j
    }
  }

  if (is.na(k)) {
    break
  }
  else {
    c <- c(c, k); rvs <- rvs[-which(rvs == k)]
  }
}

# Error on train dataset.

Error

# Selected features.

c

predictive <- RCLSMIX(x = adultest[c],
  Dataset = as.data.frame(testall[, c]), 
  Zt = test)
  
predictive  

summary(predictive)

# Plot selected features.

plot(predictive, nrow = 5, ncol = 2)
}
}
\keyword{classification}
