\name{RCLSMIX}
\alias{RCLSMIX}
\title{
Predicts Class Membership Based Upon a Model Trained by REBMIX
}
\description{
Takes a list of objects of class \code{REBMIX}, prior probabilities \code{P} and data frame \code{Dataset} with test classes and returns predictive classes.
}
\usage{
RCLSMIX(x, P = NULL, Dataset = NULL, ...)      
}
\arguments{
\item{x}{
a list of objects of class \code{REBMIX} containing outputs for a finite set of size \eqn{s} of classes \eqn{\Omega_g}.
}
\item{P}{
a vector of length \eqn{s} containing prior probabilites.
}
\item{Dataset}{
a data frame containing test classes.
}
\item{\dots}{
potential further arguments of the method.
}
}
\value{
predictive classes.
}
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

adult[["Capital.Loss"]] <- ordered(cut(adult[["Capital.Loss"]], 2000))
adult[["Capital.Gain"]] <- ordered(cut(adult[["Capital.Gain"]], 2000))

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

# Levels should start with 0 for discrete distributions except for the 
# Dirac distribution.

f <- c("Type", "Workclass", "Education", "Marital.Status", "Occupation", 
  "Relationship", "Race", "Sex", "Native.Country", "Income")

adult[, f] <- adult[, f] - 1

# Split adult dataset into two train subsets for the two Incomes
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

# Calculate prior probabilities.

P <- c(nrow(trainle50k), nrow(traingt50k))

P <- P / sum(P)

# Estimate number of components, component weights and component 
# parameters for Naive Bayes.

Variables <- c("continuous", "discrete", "continuous", "discrete",
 "discrete", "discrete", "discrete", "discrete", "discrete",
 "discrete", "discrete", "discrete", "discrete", "discrete")

pdf <- c("gamma", "Dirac", "gamma", "Dirac",
  "Dirac", "Dirac", "Dirac", "Dirac", "Dirac",
  "Dirac", "Dirac", "Dirac", "Dirac", "Dirac")

K <- list(kseq(6, 20, 0.05), 1, kseq(6, 20, 0.05), 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

ymin <- as.numeric(apply(trainall, 2, min))
ymax <- as.numeric(apply(trainall, 2, max))

adultest <- list()

for (i in c(1:14)) {
  adultest[[i]] <- REBMIX(Dataset = list(as.data.frame(trainle50k[, i]), 
    as.data.frame(traingt50k[, i])),
    Preprocessing = "histogram",
    cmax = if (Variables[i] == "continuous") 8 else 100,
    Criterion = if (Variables[i] == "continuous") "BIC" else "D",
    Variables = Variables[i],
    pdf = pdf[i],
    K = K[[i]],
    ymin = ymin[i],
    ymax = ymax[i])

    plot(adultest[[i]], pos = 1, nrow = 1, ncol = 2, what = c("den", "IC"))
    plot(adultest[[i]], pos = 2, nrow = 1, ncol = 2, what = c("den", "IC"))
}

# Best-first feature subset selection.

c <- NULL; rvs <- 1:14; Error <- 1.0

for (i in 1:14) {
  k <- NA

  for (j in rvs) {
    predictive <- RCLSMIX(adultest[c(c, j)],
      P = P,
      Dataset = trainall[, c(c, j)]) 

    CM <- table(train, predictive)

    Accuracy <- (CM[1, 1] + CM[2, 2]) / sum(CM)

    if (1.0 - Accuracy < Error) {
      Error <- 1.0 - Accuracy; k <- j
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

predictive <- RCLSMIX(adultest[c],
  P = P,
  Dataset = testall[, c]) 

CM <- table(test, predictive)

Accuracy <- (CM[1, 1] + CM[2, 2]) / sum(CM)

Error <- 1.0 - Accuracy

# Error on test dataset.

Error
}
}
\keyword{classification}
