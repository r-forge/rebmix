\name{adult}
\docType{data}
\alias{adult}
\title{Adult Dataset}
\description{
The \code{adult} dataset containing 48842 instances with 16 continuous, binary and discrete variables was extracted from the census bureau database \url{http://www.census.gov/}. Extraction was done by Barry Becker from the 1994 census bureau database.
}
\usage{
adult
}
\format{
\code{adult} is a data frame with 48842 cases (rows) and 16 variables (columns) named:
\enumerate{
\item{
\code{Type} binary \code{train} or \code{test}.
}
\item{
\code{Age} continuous.
}
\item{
\code{Workclass} one of the 8 discrete values
\code{private},
\code{self-emp-not-inc},
\code{self-emp-inc},
\code{federal-gov},
\code{local-gov},
\code{state-gov},
\code{without-pay} or
\code{never-worked}.
}
\item{
\code{Fnlwgt} stands for continuous final weight.
}
\item{
\code{Education} one of the 16 discrete values
\code{bachelors},
\code{some-college},
\code{11th},
\code{hs-grad},
\code{prof-school},
\code{assoc-acdm},
\code{assoc-voc},
\code{9th},
\code{7th-8th},
\code{12th},
\code{masters},
\code{1st-4th},
\code{10th},
\code{doctorate},
\code{5th-6th} or
\code{preschool}.
}
\item{
\code{Education.Num} continuous.
}
\item{
\code{Marital.Status} one of the 7 discrete values
\code{married-civ-spouse},
\code{divorced},
\code{never-married},
\code{separated},
\code{widowed},
\code{married-spouse-absent} or
\code{married-af-spouse}.
}
\item{
\code{Occupation} one of the 14 discrete values
\code{tech-support},
\code{craft-repair},
\code{other-service},
\code{sales},
\code{exec-managerial},
\code{prof-specialty},
\code{handlers-cleaners},
\code{machine-op-inspct},
\code{adm-clerical},
\code{farming-fishing},
\code{transport-moving},
\code{priv-house-serv},
\code{protective-serv} or
\code{armed-forces}.
}
\item{
\code{Relationship} one of the 6 discrete values
\code{wife},
\code{own-child},
\code{husband},
\code{not-in-family},
\code{other-relative} or
\code{unmarried}.
}
\item{
\code{Race} one of the 5 discrete values
\code{white},
\code{asian-pac-islander},
\code{amer-indian-eskimo},
\code{other} or
\code{black}.
}
\item{
\code{Sex} binary \code{female} or \code{male}.
}
\item{
\code{Capital.Gain} continuous.
}
\item{
\code{Capital.Loss} continuous.
}
\item{
\code{Hours.Per.Week} continuous.
}
\item{
\code{Native.Country} one of the 41 discrete values
\code{united-states},
\code{cambodia},
\code{england},
\code{puerto-rico},
\code{canada},
\code{germany},
\code{outlying-us(guam-usvi-etc)},
\code{india},
\code{japan},
\code{greece},
\code{south},
\code{china},
\code{cuba},
\code{iran},
\code{honduras},
\code{philippines},
\code{italy},
\code{poland},
\code{jamaica},
\code{vietnam},
\code{mexico},
\code{portugal},
\code{ireland},
\code{france},
\code{dominican-republic},
\code{laos},
\code{ecuador},
\code{taiwan},
\code{haiti},
\code{columbia},
\code{hungary},
\code{guatemala},
\code{nicaragua},
\code{scotland},
\code{thailand},
\code{yugoslavia},
\code{el-salvador},
\code{trinadad&tobago},
\code{peru},
\code{hong} or
\code{holand-netherlands}.
}
\item{
\code{Income} binary \code{<=50k} or \code{>50k}.
}
}
}
\source{
Frank A, Asuncion A (2010). UCI Machine Learning Repository. \url{http://archive.ics.uci.edu/ml}.
}
\references{
Frank A, Asuncion A (2010). UCI Machine Learning Repository. \url{http://archive.ics.uci.edu/ml}.
}
\examples{
data("adult")

## Find complete cases.

adult <- adult[complete.cases(adult), ]

## Map metric attributes.

adult[["Capital.Loss"]] <- ordered(cut(adult[["Capital.Loss"]], 2000))
adult[["Capital.Gain"]] <- ordered(cut(adult[["Capital.Gain"]], 2000))

## Show level attributes for binary and discrete variables.

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

## Replace levels with numbers.

adult <- as.data.frame(data.matrix(adult))

## Levels should start with 0 for discrete distributions except for the 
## Dirac distribution.

f <- c("Type", "Workclass", "Education", "Marital.Status", "Occupation", 
  "Relationship", "Race", "Sex", "Native.Country", "Income")

adult[, f] <- adult[, f] - 1

## Split adult dataset into two train subsets for the two Incomes
## and remove Type and Income columns.

trainle50k <- subset(adult, subset = (Type == 1) & (Income == 0), 
  select = c(-Type, -Income))
traingt50k <- subset(adult, subset = (Type == 1) & (Income == 1), 
  select = c(-Type, -Income))

trainall <- subset(adult, subset = Type == 1, select = c(-Type, -Income))

train <- as.factor(subset(adult, subset = Type == 1, select = c(Income))[, 1])

## Extract test dataset form adult dataset and remove Type 
## and Income columns.

testle50k <- subset(adult, subset = (Type == 0) & (Income == 0), 
  select = c(-Type, -Income))
testgt50k <- subset(adult, subset = (Type == 0) & (Income == 1), 
  select = c(-Type, -Income))

testall <- subset(adult, subset = Type == 0, select = c(-Type, -Income))

test <- as.factor(subset(adult, subset = Type == 0, select = c(Income))[, 1])

save(trainall, file = "trainall.rda")
save(testall, file = "testall.rda")
}
\keyword{datasets}
