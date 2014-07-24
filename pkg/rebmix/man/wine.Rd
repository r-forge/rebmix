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
wine
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
\code{Cultivar} continuous.
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
data("wine")

colnames(wine)

## Split wine dataset into three subsets for three Cultivars
## and remove Cultivar column from datasets.

winecolnames <- !(colnames(wine) \%in\% "Cultivar")

wine1 <- wine[wine$Cultivar == 1, winecolnames]
wine2 <- wine[wine$Cultivar == 2, winecolnames]
wine3 <- wine[wine$Cultivar == 3, winecolnames]

wine <- wine[ , winecolnames]

## Write datasets without Cultivar information into tab 
## delimited ASCII files.

write.table(wine, file = "wine.txt", sep = "\t", 
  eol = "\n", row.names = FALSE, col.names = FALSE)
write.table(wine1, file = "wine1.txt", sep = "\t", 
  eol = "\n", row.names = FALSE, col.names = FALSE)
write.table(wine2, file = "wine2.txt", sep = "\t", 
  eol = "\n", row.names = FALSE, col.names = FALSE)
write.table(wine3, file = "wine3.txt", sep = "\t", 
  eol = "\n", row.names = FALSE, col.names = FALSE)
}
\keyword{datasets}
