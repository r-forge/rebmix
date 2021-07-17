\name{steelplates}
\docType{data}
\alias{steelplates}
\title{Steel Plates Faults Recognition Data}
\description{
These data are the results of an extraction process from images of faults of steel plates. There are seven different faults: Pastry (1), Z_Scratch (2), K_Scratch (3), Stains (4), Dirtiness (5),
Bumps (6), Other faults (7).
}
\usage{
data(steelplates)
}
\format{
\code{steelplates} is a data frame with 1941 cases (rows) and 28 variables (columns) named:
\enumerate{
\item{
\code{X_Minimum} integer.
}
\item{
\code{X_Maximum} integer.
}
\item{
\code{Y_Minimum} integer.
}
\item{
\code{Y_Maximum} integer.
}
\item{
\code{Pixels_Areas} integer.
}
\item{
\code{X_Perimeter} integer.
}
\item{
\code{Y_Perimeter} integer.
}
\item{
\code{Sum_of_Luminosity} integer.
}
\item{
\code{Minimum_of_Luminosity} integer.
}
\item{
\code{Maximum_of_Luminosity} integer.
}
\item{
\code{Length_of_Conveyer} integer.
}
\item{
\code{TypeOfSteel_A300} binary.
}
\item{
\code{TypeOfSteel_A400} binary.
}
\item{
\code{Steel_Plate_Thickness} integer.
}
\item{
\code{Edges_Index} continuous.
}
\item{
\code{Empty_Index} continuous.
}
\item{
\code{Square_Index} continuous.
}
\item{
\code{Outside_X_Index} continuous.
}
\item{
\code{Edges_X_Index} continuous.
}
\item{
\code{Edges_Y_Index} continuous.
}
\item{
\code{Outside_Global_Index} continuous.
}
\item{
\code{LogOfAreas} continuous.
}
\item{
\code{Log_X_Index} continuous.
}
\item{
\code{Log_Y_Index} continuous.
}
\item{
\code{Orientation_Index} continuous.
}
\item{
\code{Luminosity_Index} continuous.
}
\item{
\code{SigmoidOfAreas} continuous.
}
\item{
\code{Class} discrete \code{1}, \code{2}, \code{3}, \code{4}, \code{5}, \code{6} or \code{7}.
}
}
}
\source{
A. Asuncion and D. J. Newman. Uci machine learning repository, 2007. \url{http://archive.ics.uci.edu/ml/}.
}
\references{
M. Buscema, S. Terzi, W. Tastle. A new meta-classifier. Annual Conference of the North American Fuzzy Information Processing Society - NAFIPS, 2010. \doi{10.1109/NAFIPS.2010.5548298}.\cr\cr
M. Buscema. MetaNet*: The theory of independent judges. Substance Use & Misuse. 33(2):439-461, 1998. \doi{10.3109/10826089809115875}.
}
\examples{
\dontrun{
data(steelplates)

# Split dataset into train (75%) and test (25%) subsets.

set.seed(3)

Steelplates <- split(p = 0.75, Dataset = steelplates, class = 28)

# Estimate number of components, component weights and component
# parameters for train subsets.

steelplatesest <- REBMIX(model = "REBMVNORM",
  Dataset = a.train(Steelplates),
  Preprocessing = "histogram",
  cmax = 15,
  Criterion = "BIC")

# Classification.

steelplatescla <- RCLSMIX(model = "RCLSMVNORM",
  x = list(steelplatesest),
  Dataset = a.test(Steelplates),
  Zt = a.Zt(Steelplates))

steelplatescla

summary(steelplatescla)
}
}
\keyword{datasets}
