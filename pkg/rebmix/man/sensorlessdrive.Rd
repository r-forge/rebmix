\name{sensorlessdrive}
\docType{data}
\alias{sensorlessdrive}
\title{Sensorless Drive Faults Detection Data}
\description{
These data are the results of a sensorless drive diagnosis procedure. Features are extracted from the electric current drive signals. The drive has intact and defective components. 
This results in 11 different classes with different conditions. Each condition has been measured several times by 12 different operating conditions, this means by different speeds, 
load moments and load forces. The current signals are measured with a current probe and an oscilloscope on two phases. The original dataset contains 49 features, however, here only 3
are used, that is, features 5, 7 and 11. First class (1) are the healthy drives and the rest are the drives with fault components.
}
\usage{
data(sensorlessdrive)
}
\format{
\code{sensorlessdrive} is a data frame with 58509 cases (rows) and 4 variables (columns) named:
\enumerate{
\item{
\code{V5} continuous.
}
\item{
\code{V7} continuous.
}
\item{
\code{V11} continuous.
}
\item{
\code{Class} discrete \code{1}, \code{2}, \code{3}, \code{4}, \code{5}, \code{6}, \code{7}, \code{8}, \code{9}, \code{10} or \code{11}.
}
}
}
\source{
A. Asuncion and D. J. Newman. Uci machine learning repository, 2007. \url{http://archive.ics.uci.edu/ml/}.
}
\references{
F. Paschke1, C. Bayer, M. Bator, U. Moenks, A. Dicks, O. Enge-Rosenblatt and V. Lohweg. Sensorlose Zustandsueberwachung an Synchronmotoren. 
23. Workshop Computational Intelligence VDI/VDE-Gesellschaft Mess- und Automatisierungstechnik (GMA), 2013.\cr\cr
M. Bator, A. Dicks, U. Moenks and V. Lohweg. Feature extraction and reduction applied to sensorless drive diagnosis. 
22. Workshop Computational Intelligence VDI/VDE-Gesellschaft Mess- und Automatisierungstechnik (GMA), 2012. \doi{10.13140/2.1.2421.5689}. 
}
\examples{
\dontrun{
data(sensorlessdrive)

# Split dataset into train (75%) and test (25%) subsets.

set.seed(3)

Drive <- split(p = 0.75, Dataset = sensorlessdrive, class = 4)

# Estimate number of components, component weights and component
# parameters for train subsets.

driveest <- REBMIX(model = "REBMVNORM",
  Dataset = a.train(Drive),
  Preprocessing = "histogram",
  cmax = 15,
  Criterion = "BIC")

# Classification.

drivecla <- RCLSMIX(model = "RCLSMVNORM",
  x = list(driveest),
  Dataset = a.test(Drive),
  Zt = a.Zt(Drive))

drivecla

summary(drivecla)
}
}
\keyword{datasets}
