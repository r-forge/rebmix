\name{bearings}
\docType{data}
\alias{bearings}
\title{Bearings Faults Detection Data}
\description{
These data are the results of the extraction process from the vibrational data of healthy and faulty bearings. Different faults are considered: faultless (1), defect on outer race (2), 
defect on inner race (3) and defect on ball (4). The extracted features are: root mean square (RMS), square root of the amplitude (SRA), kurtosis value (KV), skewness value (SV), 
peak to peak value (PPV), crest factor (CF), impulse factor (IF), margin factor (MF), shape factor (SF), kurtosis factor (KF), frequency centre (FC), root mean square frequency (RMSF) 
and root variance frequency (RVF).
}
\usage{
data(bearings)
}
\format{
\code{bearings} is a data frame with 1906 cases (rows) and 14 variables (columns) named:
\enumerate{
\item{
\code{RMS} continuous.
}
\item{
\code{SRA} continuous.
}
\item{
\code{KV} continuous.
}
\item{
\code{SV} continuous.
}
\item{
\code{PPV} continuous.
}
\item{
\code{CF} continuous.
}
\item{
\code{IF} continuous.
}
\item{
\code{MF} continuous.
}
\item{
\code{SF} continuous.
}
\item{
\code{KF} continuous.
}
\item{
\code{FC} continuous.
}
\item{
\code{RMSF} continuous.
}
\item{
\code{RVF} continuous.
}
\item{
\code{Class} discrete \code{1}, \code{2}, \code{3} or \code{4}.
}
}
}
\source{
Case Western Reserve University Bearing Data Center Website \url{https://csegroups.case.edu/bearingdatacenter/pages/welcome-case-western-reserve-university-bearing-data-center-website}.
}
\references{
B. Panic, J. Klemenc and M. Nagode. Gaussian mixture model based classification revisited: Application to the bearing fault classification. Journal of Mechanical Engineering,
66(4):215-226, 2020. \doi{http://dx.doi.org/10.5545/sv-jme.2020.6563}.
}
\examples{
\dontrun{
data(bearings)

# Split dataset into train (75%) and test (25%) subsets.

set.seed(3)

Bearings <- split(p = 0.75, Dataset = bearings, class = 14)

# Estimate number of components, component weights and component
# parameters for train subsets.

bearingsest <- REBMIX(model = "REBMVNORM",
  Dataset = a.train(Bearings),
  Preprocessing = "histogram",
  cmax = 15,
  Criterion = "BIC")

# Classification.

bearingscla <- RCLSMIX(model = "RCLSMVNORM",
  x = list(bearingsest),
  Dataset = a.test(Bearings),
  Zt = a.Zt(Bearings))

bearingscla

summary(bearingscla)
}
}
\keyword{datasets}
