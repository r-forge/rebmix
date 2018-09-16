\name{galaxy}
\docType{data}
\alias{galaxy}
\title{Galaxy Dataset}
\description{
The unfilled survey of the Corona Borealis region contains the velocities of 82 galaxies from 6 well separated conic sections of space.
}
\usage{
data("galaxy")
}
\format{
\code{galaxy} is a data frame with 82 cases (rows) and 1 continuous variable (columns) called \code{Velocity}.
}
\source{
K. Roeder. Density estimation with confidence sets exemplified by superclusters and voids
in the galaxies. Journal of American Statistical Association, 85(411):617-624, 1990. \url{http://www.jstor.org/stable/2289993}.
}
\references{
S. Richardson and P. J. Green. On bayesian analysis of mixtures with an unknown number
of components. Journal of the Royal Statistical Society B, 59(4):731-792, 1997. \url{http://www.jstor.org/stable/2985194}.\cr\cr
G. McLachlan and D. Peel. Contribution to the discussion of paper by s. richardson
and p.j. green. Journal of the Royal Statistical Society B, 59(4):779-780, 1997. \url{http://www.jstor.org/stable/2985194}.\cr\cr
M. Stephens. Bayesian analysis of mixture models with an unknown number of components -
an alternative to reversible jump methods. The Annals of Statistics, 28(1):40-74, 2000. \url{http://www.jstor.org/stable/2673981}.
}
\examples{
data("galaxy")

colnames(galaxy)

# Write dataset into tab delimited ASCII file.

write.table(galaxy, file = "galaxy.txt", sep = "\t",
  eol = "\n", row.names = FALSE, col.names = FALSE)

# Write dataset into rda file.

save(galaxy, file = "galaxy.rda")
}
\keyword{datasets}
