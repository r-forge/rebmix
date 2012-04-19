\name{galaxy}
\docType{data}
\alias{galaxy}
\title{Galaxy Dataset}
\description{
The unfilled survey of the Corona Borealis region contains the velocities of 82 galaxies from 6 well separated conic sections of space.
}
\usage{
galaxy
}
\format{
\code{galaxy} is a data frame with 82 cases (rows) and 1 continuous variable (columns) called \code{Velocity}.
}
\source{
Roeder K (1990). Density Estimation with Confidence Sets Exemplified by Superclusters
and Voids in the Galaxies. \emph{Journal of American Statistical Association}, 85, 617-624.
}
\references{
Richardson S, Green PJ (1997). On Bayesian Analysis of Mixtures with an Unknown Number
of Components. \emph{Journal of the Royal Statistical Society B}, 59, 731-792.\cr\cr
McLachlan GJ, Peel D (1997). Contribution to the Discussion of Paper by S. Richardson
and P.J. Green. \emph{Journal of the Royal Statistical Society B}, 59, 779-780.\cr\cr
Stephens M (2000). Bayesian Analysis of Mixture Models with an Unknown Number of
Components - An Alternative to Reversible Jump Methods. \emph{The Annals of Statistics}, 28, 40-74.
}
\examples{
data("galaxy")

colnames(galaxy)

## Write dataset into tab delimited ASCII file.

write.table(galaxy, file = "galaxy.txt", sep = "\t", 
  eol = "\n", row.names = FALSE, col.names = FALSE)
  
## Write dataset into rda file.

save(galaxy, file = "galaxy.rda")
}
\keyword{datasets}
