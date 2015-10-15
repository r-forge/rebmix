\name{plot.REBMIX}
\alias{plot.REBMIX}
\title{
Plots Univariate or Multivariate REBMIX Output
}
\description{
Takes an object of class \code{REBMIX} and plots empirical and predictive univariate \eqn{d = 1} mixture densities
or pairs of empirical and predictive bivariate marginal mixture densities for \eqn{d > 1}. The empirical
densities follow the histogram, Parzen window or \emph{k}-nearest neighbour approach.
}
\usage{
\method{plot}{REBMIX}(x, pos = 1, what = c("density"),
     nrow = 1, ncol = 1, npts = 200, n = 200, cex = 0.8, fg = "black",
     lty = "solid", lwd = 1, pty = "m", tcl = 0.5, 
     plot.cex = 0.8, plot.pch = 19, contour.drawlabels = FALSE, 
     contour.labcex = 0.8, contour.method = "flattest", 
     contour.nlevels = 12, ...)
}
\arguments{
\item{x}{
an object of class \code{REBMIX}.
}
\item{pos}{
a desired row number in \code{x$summary} to be ploted. The default value is \code{1}.
}
\item{what}{
a character vector, giving the plot types. One of \code{"density"} for probability density function, \code{"marginal"} for marginal probability density function, \code{"IC"} for information criterion depending on numbers of components \eqn{c}, \code{"logL"} for log likelihood, \code{"D"} for total of positive relative deviations, \code{"distribution"} for distribution function or \code{"K"} for information criterion depending on bins \eqn{v} or numbers of nearest neighbours \eqn{k}.
The default value is \code{"density"}.
}
\item{nrow}{
a desired number of rows in which the empirical and predictive densities are to be ploted. The default value is \code{1}.
}
\item{ncol}{
a desired number of columns in which the empirical and predictive densities are to be ploted. The default value is \code{1}.
}
\item{npts}{
a number of points at which the predictive densities are to be ploted. The default value is \code{200}.
}
\item{n}{
a number of observations to be ploted. The default value is \code{200}.
}
\item{cex}{
a numerical value giving the amount by which the plotting text and symbols should be magnified
relative to the default, see also \code{\link[graphics]{par}}. The default value is \code{0.8}.
}
\item{fg}{
a colour used for things like axes and boxes around plots, see also \code{\link[graphics]{par}}. 
The default value is \code{"black"}.
}
\item{lty}{
a line type, see also \code{\link[graphics]{par}}. The default value is \code{"solid"}.
}
\item{lwd}{
a line width, see also \code{\link[graphics]{par}}. The default value is \code{1}.
}
\item{pty}{
a character specifying the type of the plot region to be used. One of \code{"s"} generating a square plotting
region or \code{"m"} generating the maximal plotting region. The default value is \code{"m"}.
}
\item{tcl}{
a length of tick marks as a fraction of the height of a line of the text, see also \code{\link[graphics]{par}}. 
The default value is \code{0.5}.
}
\item{plot.cex}{
a numerical vector giving the amount by which plotting characters and symbols should be
scaled relative to the default. It works as a multiple of \code{\link[graphics]{par}("cex")}. \code{NULL} and \code{NA} are
equivalent to \code{1.0}. Note that this does not affect annotation, see also \code{\link[graphics]{plot.default}}.
The default value is \code{0.8}.
}
\item{plot.pch}{
a vector of plotting characters or symbols, see also \code{\link[graphics]{points}}. The default value is \code{19}.
}
\item{contour.drawlabels}{
logical. The contours are labelled if \code{TRUE}. The default value is \code{FALSE}.
}
\item{contour.labcex}{
\code{cex} for contour labelling. The default value is \code{0.8}. This is an absolute size, not a multiple of
\code{\link[graphics]{par}("cex")}.
}
\item{contour.method}{
a character string specifying where the labels will be located. The possible values
are \code{"simple"}, \code{"edge"} and default \code{"flattest"}, see also \code{\link[graphics]{contour}}.
}
\item{contour.nlevels}{
a number of desired contour levels. The default value is \code{12}.
}
\item{\dots}{
further arguments to \code{\link[graphics]{par}}.
}
}
\value{
\code{par} is returned in an invisible named list. Such a list
can be passed as an argument to \code{\link[graphics]{par}} to restore the parameter values.
}
\references{
C. M. Bishop. Neural Networks for Pattern Recognition. Clarendon Press, Oxford, 1995.
}
\examples{
\dontrun{
devAskNewPage(ask = TRUE)

data("wine")

colnames(wine)

# Remove Cultivar column from wine dataset.

winecolnames <- !(colnames(wine) \%in\% "Cultivar")

wine <- wine[ , winecolnames]

# Determine number of dimensions d and wine dataset size n.

d <- ncol(wine) 
n <- nrow(wine)

# Estimate number of components, component weights and component parameters.

Sturges <- as.integer(1 + log2(n)) # Minimum v follows the Sturges rule.
RootN <- as.integer(2 * n^0.5) # Maximum v follows the RootN rule.

wineest <- REBMIX(Dataset = list(wine = wine), 
  Preprocessing = "Parzen window", 
  Criterion = "ICL-BIC", 
  Variables = rep("continuous", d),
  pdf = rep("normal", d), 
  K = Sturges:RootN)

# Plot finite mixture.

plot(wineest, what = c("density", "marginal", "IC", "logL", "D"),
  nrow = 2, ncol = 2, pty = "s")
}
}
\keyword{plot}

