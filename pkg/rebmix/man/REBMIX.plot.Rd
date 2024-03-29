\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,REBMIX,missing-method}
\alias{plot,REBMVNORM,missing-method}
\alias{plot,RNGMIX,missing-method}
\alias{plot,RNGMVNORM,missing-method}
\alias{plot,RCLRMIX,missing-method}
\alias{plot,RCLRMVNORM,missing-method}
\alias{plot,RCLSMIX,missing-method}
\alias{plot,RCLSMVNORM,missing-method}
\title{
Plots RNGMIX, REBMIX, RCLRMIX and RCLSMIX Output
}
\description{
Plots true clusters if \code{x} equals \code{"RNGMIX"}. Plots the REBMIX output
depending on \code{what} argument if \code{x} equals \code{"REBMIX"}.
Plots predictive clusters if \code{x} equals \code{"RCLRMIX"}.
Wrongly clustered observations are plotted only if \code{x@Zt} is available.
Plots predictive classes and wrongly classified observations if \code{x} equals \code{"RCLSMIX"}.
}
\usage{
\S4method{plot}{RNGMIX,missing}(x, y, pos = 1, nrow = 1, ncol = 1, cex = 0.8,
     fg = "black", lty = "solid", lwd = 1, pty = "m", tcl = 0.5,
     plot.cex = 0.8, plot.pch = 19, ...)
\S4method{plot}{REBMIX,missing}(x, y, pos = 1, what = c("pdf"),
     nrow = 1, ncol = 1, npts = 200, n = 200, cex = 0.8, fg = "black",
     lty = "solid", lwd = 1, pty = "m", tcl = 0.5,
     plot.cex = 0.8, plot.pch = 19, contour.drawlabels = FALSE,
     contour.labcex = 0.8, contour.method = "flattest",
     contour.nlevels = 12, log = "", ...)
\S4method{plot}{RCLRMIX,missing}(x, y, s = expression(c), nrow = 1, ncol = 1, cex = 0.8,
     fg = "black", lty = "solid", lwd = 1, pty = "m", tcl = 0.5,
     plot.cex = 0.8, plot.pch = 19, ...)
\S4method{plot}{RCLSMIX,missing}(x, y, nrow = 1, ncol = 1, cex = 0.8,
     fg = "black", lty = "solid", lwd = 1, pty = "m", tcl = 0.5,
     plot.cex = 0.8, plot.pch = 19, ...)
## ... and for other signatures
}
\arguments{
\item{x}{
see Methods section below.
}
\item{y}{
currently not used.
}
\item{pos}{
a desired row number in \code{x@summary} to be plotted. The default value is \code{1}.
}
\item{s}{
a desired number of clusters to be plotted. The default value is \code{expression(c)}.
}
\item{what}{
a character vector giving the plot types. One of \code{"pdf"} for probability density function, \code{"marginal pdf"} 
for marginal probability density function, \code{"IC"} for information criterion depending on numbers of components 
\eqn{c}, \code{"logL"} for log likelihood, \code{"D"} for total of positive relative deviations, \code{"K"} for information criterion depending 
on bins \eqn{v} or numbers of nearest neighbours \eqn{k}, \code{"cdf"} for univariate distribution function or
\code{"marginal cdf"} for marginal distribution function. The default value is \code{"pdf"}.
}
\item{nrow}{
a desired number of rows in which the empirical and predictive densities are to be plotted. The default value is \code{1}.
}
\item{ncol}{
a desired number of columns in which the empirical and predictive densities are to be plotted. The default value is \code{1}.
}
\item{npts}{
a number of points at which the predictive densities are to be plotted. The default value is \code{200}.
}
\item{n}{
a number of observations to be plotted. The default value is \code{200}.
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
a character specifying where the labels will be located. The possible values
are \code{"simple"}, \code{"edge"} and default \code{"flattest"}, see also \code{\link[graphics]{contour}}.
}
\item{contour.nlevels}{
a number of desired contour levels. The default value is \code{12}.
}
\item{log}{
a character which contains \code{"x"} if the x axis is to be logarithmic, \code{"y"} if the
y axis is to be logarithmic and \code{"xy"} or \code{"yx"} if both axes are to be logarithmic. The default value is \code{""}.
}
\item{\dots}{
further arguments to \code{\link[graphics]{par}}.
}
}
\value{
Returns (invisibly) a list containing graphical parameters \code{par}. Such a list can be passed as an argument to \code{\link[graphics]{par}} to restore the parameter values.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "RNGMIX", y = "missing")}}{an object of class \code{RNGMIX}.}
\item{\code{signature(x = "RNGMVNORM", y = "missing")}}{an object of class \code{RNGMVNORM}.}
\item{\code{signature(x = "REBMIX", y = "missing")}}{an object of class \code{REBMIX}.}
\item{\code{signature(x = "REBMVNORM", y = "missing")}}{an object of class \code{REBMVNORM}.}
\item{\code{signature(x = "RCLRMIX", y = "missing")}}{an object of class \code{RCLRMIX}.}
\item{\code{signature(x = "RCLRMVNORM", y = "missing")}}{an object of class \code{RCLRMVNORM}.}
\item{\code{signature(x = "RCLSMIX", y = "missing")}}{an object of class \code{RCLSMIX}.}
\item{\code{signature(x = "RCLSMVNORM", y = "missing")}}{an object of class \code{RCLSMVNORM}.}
}
}
\author{Marko Nagode}
\references{
C. M. Bishop. Neural Networks for Pattern Recognition. Clarendon Press, Oxford, 1995.
}
\examples{
\dontrun{
devAskNewPage(ask = TRUE)

data(wine)

colnames(wine)

# Remove Cultivar column from wine dataset.

winecolnames <- !(colnames(wine) %in% "Cultivar")

wine <- wine[, winecolnames]

# Determine number of dimensions d and wine dataset size n.

d <- ncol(wine)
n <- nrow(wine)

wineest <- REBMIX(model = "REBMVNORM",
  Dataset = list(wine = wine),
  Preprocessing = "kernel density estimation",
  Criterion = "ICL-BIC",
  EMcontrol = new("EM.Control", strategy = "best"))

# Plot finite mixture.

plot(wineest, what = c("pdf", "IC", "logL", "D"),
  nrow = 2, ncol = 2, pty = "s")
}
}
\keyword{plot}

