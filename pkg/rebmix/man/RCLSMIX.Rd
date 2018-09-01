\name{RCLSMIX-methods}
\docType{methods}
\alias{RCLSMIX}
\alias{RCLSMIX-methods}
\alias{RCLSMIX,RCLSMIX-method}
\alias{RCLSMIX,RCLSMVNORM-method}
\alias{show,RCLSMIX-method}
\alias{show,RCLSMVNORM-method}
\alias{summary,RCLSMIX-method}
\alias{summary,RCLSMVNORM-method}
\alias{g.o}
\alias{g.s}
\alias{g.ntrain}
\alias{g.ntest}
\alias{g.CM}
\alias{g.Accuracy}
\alias{g.Error}
\alias{g.Precision}
\alias{g.Sensitivity}
\alias{g.Specificity}
\alias{g.Chunks}
\alias{g.o,RCLSMIX-method}
\alias{g.Dataset,RCLSMIX-method}
\alias{g.s,RCLSMIX-method}
\alias{g.ntrain,RCLSMIX-method}
\alias{g.P,RCLSMIX-method}
\alias{g.ntest,RCLSMIX-method}
\alias{g.Zt,RCLSMIX-method}
\alias{g.Zp,RCLSMIX-method}
\alias{g.CM,RCLSMIX-method}
\alias{g.Accuracy,RCLSMIX-method}
\alias{g.Error,RCLSMIX-method}
\alias{g.Precision,RCLSMIX-method}
\alias{g.Sensitivity,RCLSMIX-method}
\alias{g.Specificity,RCLSMIX-method}
\alias{g.Chunks,RCLSMIX-method}
\alias{g.o,RCLSMVNORM-method}
\alias{g.Dataset,RCLSMVNORM-method}
\alias{g.s,RCLSMVNORM-method}
\alias{g.ntrain,RCLSMVNORM-method}
\alias{g.P,RCLSMVNORM-method}
\alias{g.ntest,RCLSMVNORM-method}
\alias{g.Zt,RCLSMVNORM-method}
\alias{g.Zp,RCLSMVNORM-method}
\alias{g.CM,RCLSMVNORM-method}
\alias{g.Accuracy,RCLSMVNORM-method}
\alias{g.Error,RCLSMVNORM-method}
\alias{g.Precision,RCLSMVNORM-method}
\alias{g.Sensitivity,RCLSMVNORM-method}
\alias{g.Specificity,RCLSMVNORM-method}
\alias{g.Chunks,RCLSMVNORM-method}
\title{
Predicts Class Membership Based Upon a Model Trained by REBMIX
}
\description{
Returns as default the RCLSMIX algorithm output for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson, Dirac or von Mises component densities. If \code{model} equals \code{"RCLSMVNORM"} output for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices is returned.
}
\usage{
\S4method{RCLSMIX}{RCLSMIX}(model = "RCLSMIX", x = list(), Dataset = data.frame(),
        Zt = factor(), ...)
## ... and for other signatures
\S4method{summary}{RCLSMIX}(object, ...)
## ... and for other signatures     
}
\arguments{
\item{model}{
see Methods section below.
}
\item{x}{
a list of objects of class \code{REBMIX} of length \eqn{o} obtained by running \code{\link[rebmix]{REBMIX}} on \eqn{g = 1, \ldots, s} train datasets \eqn{Y_{\mathrm{train}g}} all of length \eqn{n_{\mathrm{train}g}}.
For the train datasets the corresponding class membership \eqn{\bm{\Omega}_{g}} is known. This yields
\eqn{n_{\mathrm{train}} = \sum_{g = 1}^{s} n_{\mathrm{train}g}}, while \eqn{Y_{\mathrm{train}q} \cap Y_{\mathrm{train}g} = \emptyset} for all \eqn{q \neq g}.
Each object in the list corresponds to one chunk, e.g., \eqn{(y_{1j}, y_{3j})^{\top}}. The default value is \code{list()}.
}
\item{Dataset}{
a data frame containing test dataset \eqn{Y_{\mathrm{test}}} of length \eqn{n_{\mathrm{test}}}. For the test dataset the corresponding class membership \eqn{\bm{\Omega}_{g}} is not known.
The default value is \code{data.frame()}.
}
\item{Zt}{
a factor of true class membership \eqn{\bm{\Omega}_{g}} for the test dataset. The default value is \code{factor()}.
}
\item{object}{
see Methods section below.
}
\item{\dots}{
currently not used; additional arguments affecting the summary produced.
}
}
\value{
Returns an object of class \code{RCLSMIX} or \code{RCLSMVNORM}.
}
\section{Methods}{
\describe{
\item{\code{signature(model = "RCLSMIX")}}{a character giving the default class name \code{"RCLSMIX"} for mixtures of conditionally independent normal, lognormal, Weibull, gamma, binomial, Poisson, Dirac or von Mises component densities.}
\item{\code{signature(model = "RCLSMVNORM")}}{a character giving the class name \code{"RCLSMVNORM"} for mixtures of multivariate normal component densities with unrestricted variance-covariance matrices.}
\item{\code{signature(object = "RCLSMIX")}}{an object of class \code{RCLSMIX}.}
\item{\code{signature(object = "RCLSMVNORM")}}{an object of class \code{RCLSMVNORM}.}
}
}
\author{Marko Nagode}
\references{
R. O. Duda and P. E. Hart. Pattern Classification and Scene Analysis. John Wiley & Sons, New
York, 1973.
}
\examples{
\dontrun{
devAskNewPage(ask = TRUE)

data("adult")

# Find complete cases.

adult <- adult[complete.cases(adult),]

# Replace levels with numbers.

adult <- as.data.frame(data.matrix(adult))

# Find numbers of levels.

cmax <- unlist(lapply(apply(adult[, c(-1, -16)], 2, unique), length))

cmax

# Split adult dataset into train and test subsets for two Incomes
# and remove Type and Income columns.

Adult <- split(p = list(type = 1, train = 2, test = 1),
  Dataset = adult, class = 16)

# Estimate number of components, component weights and component parameters
# for the set of chunks 1:14.

adultest <- list()

for (i in 1:14) {
  adultest[[i]] <- REBMIX(Dataset = chunk(Adult, i)@train,
    Preprocessing = "histogram",
    cmax = min(120, cmax[i]),
    Criterion = "BIC",
    pdf = "Dirac",
    K = 1)
}

# Class membership prediction based upon the best first search algorithm.

adultcla <- BFSMIX(x = adultest, 
  Dataset = Adult@test, 
  Zt = Adult@Zt)

adultcla

summary(adultcla)

# Plot selected chunks.

plot(adultcla, nrow = 5, ncol = 2)
}
}
\keyword{classification}
