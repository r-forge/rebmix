\name{AIC-methods}
\docType{methods}
\alias{AIC}
\alias{AIC3}
\alias{AIC4}
\alias{AICc}
\alias{CAIC}
\alias{AIC-methods}
\alias{AIC,REBMIX-method}
\alias{AIC,REBMVNORM-method}
\alias{AIC3-methods}
\alias{AIC3,REBMIX-method}
\alias{AIC3,REBMVNORM-method}
\alias{AIC4-methods}
\alias{AIC4,REBMIX-method}
\alias{AIC4,REBMVNORM-method}
\alias{AICc-methods}
\alias{AICc,REBMIX-method}
\alias{AICc,REBMVNORM-method}
\alias{CAIC-methods}
\alias{CAIC,REBMIX-method}
\alias{CAIC,REBMVNORM-method}
\title{
Akaike Information Criterion
}
\description{
Returns the Akaike information criterion at \code{pos}.
}
\usage{
\S4method{AIC}{REBMIX}(x = NULL, pos = 1, ...)
\S4method{AIC3}{REBMIX}(x = NULL, pos = 1, ...)
\S4method{AIC4}{REBMIX}(x = NULL, pos = 1, ...)
\S4method{AICc}{REBMIX}(x = NULL, pos = 1, ...)
\S4method{CAIC}{REBMIX}(x = NULL, pos = 1, ...)
## ... and for other signatures
}
\arguments{
\item{x}{
see Methods section below.
}
\item{pos}{
a desired row number in \code{x@summary} for which the information criterion is calculated. The default value is \code{1}.
}
\item{\dots}{
currently not used.
}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "REBMIX")}}{an object of class \code{REBMIX}.}
\item{\code{signature(x = "REBMVNORM")}}{an object of class \code{REBMVNORM}.}
}
}
\author{Marko Nagode}
\references{
H. Akaike. A new look at the statistical model identification. IEEE Transactions on Automatic Control, 19(51):716-723, 1974.\cr\cr
A. F. M. Smith and D. J. Spiegelhalter. Bayes factors and choice criteria for linear
models. Journal of the Royal Statistical Society. Series B, 42(2):213-220, 1980. https://www.jstor.org/stable/2984964.\cr\cr
H. Bozdogan. Model selection and akaike's information criterion (aic): The general theory and its
analytical extensions. Psychometrika, 52(3):345-370, 1987. \doi{10.1007/BF02294361}.\cr\cr
C. M. Hurvich and C.-L. Tsai. Regression and time series model selection in small samples. Biometrika,
76(2):297-307, 1989. https://www.jstor.org/stable/2336663.
}
\keyword{information criterion}
