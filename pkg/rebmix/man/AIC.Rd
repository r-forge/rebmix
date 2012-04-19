\name{AIC}
\alias{AIC.REBMIX}
\alias{AIC3}
\alias{AIC3.default}
\alias{AIC3.REBMIX}
\alias{AIC4}
\alias{AIC4.default}
\alias{AIC4.REBMIX}
\alias{AICc}
\alias{AICc.default}
\alias{AICc.REBMIX}
\alias{CAIC}
\alias{CAIC.default}
\alias{CAIC.REBMIX}
\title{
Akaike Information Criterion
}
\description{
Returns the Akaike information criterion at \code{pos} for class REBMIX.
}
\usage{
\method{AIC}{REBMIX}(object, pos = 1, ...)
\method{AIC3}{REBMIX}(object, pos = 1, ...)
\method{AIC4}{REBMIX}(object, pos = 1, ...)
\method{AICc}{REBMIX}(object, pos = 1, ...)
\method{CAIC}{REBMIX}(object, pos = 1, ...)
}
\arguments{
\item{object}{
an object of class \code{REBMIX}.
}
\item{pos}{
a desired row number in \code{summary} for which the information criterion is calculated. The default value is \code{1}.
}
\item{\dots}{
currently not used.
}
}
\references{
Akaike H (1974). A New Look at the Statistical Model Identification. \emph{IEEE Transactions
on Automatic Control}, 19, 716-723.\cr\cr
Smith AFM, Spiegelhalter DJ (1980). Bayes Factors and Choice Criteria for Linear Models. \emph{Journal 
of the Royal Statistical Society B}, 42, 213-220.\cr\cr
Hurvich CM, Tsai CL (1989). Regression and Time Series Model Selection in Small Samples. \emph{Biometrika}, 76, 297-307.\cr\cr
Bozdogan H (1987). Model Selection and Akaike's Information Criterion (AIC): The General Theory
and its Analytical Extensions. \emph{Psychometrika}, 52, 345-370.
}
\keyword{information criterion}
