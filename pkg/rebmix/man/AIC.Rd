\name{AIC}
\alias{AIC}
\alias{AIC3}
\alias{AIC4}
\alias{AICc}
\alias{CAIC}
\title{
Akaike Information Criterion
}
\description{
Returns the Akaike information criterion at \code{pos} for class REBMIX.
}
\usage{
AIC(x, pos = 1, ...)
AIC3(x, pos = 1, ...)
AIC4(x, pos = 1, ...)
AICc(x, pos = 1, ...)
CAIC(x, pos = 1, ...)
}
\arguments{
\item{x}{
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
H. Akaike. A new look at the statistical model identification. IEEE Transactions on Automatic Control, 19(51):716-723, 1974.\cr\cr
A. F. M. Smith and D. J. Spiegelhalter. Bayes factors and choice criteria for linear
models. Journal of the Royal Statistical Society. Series B, 42(2):213-220, 1980. \url{http://www.jstor.org/stable/2984964}.\cr\cr
C. M. Hurvich and C.-L. Tsai. Regression and time series model selection in small samples. Biometrika,
76(2):297-307, 1989. \url{http://www.jstor.org/stable/2336663}.\cr\cr
H. Bozdogan. Model selection and akaike's information criterion (aic): The general theory and its
analytical extensions. Psychometrika, 52(3):345-370, 1987. \url{http://dx.doi.org/10.1007/BF02294361}.
}
\keyword{information criterion}
