\name{labelmoments-methods}
\alias{labelmoments}
\alias{labelmoments-methods}
\alias{labelmoments,array-method}
\title{
Label Image Moments
}
\description{
Returns the list with the data frame \code{Mij} containing the cluster levels \eqn{l}, the numbers of pixels \eqn{n} and the cluster moments
\eqn{\bm{M} = (M_{\mathrm{10}}, M_{\mathrm{01}}, M_{\mathrm{11}})^{\top}} for 2D images or the data frame \code{Mijk} containing the cluster levels \eqn{l}, the numbers of voxels \eqn{n} and the cluster moments \eqn{\bm{M} = (M_{\mathrm{100}}, M_{\mathrm{010}}, M_{\mathrm{001}}, M_{\mathrm{111}})^{\top}} 
for 3D images and the adjacency matrix \code{A} of size \eqn{c_{\mathrm{max}} \times c_{\mathrm{max}}}. It may have some \code{NA} rows and columns. To calculate the adjacency matrix \eqn{A(i,j) = \exp{(-\|\bm{M}_{i} - \bm{M}_{j}\|^2 / 2 \sigma^2)}}, the raw cluster moments are first converted into z-scores.
}
\usage{
\S4method{labelmoments}{array}(Zp = array(), cmax = integer(), Sigma = 1.0, ...)
## ... and for other signatures
}
\arguments{
\item{Zp}{
a 2D array of size \eqn{width \times height} or a 3D array of size \eqn{width \times height \times depth} containing the predictive cluster membership \eqn{\bm{\Omega}_{l}}, where \eqn{l \in \{ 0, 1, \ldots, c\}}. The cluster \eqn{l = 0} may contain the pixels for 2D images or the voxels for 3D images, which are ignored by this method.
}
\item{cmax}{
maximum number of clusters \eqn{c_{\mathrm{max}} \geq c}.
}
\item{Sigma}{
scale parameter \eqn{\sigma} for the Gaussian kernel function \eqn{\exp{(-\|\bm{M}_{i} - \bm{M}_{j}\|^2 / 2 \sigma^2)}}. The default value is \code{1.0}.
}
\item{\dots}{
currently not used.
}
}
\section{Methods}{
\describe{
\item{\code{signature(Zp = "array")}}{an array.}
}
}
\author{Marko Nagode, Branislav Panic}
\references{
A. Ng, M. Jordan and Y. Weiss. On spectral clustering: Analysis and an algorithm. Advances in Neural Information Processing Systems 14 (NIPS 2001).
}
\examples{
Zp <- matrix(rep(0, 100), nrow = 10, ncol = 10)

Zp[2, 2:4] <- 1; 
Zp[2:4, 5] <- 2; 
Zp[8, 7:10] <- 3; 
Zp[9, 6] <- 4; Zp[10, 5] <- 4
Zp[10, 1:4] <- 5
Zp[6:9, 1] <- 6

labelmoments <- labelmoments(Zp, cmax = 6, Sigma = 1.0)

set.seed(12)

mergelabels <- mergelabels(list(labelmoments$A), w = 1.0, k = 2, nstart = 3)

Zp

mergelabels
}
\keyword{auxiliary}
