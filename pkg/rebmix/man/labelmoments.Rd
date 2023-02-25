\name{labelmoments-methods}
\alias{labelmoments}
\alias{labelmoments-methods}
\alias{labelmoments,array-method}
\title{
Label Image Moments
}
\description{
Returns the list with the data frame \code{Mij} containing the cluster levels \eqn{l}, the numbers of pixels \eqn{n} and the rough cluster image moments
\eqn{\bm{M} = (M_{\mathrm{10}}, M_{\mathrm{01}}, M_{\mathrm{11}})^{\top}} for 2D images or the data frame \code{Mijk} containing the cluster levels \eqn{l}, the numbers of voxels \eqn{n} and the rough cluster image moments \eqn{\bm{M} = (M_{\mathrm{100}}, M_{\mathrm{010}}, M_{\mathrm{001}}, M_{\mathrm{111}})^{\top}} 
for 3D images and the normalised affinity matrix \code{L}. The normalised affinity matrix \eqn{L = D^{-1/2} A D^{-1/2}} depends on the affinity matrix \eqn{A_{i,j} = \exp{(-\|\bm{M}_{i} - \bm{M}_{j}\|^2 / 2 \sigma^2)}}
and the degree matrix \eqn{D_{i,j} = \sum_{j} A_{i,j}}. To calculate the affinity matrix, the rough cluster moments are first converted into z-scores.
}
\usage{
\S4method{labelmoments}{array}(Zp = array(), ...)
## ... and for other signatures
}
\arguments{
\item{Zp}{
a 2D array of size \eqn{width \times height} or a 3D array of size \eqn{width \times height \times depth} containing the predictive cluster membership \eqn{\bm{\Omega}_{l}}, where \eqn{l \in \{ 0, 1, \ldots, c\}}. The cluster \eqn{l = 0} may contain the pixels for 2D images or the voxels for 3D images, which are ignored by this method.
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
i <- 1
}
\keyword{auxiliary}
