\name{RCLRMIX-class}
\docType{class}
\alias{RCLRMIX-class}
\alias{RCLRMVNORM-class}
\title{Class \code{"RCLRMIX"}}
\description{Object of class \code{RCLRMIX}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RCLRMIX", ...)}.
Accessor methods for the slots are \code{a.pos(x = NULL)}, \code{a.Zt(x = NULL)}, \code{a.Zp(x = NULL, s = expression(c))}, \code{a.c(x = NULL)},
\code{a.p(x = NULL, s = expression(c))}, \code{a.pi(x = NULL, s = expression(c))}, \cr
\code{a.P(x = NULL, s = expression(c))}, \code{a.tau(x = NULL, s = expression(c))},
\code{a.prob(x = NULL)}, \code{a.from(x = NULL)}, \code{a.to(x = NULL)}, \code{a.EN(x = NULL)} and \code{a.ED(x = NULL)}, where \code{x} stands for an object of class \code{RCLRMIX} and \code{s}
a desired number of clusters for which the slot is calculated.
}
\section{Slots}{
\describe{
\item{\code{x}:}{
an object of class \code{REBMIX}.
}
\item{\code{pos}:}{
a desired row number in \code{x@summary} to be clustered. The default value is \code{1}.
}
\item{\code{Zt}:}{
a factor of true cluster membership.
}
\item{\code{Zp}:}{
a factor of predictive cluster membership.
}
\item{\code{c}:}{
number of clusters.
}
\item{\code{p}:}{
a vector of length \eqn{c} containing prior probabilities of cluster memberships \eqn{p_{l}} summing to 1. The value is returned only if all variables in slot \code{x} follow either binomial or Dirac parametric families. The default value is \code{numeric()}.
}
\item{\code{pi}:}{
a list of length \eqn{d} of matrices of size \eqn{c \times K_{i}} containing cluster conditional probabilities \eqn{\pi_{ilk}}. Let \eqn{\pi_{ilk}}
denote the cluster conditional probability that an observation in cluster \eqn{l = 1, \ldots, c} produces the \eqn{k}th outcome on the \eqn{i}th variable.
Suppose we observe \eqn{i = 1, \ldots, d} polytomous categorical variables (the manifest variables), each of which contains \eqn{K_{i}} possible outcomes for observations \eqn{j = 1, \ldots, n}.
A manifest variable is a variable that can be measured or observed directly. It must be coded as whole number starting at zero for the first outcome and increasing to the possible number of outcomes minus one.
It is presumed here that all variables are statistically independentand within clusters and that \eqn{\bm{y}_{1}, \ldots, \bm{y}_{n}}
stands for an observed \eqn{d} dimensional dataset of size \eqn{n} of vector observations \eqn{\bm{y}_{j} = (y_{1j}, \ldots, y_{ij}, \ldots, y_{dj})^\top}.
The value is returned only if all variables in slot \code{x} follow either binomial or Dirac parametric families. The default value is \code{list()}.
}
\item{\code{P}:}{
a data frame containing true \eqn{N_{\mathrm{t}}(\bm{y}_{\tilde{\jmath}})} and predictive \eqn{N_{\mathrm{p}}(\bm{y}_{\tilde{\jmath}})} frequencies calculated for unique \eqn{\bm{y}_{\tilde{\jmath}} \in \{ \bm{y}_{1}, \ldots, \bm{y}_{n} \}}, where \eqn{\tilde{\jmath} = 1, \ldots, \tilde{n}$ and $\tilde{n} \leq n}.
}
\item{\code{tau}:}{
a matrix of size \eqn{n \times c} containing conditional probabilities \eqn{\tau_{jl}} that observations \eqn{\bm{y}_{1}, \ldots, \bm{y}_{n}} arise from clusters \eqn{1, \ldots, c}.
}
\item{\code{prob}:}{
a vector of length \eqn{c} containing probabilities of correct clustering for \eqn{s = 1, \ldots, c}.
}
\item{\code{from}:}{
a vector of length \eqn{c - 1} containing clusters merged to \code{to} clusters.
}
\item{\code{to}:}{
a vector of length \eqn{c - 1} containing clusters originating from \code{from} clusters.
}
\item{\code{EN}:}{
a vector of length \eqn{c - 1} containing entropies for combined clusters.
}
\item{\code{ED}:}{
a vector of length \eqn{c - 1} containing decrease of entropies for combined clusters.
}
}
}
\author{Marko Nagode}
\references{
J. P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo. Combining mixture components for clustering.
Journal of Computational and Graphical Statistics, 19(2):332-353, 2010. \url{http://dx.doi.org/10.1198/jcgs.2010.08111}
}
\examples{
devAskNewPage(ask = TRUE)

# Generate normal dataset.

n <- c(500, 200, 400)

Theta <- new("RNGMVNORM.Theta", c = 3, d = 2)

a.theta1(Theta, 1) <- c(3, 10)
a.theta1(Theta, 2) <- c(8, 6)
a.theta1(Theta, 3) <- c(12, 11)
a.theta2(Theta, 1) <- c(3, 0.3, 0.3, 2)
a.theta2(Theta, 2) <- c(5.7, -2.3, -2.3, 3.5)
a.theta2(Theta, 3) <- c(2, 1, 1, 2)

normal <- RNGMIX(model = "RNGMVNORM", Dataset.name = "normal_1", n = n, Theta = a.Theta(Theta))

# Estimate number of components, component weights and component parameters.

normalest <- REBMIX(model = "REBMVNORM",
  Dataset = a.Dataset(normal),
  Preprocessing = "histogram",
  cmax = 6,
  Criterion = "BIC",
  pdf = rep("normal", 2))

summary(normalest)

# Plot finite mixture.

plot(normalest)

# Cluster dataset.

normalclu <- RCLRMIX(model = "RCLRMVNORM", x = normalest, Zt = a.Zt(normal))

# Plot clusters.

plot(normalclu)

summary(normalclu)
}
\keyword{classes}
