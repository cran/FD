\encoding{latin1}
\name{mahaldis}
\alias{mahaldis}
\title{Mahalanobis Distance}
\description{
\code{mahaldis} measures the pairwise Mahalanobis (1936) distances between individual objects.
}
\usage{
mahaldis(x)
}
\arguments{
  \item{x}{matrix containing the variables. \code{\link{NA}}s are not tolerated.}
}
\details{
\code{mahaldis} computes the Mahalanobis (1936) distances between individual objects. The Mahalanobis distance takes into account correlations among variables and does not depend on the scales of the variables.

\code{mahaldis} builds on the fact that type-II principal component analysis (PCA) preserves the Mahalanobis distance among objects (Legendre and Legendre 2012). Therefore, \code{mahaldis} first performs a type-II PCA on standardized variables, and then computes the Euclidean distances among (repositioned) objects whose positions are given in the matrix \eqn{\mathbf{G}}{G}. This is equivalent to the Mahalanobis distances in the space of the original variables (Legendre and Legendre 2012). 
}

\value{
an object of class \code{\link{dist}}.
}
\references{
Legendre, P. and L. Legendre (2012) \emph{Numerical Ecology}. 3nd English edition. Amsterdam: Elsevier.
}
\author{
Pierre Legendre \email{pierre.legendre@umontreal.ca}

\url{http://www.bio.umontreal.ca/legendre/indexEn.html}

Ported to \pkg{FD} by Etienne \enc{Laliberté}{Laliberte}.
}
\seealso{\code{\link[stats]{mahalanobis}} computes the Mahalanobis distances among groups of objects, not individual objects.}
\examples{
mat <- matrix(rnorm(100), 50, 20)

ex1 <- mahaldis(mat)

# check attributes
attributes(ex1)
}
\keyword{multivariate}
