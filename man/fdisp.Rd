\encoding{latin1}
\name{fdisp}
\alias{fdisp}
\concept{multivariate dispersion}
\concept{functional dispersion}
\concept{functional diversity}
\title{Functional Dispersion}
\description{
\code{fdisp} measures the functional dispersion (FDis) of a set of communities, as described by \enc{Laliberté}{Laliberte} and Legendre (submitted).
}
\usage{
fdisp(d, a, tol = 1e-07)
}
\arguments{
  \item{d}{a species-by- species distance matrix computed from functional traits, such as that returned by \code{\link{dist}} or \code{\link{gowdis}}. \code{\link{NA}s} are not allowed.}
  \item{a}{matrix containing the abundances of the species in \code{d} (or presence-absence, i.e. 0 or 1). Rows are sites and species are columns. Can be missing, in which case \code{fdisp} assumes that there is only one community with equal abundances of all species. \code{\link{NA}s} will be replaced by 0.
The number of species (columns) in \code{a} must match the number of species in \code{d}. In addition, the species labels in \code{a} and \code{d} must be identical and in the same order.}
  \item{tol}{tolerance threshold to test whether the distance matrix is Euclidean : an eigenvalue is considered positive if it is larger than -tol*\eqn{{\lambda}_1}{lambda1}, where \eqn{{\lambda}_1}{lambda1} is the largest eigenvalue.}
}
\details{
\code{fdisp} computes, for a set of communities, the average distance of individual objects (species) in PCoA space from any distance or dissimilarity measure, as described by Anderson (2006). The average distance to the centroid is a measure of multivariate dispersion and as been suggested as an index of beta diversity (Anderson et al. 2006). However, in \code{fdisp} both the centroid and the average distance to this centroid can be weighted by individual objects. In other words, \code{fdisp} returns the weighted average distance to the weighted centroid. This was suggested so that multivariate dispersion could be used as a multidimensional functional diversity (FD) index that can be weighted by species abundances. This FD index has been called functional dispersion (FDis) and is described by \enc{Laliberté}{Laliberte} and Legendre (submitted).

In sum, FDis can account for relative abundances by shifting the position of the centroid towards the most abundant species, and then computing a weighted average distance to this new centroid, using again the relative abundances as weights (\enc{Laliberté}{Laliberte} and Legendre submitted). FDis has no upper limit and requires at least two species to be computed. For communities composed of only one species, \code{dbFD} returns a FDis value of 0. FDis is by construction unaffected by species richness, it can be computed from any distance or dissimilarity measure (Anderson et al. 2006), it can handle any number and type of traits (including more traits than species), and it is not strongly influenced by outliers.

FDis is conceptually similar to Rao's quadratic entropy Q (\enc{Botta-Dukát}{Botta-Dukat} 2005), and simulations (via \code{\link{simul.dbFD}}) have shown high positive correlations between the two indices (\enc{Laliberté}{Laliberte} and Legendre submitted). Still, one potential advantage of FDis over Rao's Q is that in the unweighted case (i.e. with presence-absence data), it opens possibilities for formal statistical tests for differences in FD between two or more communities through a distance-based test for homogeneity of multivariate dispersions (Anderson 2006); see \code{\link[vegan]{betadisper}} for more details.

Corrections for PCoA axes corresponding to negative eigenvalues are applied following Anderson (2006); see also \code{\link[vegan]{betadisper}} for more details on these corrections.
}
\value{
  \item{FDis}{vector listing the FDis of each community}
  \item{eig}{vector listing the eigenvalues of the PCoA}
  \item{vectors}{matrix containing the PCoA axes}
}
\references{
Anderson, M. J. (2006) Distance-based tests for homogeneity of multivariate dispersions. \emph{Biometrics} \bold{62}:245-253.

Anderson, M. J., K. E. Ellingsen and B. H. McArdle (2006) Multivariate dispersion as a measure of beta diversity. \emph{Ecology Letters} \bold{9}:683-693.

\enc{Botta-Dukát}{Botta-Dukat}, Z. (2005) Rao's quadratic entropy as a measure of functional diversity based on multiple traits. \emph{Journal of Vegetation Science} \bold{16}:533-540.

\enc{Laliberté}{Laliberte}, E. and P. Legendre (submitted) A distance-based framework for measuring functional diversity from multiple traits. \emph{Ecology}.
}
\author{Etienne \enc{Laliberté}{Laliberte} \email{etiennelaliberte@gmail.com} \url{http://www.elaliberte.info}}
\note{
\code{fdisp} is implemented in \code{\link{dbFD}} and is used to compute the functional dispersion (FDis) index.
}
\seealso{\code{\link{dbFD}} for computing multidimensional FD indices and \code{\link[vegan]{betadisper}} from which \code{fdisp} borrows some code.}
\examples{
# dummy dataset
dummy.dist <- gowdis(dummy$trait)
ex1 <- fdisp(dummy.dist, dummy$abun)
ex1

# example with real data from New Zealand short-tussock grasslands
ex2 <- fdisp(gowdis(tussock$trait), tussock$abun)
ex2
}
\keyword{multivariate}
