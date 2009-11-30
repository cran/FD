\encoding{latin1}
\name{FD-package}
\alias{FD-package}
\alias{FD}
\concept{functional diversity}
\concept{functional composition}
\concept{functional trait}
\concept{functional richness}
\concept{FRic}
\concept{functional evenness}
\concept{FEve}
\concept{functional divergence}
\concept{FDiv}
\concept{functional dispersion}
\concept{FDis}
\concept{quadratic entropy}
\concept{functional group richness}
\concept{FGR}
\concept{community weighted means}
\concept{CWM}
\docType{package}
\title{
Measuring Functional Diversity from Multiple Traits, and Other Tools for Functional Ecology
}
\description{
\pkg{FD} is a package to compute different multidimensional functional diversity (FD) indices. It implements a distance-based framework to measure FD that allows any number and type of functional traits, and can also consider species relative abundances. It also contains other tools for functional ecologists (e.g. \code{\link{maxent}}).
}
\details{
\tabular{ll}{
Package: \tab FD\cr
Type: \tab Package\cr
Version: \tab 1.0-5\cr
Date: \tab 2009-11-30\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes \cr
}


\pkg{FD} computes different multidimensional FD indices. To compute FD indices, a species-by-trait(s) matrix is required (or at least a species-by-species distance matrix). \code{\link{gowdis}} computes the Gower dissimilarity from different trait types (continuous, ordinal, nominal, or binary), and tolerates \code{\link{NA}}s. It can treat ordinal variables as described by Podani (1999), and can handle asymetric binary variables and variable weights. \code{\link{gowdis}} is called by \code{\link{dbFD}}, the main function of \pkg{FD}.  

\code{\link{dbFD}} uses principal coordinates analysis (PCoA) to return PCoA axes, which are then used as \sQuote{traits} to compute FD. \code{\link{dbFD}} computes several multidimensional FD indices, including the three indices of \enc{Villéger}{Villeger} et al. (2008): functional richness (FRic), functional evenness (FEve), and functional divergence (FDiv). It also computes functional dispersion (FDis) (\enc{Laliberté}{Laliberte} and Legendre 2009), Rao's quadratic entropy (Q) (\enc{Botta-Dukát}{Botta-Dukat} 2005), a posteriori functional group richness (FGR), and the community-level weighted means of trait values (CWM), an index of functional composition. Some of these indices can be weighted by species abundances. \code{\link{dbFD}} includes several options for flexibility.  
}
\author{
Etienne \enc{Laliberté}{Laliberte} and Bill Shipley

Maintainer: Etienne \enc{Laliberté}{Laliberte} <etiennelaliberte@gmail.com> \url{http://www.elaliberte.info}
}
\references{
\enc{Botta-Dukát}{Botta-Dukat}, Z. (2005) Rao's quadratic entropy as a measure of functional diversity based on multiple traits. \emph{Journal of Vegetation Science} \bold{16}:533-540.

\enc{Laliberté}{Laliberte}, E. and P. Legendre (2009) A distance-based framework for measuring functional diversity from multiple traits. \emph{Ecology} \bold{in press}.

Podani, J. (1999) Extending Gower's general coefficient of similarity to ordinal characters. \emph{Taxon} \bold{48}:331-340.

\enc{Villéger}{Villeger}, S., N. W. H. Mason and D. Mouillot (2008) New multidimensional functional diversity indices for a multifaceted framework in functional ecology. \emph{Ecology} \bold{89}:2290-2301.


}
\keyword{package}
\examples{
# examples with a dummy dataset

ex1 <- gowdis(dummy$trait)
ex1

ex2 <- functcomp(dummy$trait, dummy$abun)
ex2

ex3 <- dbFD(dummy$trait, dummy$abun)
ex3

# examples with real data from New Zealand short-tussock grasslands
# these examples may take a few seconds to a few minutes each to run

ex4 <- gowdis(tussock$trait)

ex5 <- functcomp(tussock$trait, tussock$abun)

# 'lingoes' correction used because 'sqrt' does not work in that case
ex6 <- dbFD(tussock$trait, tussock$abun, corr = "lingoes")

\dontrun{
# ward clustering to compute FGR, cailliez correction
ex7 <- dbFD(tussock$trait, tussock$abun, corr = "cailliez",
calc.FGR = TRUE, clust.type = "ward")
# choose 'g' for number of groups
# 6 groups seems to make good ecological sense
ex7

# however, calinksi criterion in 'kmeans' suggests
# that 6 groups may not be optimal
ex8 <- dbFD(tussock$trait, tussock$abun, corr = "cailliez",
calc.FGR = TRUE, clust.type = "kmeans", km.sup.gr = 10)
}
}
