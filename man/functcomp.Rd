\encoding{latin1}
\name{functcomp}
\alias{functcomp}
\concept{functional diversity}
\concept{functional composition}
\concept{functional trait}
\concept{functional identity}
\concept{community weighted means}
\concept{CWM}
\title{Functional Composition}
\description{
\code{functcomp} returns the functional composition of a set of communities, as measured by the community-level weighted means of trait values (CWM; e.g. Lavorel et al. 2008).  
}
\usage{
functcomp(x, a, CWM.type = c("dom", "all"), bin.num = NULL)
}
\arguments{
  \item{x}{matrix or data frame containing the functional traits. Traits can be \code{\link{numeric}}, \code{\link{ordered}}, or \code{\link{factor}}. Binary traits should be \code{\link{numeric}} and only contain 0 and 1. \code{\link{character}} traits will be converted to \code{\link{factor}}. For a given trait, species with \code{\link{NA}} are excluded.}
  \item{a}{matrix containing the abundances of the species in \code{x} (or presence-absence, i.e. 0 or 1). Rows are sites and species are columns. The number of species (columns) in \code{a} must match the number of species (rows) in \code{x}. In addition, the species labels in \code{a} and \code{x} must be identical and in the same order. \code{\link{NA}s} will be replaced by 0.}
  \item{CWM.type}{character string indicating how nominal, binary and ordinal traits should be handled. See \sQuote{details}.}
  \item{bin.num}{vector indicating binary traits to be treated as continuous.}
}
\details{
\code{functcomp} computes the community-level weighted means of trait values for a set of communities (i.e. sites). For a continuous trait, CWM is the mean trait value of all species present in the community (after excluding species with \code{\link{NA}s}), weighted by their relative abundances.

For ordinal, nominal and binary traits, either the dominant class is returned (when \code{CWM.type} is \code{"dom"}), or the abundance of each individual class is returned (when \code{CWM.type} is \code{"all"}).

The default behaviour of binary traits being treated as nominal traits can be over-ridden by specifying \code{bin.num}, in which case they are treated as numeric traits.

When \code{CWM.type = "dom"}, if the maximum abundance value is shared between two or more classes, then one of these classes is randomly selected for CWM. Because species with \code{\link{NA}s} for a given trait are excluded for that trait, it is possible that when \code{CWM.type} is set to \code{"all"}, the sum of the abundances of all classes for a given ordinal/nominal/binary trait does not equal the sum of the species abundances. Thus, it is definitely not recommended to have \code{\link{NA}s} for very abundant species, as this will lead to biased estimates of functional composition.  
}
\value{
a data frame containing the CWM values of each trait for each community.
}
\references{
Lavorel, S., K. Grigulis, S. McIntyre, N. S. G. Williams, D. Garden, J. Dorrough, S. Berman, F. \enc{Quétier}{Quetier}, A. Thébault and A. Bonis (2008) Assessing functional diversity in the field - methodology matters! \emph{Functional Ecology} \bold{22}:134-147.
}
\author{Etienne \enc{Laliberté}{Laliberte} \email{etiennelaliberte@gmail.com} \url{https://www.elaliberte.info/}}
\note{
\code{functcomp} is implemented in \code{\link{dbFD}} and will be returned if \code{calc.CWM} is \code{TRUE}. 
}
\seealso{\code{\link{dbFD}} for measuring distance-based multidimensional functional diversity indices, including CWM.}
\examples{
# for ordinal, nominal and binary variables
# returns only the most frequent class
ex1 <- functcomp(dummy$trait, dummy$abun)
ex1

# returns the frequencies of each class
ex2 <- functcomp(dummy$trait, dummy$abun, CWM.type = "all")
ex2

# example with real data from New Zealand short-tussock grasslands
ex3 <- functcomp(tussock$trait, tussock$abun)
ex3
}
\keyword{multivariate}
