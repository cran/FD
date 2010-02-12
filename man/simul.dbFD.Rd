\encoding{latin1}
\name{simul.dbFD}
\alias{simul.dbFD}
\concept{functional diversity}
\title{Simulations to Explore Relationships Between Functional Diversity Indices}
\description{
\code{simul.dbFD} generates artificial communities of species with artificial functional traits. Different functional diversity (FD) indices are computed from these communities using \code{\link{dbFD}} to explore their inter-relationships.
}
\usage{
simul.dbFD(s = c(5, 10, 15, 20, 25, 30, 35, 40), t = 3, 
          r = 10, p = 100, tr.method = c("unif", "norm", "lnorm"),
          abun.method = c("lnorm", "norm", "unif"), w.abun = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{vector listing the different levels of species richness used in the simulations}
  \item{t}{number of traits}
  \item{r}{number of replicates per species richness level}
  \item{p}{number of species in the common species pool}
  \item{tr.method}{character string indicating the sampling distribution for the traits. \code{"unif"} is a uniform distribution, \code{"norm"} is a normal distribution, and \code{"lnorm"} is a lognormal distribution.}
  \item{abun.method}{character string indicating the sampling distribution for the species abundances. Same as for \code{tr.method}.}
  \item{w.abun}{logical; should FDis, FEve, FDiv, and Rao's quadratic entropy (Q) be weighted by species abundances?}
}
\value{
A list contaning the following elements:

  \item{results}{data frame containing the results of the simulations}
  \item{traits}{matrix containing the traits}
  \item{abun}{matrix containing the abundances}
  \item{abun.gamma}{species abundances from the pooled set of communities}
  \item{FDis.gamma}{FDis of the pooled set of communities}
  \item{FDis.mean}{mean FDis from all communities}
\code{FDis.gamma} and \code{FDis.mean} can be used to explore the set concavity criterion (Ricotta 2005) for FDis.

A graph plotting the results of the simulations is also returned.
}
\references{
\enc{Laliberté}{Laliberte}, E. and P. Legendre (2010) A distance-based framework for measuring functional diversity from multiple traits. \emph{Ecology} \bold{91}299:305.

Ricotta, C. (2005) A note on functional diversity measures. \emph{Basic and Applied Ecology} \bold{6}:479-486.
}
\author{Etienne \enc{Laliberté}{Laliberte} \email{etiennelaliberte@gmail.com} \url{http://www.elaliberte.info}}
\section{Warning}{The simulations performed by \code{simul.dbFD} can take several hours if \code{length(s)} and/or \code{r} is large. Run a test with the default parameters first.} 
\seealso{\code{\link{dbFD}}, the function called in \code{simul.dbFD}}
\examples{
# this should take just a few minutes
\dontrun{
ex1 <- simul.dbFD(s = c(10, 20, 30, 40, 50), r = 5)
ex1
}
}
\keyword{datagen}
