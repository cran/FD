\encoding{latin1}
\name{dbFD}
\alias{dbFD}
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
\title{Distance-Based Functional Diversity Indices}
\description{
\code{dbFD} implements a flexible distance-based framework to compute multidimensional functional diversity (FD) indices. \code{dbFD} returns the three FD indices of \enc{Villéger}{Villeger} et al. (2008): functional richness (FRic), functional evenness (FEve), and functional divergence (FDiv), as well functional dispersion (FDis; \enc{Laliberté}{Laliberte} and Legendre 2010), Rao's quadratic entropy (Q) (\enc{Botta-Dukát}{Botta-Dukat} 2005), a posteriori functional group richness (FGR) (Petchey and Gaston 2006), and the community-level weighted means of trait values (CWM; e.g. Lavorel et al. 2008). Some of these FD indices consider species abundances. \code{dbFD} includes several options for flexibility.
}
\usage{
dbFD(x, a, w, w.abun = TRUE, stand.x = TRUE,
    ord = c("podani", "metric"), asym.bin = NULL,
    corr = c("sqrt", "cailliez", "lingoes", "none"),
    calc.FRic = TRUE, m = "max", stand.FRic = FALSE,
    scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward",
    km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100,
    km.crit = c("calinski", "ssi"), calc.CWM = TRUE,
    CWM.type = c("dom", "all"), calc.FDiv = TRUE, dist.bin = 2, 
    print.pco = FALSE, messages = TRUE)
}
\arguments{
  \item{x}{matrix or data frame of functional traits. Traits can be \code{\link{numeric}}, \code{\link{ordered}}, or \code{\link{factor}}. Binary traits should be \code{\link{numeric}} and only contain 0 and 1. \code{\link{character}} traits will be converted to \code{\link{factor}}. \code{\link{NA}}s are tolerated.

\code{x} can also be a species-by-species distance matrix of class \code{\link{dist}}, in which case \code{\link{NA}s} are not allowed.

When there is only one trait, \code{x} can be also be a \code{\link{numeric}} vector, an \code{\link{ordered}} factor, or a unordered \code{\link{factor}}.

In all cases, species labels are required.
}.
  \item{a}{matrix containing the abundances of the species in \code{x} (or presence-absence, i.e. 0 or 1). Rows are sites and species are columns. Can be missing, in which case \code{dbFD} assumes that there is only one community with equal abundances of all species. \code{\link{NA}s} will be replaced by 0.
The number of species (columns) in \code{a} must match the number of species (rows) in \code{x}. In addition, the species labels in \code{a} and \code{x} must be identical and in the same order. 
}
  \item{w}{vector listing the weights for the traits in \code{x}. Can be missing, in which case all traits have equal weights.}
  \item{w.abun}{logical; should FDis, Rao's Q, FEve, FDiv, and CWM be weighted by the relative abundances of the species?}
  \item{stand.x}{logical; if all traits are \code{\link{numeric}}, should they be standardized to mean 0 and unit variance? If not all traits are \code{numeric}, Gower's (1971) standardization by the range is automatically used; see \code{\link{gowdis}} for more details.}
  \item{ord}{character string specifying the method to be used for ordinal traits (i.e. \code{\link{ordered}}). \code{"podani"} refers to Eqs. 2a-b of Podani (1999), while \code{"metric"} refers to his Eq. 3. Can be abbreviated. See \code{\link{gowdis}} for more details.}
  \item{asym.bin}{vector listing the asymmetric binary variables in \code{x}. See \code{\link{gowdis}} for more details.}
  \item{corr}{character string specifying the correction method to use when the species-by-species distance matrix cannot be represented in a Euclidean space. Options are \code{"sqrt"}, \code{"cailliez"}, \code{"lingoes"}, or \code{"none"}. Can be abbreviated. Default is \code{"sqrt"}. See \sQuote{details} section.}
  \item{calc.FRic}{logical; should FRic be computed?}
  \item{m}{the number of PCoA axes to keep as \sQuote{traits} for calculating FRic (when FRic is measured as the convex hull volume) and FDiv. Options are: any integer \eqn{>1}, \code{"min"} (maximum number of traits that allows the \eqn{s \geq 2^t}{s >= 2^t} condition to be met, where \eqn{s} is the number of species and \eqn{t} the number of traits), or \code{"max"} (maximum number of axes that allows the \eqn{s > t} condition to be met). See \sQuote{details} section.}
  \item{stand.FRic}{logical; should FRic be standardized by the \sQuote{global} FRic that include all species, so that FRic is constrained between 0 and 1?}
  \item{scale.RaoQ}{logical; should Rao's Q be scaled by its maximal value over all frequency distributions? See \code{\link[ade4]{divc}}.}
  \item{calc.FGR}{logical; should FGR be computed?}
  \item{clust.type}{character string specifying the clustering method to be used to create the dendrogram of species for FGR. Options are \code{"ward"}, \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"median"}, \code{"centroid"}, and \code{"kmeans"}. For \code{"kmeans"}, other arguments also apply (\code{km.inf.fr}, \code{km.sup.gr}, \code{km.iter}, and \code{km.crit}). See \code{\link{hclust}} and \code{\link[vegan]{cascadeKM}} for more details.}
  \item{km.inf.gr}{the number of groups for the partition with the smallest number of groups of the cascade (min). Only applies if \code{calc.FGR} is \code{TRUE} and \code{clust.type} is \code{"kmeans"}. See \code{\link[vegan]{cascadeKM}} for more details.}
  \item{km.sup.gr}{the number of groups for the partition with the largest number of groups of the cascade (max). Only applies if \code{calc.FGR} is \code{TRUE} and \code{clust.type} is \code{"kmeans"}. See \code{\link[vegan]{cascadeKM}} for more details.}
  \item{km.iter}{the number of random starting configurations for each value of \eqn{K}. Only applies if \code{calc.FGR} is \code{TRUE} and \code{clust.type} is \code{"kmeans"}. See \code{\link[vegan]{cascadeKM}} for more details.}
  \item{km.crit}{criterion used to select the best partition. The default value is \code{"calinski"} (Calinski-Harabasz 1974). The simple structure index \code{"ssi"} is also available. Only applies if \code{calc.FGR} is \code{TRUE} and \code{clust.type} is \code{"kmeans"}. Can be abbreviated. See \code{\link[vegan]{cascadeKM}} for more details.}
  \item{calc.CWM}{logical; should the community-level weighted means of trait values (CWM) be calculated? Can be abbreviated. See \code{\link{functcomp}} for more details.}
  \item{CWM.type}{character string indicating how nominal, binary and ordinal traits should be handled for CWM. See \code{\link{functcomp}} for more details.}
  \item{calc.FDiv}{logical; should FDiv be computed?}
  \item{dist.bin}{only applies when \code{x} is a single unordered \code{\link{factor}}, in which case \code{x} is coded using dummy variables. \code{dist.bin} is an integer between 1 and 10 specifying the appropriate distance measure for binary data. 2 (the default) refers to the simple matching coefficient (Sokal and Michener 1958). See \code{\link[ade4]{dist.binary}} for the other options.}
  \item{print.pco}{logical; should the eigenvalues and PCoA axes be returned?}
  \item{messages}{logical; should warning messages be printed in the console?}
}
\details{

Typical usage is
\preformatted{
dbFD(x, a, \dots)
}

If \code{x} is a matrix or a data frame that contains only continuous traits, no \code{\link{NA}s}, and that no weights are specified (i.e. \code{w} is missing), a species-species Euclidean distance matrix is computed via \code{\link{dist}}. Otherwise, a Gower dissimilarity matrix is computed via \code{\link{gowdis}}. If \code{x} is a distance matrix, it is taken as is.

When \code{x} is a single trait, species with \code{\link{NA}s} are first excluded to avoid \code{\link{NA}s} in the distance matrix. If \code{x} is a single continuous trait (i.e. of class \code{\link{numeric}}), a species-species Euclidean distance matrix is computed via \code{\link{dist}}. If \code{x} is a single ordinal trait (i.e. of class \code{\link{ordered}}), \code{\link{gowdis}} is used and argument \code{ord} applies. If \code{x} is a single nominal trait (i.e. an unordered \code{\link{factor}}), the trait is converted to dummy variables and a distance matrix is computed via \code{\link[ade4]{dist.binary}}, following argument \code{dist.bin}.

Once the species-species distance matrix is obtained, \code{dbFD} checks whether it is Euclidean. This is done via \code{\link[ade4]{is.euclid}}. PCoA axes corresponding to negative eigenvalues are imaginary axes that cannot be represented in a Euclidean space, but simply ignoring these axes would lead to biased estimations of FD. Hence in \code{dbFD} one of four correction methods are used, following argument \code{corr}. \code{"sqrt"} simply takes the square root of the distances. However, this approach does not always work for all coefficients, in which case \code{dbFD} will stop and tell the user to select another correction method. \code{"cailliez"} refers to the approach described by Cailliez (1983) and is implemented via \code{\link[ade4]{cailliez}}. \code{"lingoes"} refers to the approach described by Lingoes (1971) and is implemented via \code{\link[ade4]{lingoes}}. \code{"none"} creates a distance matrix with only the positive eigenvalues of the Euclidean representation via \code{\link[ade4]{quasieuclid}}. See Legendre and Legendre (1998) and Legendre and Anderson (1999) for more details on these corrections.

Principal coordinates analysis (PCoA) is then performed (via \code{\link[ade4]{dudi.pco}}) on the \emph{corrected} species-species distance matrix. The resulting PCoA axes are used as the new \sQuote{traits} to compute the three indices of \enc{Villéger}{Villeger} et al. (2008): FRic, FEve, and FDiv. For FEve, there is no limit on the number of traits that can be used, so all PCoA axes are used. On the other hand, FRic and FDiv both rely on finding the minimum convex hull that includes all species (\enc{Villéger}{Villeger} et al. 2008). This requires more species than traits. To circumvent this problem, \code{dbFD} takes only a subset of the PCoA axes as traits via argument \code{m}. This, however, comes at a cost of loss of information. The quality of the resulting reduced-space representation is returned by \code{qual.FRic}, which is computed as described by Legendre and Legendre (1998) and can be interpreted as a \eqn{R^2}{R^2}-like ratio.

In \code{dbFD}, FRic is generally measured as the convex hull volume, but when there is only one continuous trait it is measured as the range (or the range of the ranks for an ordinal trait). Conversely, when only nominal and ordinal traits are present, FRic is measured as the number of unique trait value combinations in a community. FEve and FDiv, but not FRic, can account for species relative abundances, as described by \enc{Villéger}{Villeger} et al. (2008).

Functional dispersion (FDis; \enc{Laliberté}{Laliberte} and Legendre 2010) is computed from the \emph{uncorrected} species-species distance matrix via \code{\link{fdisp}}. Axes with negatives eigenvalues are corrected following the approach of Anderson (2006). When all species have equal abundances (i.e. presence-absence data), FDis is simply the average distance to the centroid (i.e. multivariate dispersion) as originally described by Anderson (2006). Multivariate dispersion has been proposed as an index of beta diversity (Anderson et al. 2006). However, \enc{Laliberté}{Laliberte} and Legendre (2010) have extended it to a FD index. FDis can account for relative abundances by shifting the position of the centroid towards the most abundant species, and then computing a weighted average distance to this new centroid, using again the relative abundances as weights (\enc{Laliberté}{Laliberte} and Legendre 2010). FDis has no upper limit and requires at least two species to be computed. For communities composed of only one species, \code{dbFD} returns a FDis value of 0. FDis is by construction unaffected by species richness, it can be computed from any distance or dissimilarity measure (Anderson et al. 2006), it can handle any number and type of traits (including more traits than species), and it is not strongly influenced by outliers.

Rao's quadratic entropy (Q) is computed from the \emph{uncorrected} species-species distance matrix via \code{\link[ade4]{divc}}. See \enc{Botta-Dukát}{Botta-Dukat} (2005) for details. Rao's Q is conceptually similar to FDis, and simulations (via \code{\link{simul.dbFD}}) have shown high positive correlations between the two indices (\enc{Laliberté}{Laliberte} and Legendre 2010). Still, one potential advantage of FDis over Rao's Q is that in the unweighted case (i.e. with presence-absence data), it opens possibilities for formal statistical tests for differences in FD between two or more communities through a distance-based test for homogeneity of multivariate dispersions (Anderson 2006); see  \code{\link[vegan]{betadisper}} for more details.

Functional group richness (FGR) is based on the classification of the species by the user from visual inspection of a dengrogram. Method \code{"kmeans"} is also available by calling \code{\link[vegan]{cascadeKM}}. In that case, the Calinski-Harabasz (1974) criterion or the simple structure index (SSI) can be used to estimate the number of functional groups; see \code{\link[vegan]{cascadeKM}} for more details. FGR returns the number of functional groups per community, as well as the abundance of each group in each community. 

The community-level means of trait values (CWM) is an index of functional composition (Lavorel et al. 2008), and is computed via \code{\link{functcomp}}. Species with \code{\link{NA}s} for a given trait are excluded for that trait.
}
\value{
  \item{nbsp}{vector listing the number of species in each community}
  \item{sing.sp}{vector listing the number of functionally singular species in each community. If all species are functionally different, \code{sing.sp} will be identical to \code{nbsp}.}
  \item{FRic}{vector listing the FRic of each community}
  \item{qual.FRic}{quality of the reduced-space representation required to compute FRic and FDiv.}
  \item{FEve}{vector listing the FEve of each community}
  \item{FDiv}{vector listing the FDiv of each community. Only returned if \code{calc.FDiv} is \code{TRUE}.}
  \item{FDis}{vector listing the FDis of each community}
  \item{RaoQ}{vector listing the Rao's quadratic entropy (Q) of each community}
  \item{FGR}{vector listing the FGR of each community. Only returned if \code{calc.FGR} is \code{TRUE}.}
  \item{spfgr}{vector specifying functional group membership for each species. Only returned if \code{calc.FGR} is \code{TRUE}.}
  \item{gr.abun}{matrix containing the abundances of each functional group in each community. Only returned if \code{calc.FGR} is \code{TRUE}.}
  \item{CWM}{data frame containing the community-level weighted trait means (CWM). Only returned if \code{calc.CWM} is \code{TRUE}.}
  \item{x.values}{eigenvalues from the PCoA. Only returned if \code{print.pco} is \code{TRUE}.}
  \item{x.axes}{PCoA axes. Only returned if \code{print.pco} is \code{TRUE}.}
}
\references{
Anderson, M. J. (2006) Distance-based tests for homogeneity of multivariate dispersions. \emph{Biometrics} \bold{62}:245-253.

Anderson, M. J., K. E. Ellingsen and B. H. McArdle (2006) Multivariate dispersion as a measure of beta diversity. \emph{Ecology Letters} \bold{9}:683-693.

\enc{Botta-Dukát}{Botta-Dukat}, Z. (2005) Rao's quadratic entropy as a measure of functional diversity based on multiple traits. \emph{Journal of Vegetation Science} \bold{16}:533-540.

Cailliez, F. (1983) The analytical solution of the additive constant problem. \emph{Psychometrika} \bold{48}:305-310.

Calinski, T. and J. Harabasz (1974) A dendrite method for cluster analysis. \emph{Communications in Statistics} \bold{3}:1-27.

Gower, J. C. (1971) A general coefficient of similarity and some of its properties. \emph{Biometrics} \bold{27}:857-871.

\enc{Laliberté}{Laliberte}, E. and P. Legendre (2010) A distance-based framework for measuring functional diversity from multiple traits. \emph{Ecology} \bold{91}:299-305.

Lavorel, S., K. Grigulis, S. McIntyre, N. S. G. Williams, D. Garden, J. Dorrough, S. Berman, F. \enc{Quétier}{Quetier}, A. Thebault and A. Bonis (2008) Assessing functional diversity in the field - methodology matters! \emph{Functional Ecology} \bold{22}:134-147.

Legendre, P. and M. J. Anderson (1999) Distance-based redundancy analysis: testing multispecies responses in multifactorial ecological experiments. \emph{Ecological Monographs} \bold{69}:1-24.

Legendre, P. and L. Legendre (1998) \emph{Numerical Ecology}. 2nd English edition. Amsterdam: Elsevier.

Lingoes, J. C. (1971) Some boundary conditions for a monotone analysis of symmetric matrices. \emph{Psychometrika} \bold{36}:195-203.

Podani, J. (1999) Extending Gower's general coefficient of similarity to ordinal characters. \emph{Taxon} \bold{48}:331-340.

Sokal, R. R. and C. D. Michener (1958) A statistical method for evaluating systematic relationships. \emph{The University of Kansas Scientific Bulletin} \bold{38}:1409-1438.

\enc{Villéger}{Villeger}, S., N. W. H. Mason and D. Mouillot (2008) New multidimensional functional diversity indices for a multifaceted framework in functional ecology. \emph{Ecology} \bold{89}:2290-2301.
}
\author{Etienne \enc{Laliberté}{Laliberte} \email{etiennelaliberte@gmail.com} \url{https://www.elaliberte.info/}}
\note{\code{dbFD} borrows code from the \code{F_RED} function of \enc{Villéger}{Villeger} et al. (2008).}
\section{Warning}{Users often report that \code{dbFD} crashed during their analysis. Generally this occurs under Windows, and is almost always due to the computation of convex hull volumes. Possible solutions are to choose \code{calc.FRic = "FALSE"}, or to reduce the dimensionality of the trait matrix using the \code{"m"} argument.} 

\seealso{\code{\link{gowdis}}, \code{\link{functcomp}}, \code{\link{fdisp}}, \code{\link{simul.dbFD}}, \code{\link[ade4]{divc}}, \code{\link[vegan]{treedive}}, \code{\link[vegan]{betadisper}}}
\examples{
# mixed trait types, NA's
ex1 <- dbFD(dummy$trait, dummy$abun)
ex1

# add variable weights
# 'cailliez' correction is used because 'sqrt' does not work
w<-c(1, 5, 3, 2, 5, 2, 6, 1)
ex2 <- dbFD(dummy$trait, dummy$abun, w, corr="cailliez")

# if 'x' is a distance matrix
trait.d <- gowdis(dummy$trait)
ex3 <- dbFD(trait.d, dummy$abun)
ex3

# one numeric trait, one NA
num1 <- dummy$trait[,1] ; names(num1) <- rownames(dummy$trait)
ex4 <- dbFD(num1, dummy$abun)
ex4

# one ordered trait, one NA
ord1 <- dummy$trait[,5] ; names(ord1) <- rownames(dummy$trait)
ex5 <- dbFD(ord1, dummy$abun)
ex5

# one nominal trait, one NA
fac1 <- dummy$trait[,3] ; names(fac1) <- rownames(dummy$trait)
ex6 <- dbFD(fac1, dummy$abun)
ex6

# example with real data from New Zealand short-tussock grasslands
# 'lingoes' correction used because 'sqrt' does not work in that case
ex7 <- dbFD(tussock$trait, tussock$abun, corr = "lingoes")


\dontrun{
# calc.FGR = T, 'ward'
ex7 <- dbFD(dummy$trait, dummy$abun, calc.FGR = T)
ex7

# calc.FGR = T, 'kmeans'
ex8 <- dbFD(dummy$trait, dummy$abun, calc.FGR = T,
clust.type = "kmeans")
ex8

# ward clustering to compute FGR
ex9 <- dbFD(tussock$trait, tussock$abun,
corr = "cailliez", calc.FGR = TRUE, clust.type = "ward")
# choose 'g' for number of groups
# 6 groups seems to make good ecological sense
ex9

# however, calinksi criterion in 'kmeans' suggests
# that 6 groups may not be optimal
ex10 <- dbFD(tussock$trait, tussock$abun, corr = "cailliez",
calc.FGR = TRUE, clust.type = "kmeans", km.sup.gr = 10)
}
}
\keyword{multivariate}
