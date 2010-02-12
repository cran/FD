\encoding{latin1}
\name{maxent.test}
\Rdversion{1.1}
\alias{maxent.test}
\concept{maximum entropy}
\concept{statistical mechanics}
\title{
Inferential Permutation Tests for Maximum Entropy Models
}
\description{
\code{maxent.test} performs the permutation tests proposed by Shipley (2010) for maximum entropy models.
Two different null hypotheses can be tested:
1) the information encoded in the \emph{entire set} of constraints \eqn{\mathbf{C}}{C} is irrelevant for predicting the probabilities, and
2) the information encoded in \emph{subset} \eqn{\mathbf{B}}{B} of the entire set of constraints \eqn{\mathbf{C}=\mathbf{A}\cup{}\mathbf{B}}{C = {A U B}} is irrelevant for predicting the probabilities.
A plot can be returned to facilitate interpretation.
}

\usage{
maxent.test(model, obs, sub.c, nperm = 99, quick = FALSE,
alpha = 0.05, plot = TRUE)
}

\arguments{
  \item{model}{
list returned by \code{\link{maxent}}.
}
  \item{obs}{
vector, matrix or data frame of observed probabilities of the states (columns).
}
  \item{sub.c}{
character or numeric vector specifying the subset of constraints \eqn{\mathbf{B}}{B} associated with null hypothesis 2.
If missing, null hypothesis 1 is tested.
}
  \item{nperm}{
number of permutations for the test.
}
  \item{quick}{
if \code{TRUE}, the algorithm stops when the upper confidence interval of the \emph{P}-value
is equal to or smaller than \code{alpha}. Can be useful to speed up the routine.
}

  \item{alpha}{
desired alpha-level for the test. Only relevant if \code{quick} is \code{TRUE}.
}

  \item{plot}{
if \code{TRUE}, a plot is returned to facilitate interpretation.
}

}
\details{
\code{maxent.test} is a direct translation of the permutation tests described by Shipley (2010). Please refer to this article for details.

If \code{\link{maxent}} was run with multiple data sets (i.e. if \code{constr} had more than one row), then  \code{maxent.test} performs the test for all sets simultaneously, following the \sQuote{omnibus} procedure described by Shipley (2010).

The following measure of fit between observed and predicted probabilities is returned:

\deqn{\mathrm{fit}=1-\frac{{\displaystyle \sum_{j=1}^{D} \sum_{i=1}^{S}\left(o_{ij}-p_{ij}\right)^{2}}}{\displaystyle \sum_{j=1}^{D} \sum_{i=1}^{S}\left(o_{ij}-q_{ij}\right)^{2}}}{%
fit = 1 - (sum(o_ij - p_ij)^2 / sum(o_ij - q_ij)^2)}

where \eqn{o_{ij}}{o_ij}, \eqn{p_{ij}}{o_ij}, and \eqn{q_{ij}}{q_ij} are the observed, predicted
and prior probabilities of state \eqn{i}{i} from data set \eqn{j}{j}, respectively, \eqn{S}{S} is
the number of states, and \eqn{D}{D} the number of data sets (i.e. rows in \code{obs}). A value of 1 indicates perfect predictive capacity, while a value near zero
indicates that the constraints provide no additional information beyond what is already contained in the
prior \eqn{\mathbf{q}}{q} (Sonnier et al. 2009).
}

\value{
  \item{fit}{measure of fit giving the predictive ability of the entire set of constraints \eqn{\mathbf{C}}{C},
beyond that already provided by the prior distribution. 
}
  \item{fit.a}{measure of fit giving the predictive ability of the subset of constraints \eqn{\mathbf{A}}{A},
beyond that already provided by the prior distribution;  only returned if \code{sub.c} is specified}
  \item{r2}{Pearson \eqn{r^2}{r^2} between observed and predicted probabilities, using
the entire set of constraints \eqn{\mathbf{C}}{C}}
  \item{r2.a}{Pearson \eqn{r^2}{r^2} between observed and predicted probabilities, using
the subset of constraints \eqn{\mathbf{A}}{A}; only returned if \code{sub.c} is specified}
  \item{r2.q}{Pearson \eqn{r^2}{r^2} between observed and prior probabilities;
only returned when \code{sub.c} is missing}
  \item{obs.stat}{observed statistic used for the permutation test; see Shipley (2010)}
  \item{nperm}{number of permutations; can be smaller than the specified \code{nperm} when \code{quick} is \code{TRUE}}
  \item{pval}{\emph{P}-value}
  \item{ci.pval}{approximate confidence intervals of the \emph{P}-value}
}
\references{

Sonnier, G., Shipley, B., and M. L. Navas. 2009. Plant traits, species pools and the prediction of relative abundance in plant communities: a maximum entropy approach. \emph{Journal of Vegetation Science} \bold{in press}.

Shipley, B. (2010) Inferential permutation tests for maximum entropy models in ecology. \emph{Ecology} \bold{in press}.
}

\author{
Etienne \enc{Laliberté}{Laliberte} \email{etiennelaliberte@gmail.com}

\url{http://www.elaliberte.info}
}

\section{Warning}{\code{maxent.test} is a computationally intensive function. The tests can take a long time when \code{nperm} is large.} 

\seealso{
\code{\link{maxent}} to run the maximum entropy model that is required by \code{maxent.test}. 
}
\examples{

# example with tussock dataset
traits <- tussock$trait[, c(2:7, 11)] # use only continuous traits
traits <- na.omit(traits) # remove 2 species with NA's
abun <- tussock$abun[, rownames(traits)] # abundance matrix
abun <- t(apply(abun, 1, function(x) x / sum(x) )) # relative abundances
agg <- functcomp(traits, abun) # community-aggregated traits
traits <- t(traits) # transpose matrix

# run maxent on site 1 (first row of abun), all species
pred.abun <- maxent(agg[1,], traits)

\dontrun{
# do the constraints give predictive ability?
maxent.test(pred.abun, obs = abun[1,], nperm = 49)

# are height, LDMC, and leaf [N] important constraints?
maxent.test(pred.abun, obs = abun[1,], sub.c = c("height",
"LMDC", "leafN"), nperm = 49)
}
}

\keyword{distribution}
\keyword{math}
\keyword{models}
