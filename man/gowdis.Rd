\encoding{latin1}
\name{gowdis}
\alias{gowdis}
\title{Gower Dissimilarity}
\description{
\code{gowdis} measures the Gower (1971) dissimilarity for mixed variables, including asymmetric binary variables. Variable weights can be specified. \code{gowdis} implements Podani's (1999) extension to ordinal variables.
}
\usage{
gowdis(x, w, asym.bin = NULL, ord = c("podani", "metric", "classic"))
}
\arguments{
  \item{x}{matrix or data frame containing the variables. Variables can be \code{\link{numeric}}, \code{\link{ordered}}, or \code{\link{factor}}. Symmetric or asymmetric binary variables should be \code{\link{numeric}} and only contain 0 and 1. \code{\link{character}} variables will be converted to \code{\link{factor}}. \code{\link{NA}}s are tolerated.}
  \item{w}{vector listing the weights for the variables in \code{x}. Can be missing, in which case all variables have equal weights.}
  \item{asym.bin}{vector listing the asymmetric binary variables in \code{x}.}
  \item{ord}{character string specifying the method to be used for ordinal variables (i.e. \code{\link{ordered}}). \code{"podani"} refers to Eqs. 2a-b of Podani (1999), while \code{"metric"} refers to his Eq. 3 (see \sQuote{details}); both options convert ordinal variables to ranks. \code{"classic"} simply treats ordinal variables as continuous variables. Can be abbreviated.}
}
\details{
\code{gowdis} computes the Gower (1971) similarity coefficient exactly as described by Podani (1999), then converts it to a dissimilarity coefficient by using \eqn{D = 1 - S}. It integrates variable weights as described by Legendre and Legendre (1998).

Let \eqn{\mathbf{X} = \{x_{ij}\} }{X = {Xij}} be a matrix containing \eqn{n} objects (rows) and \eqn{m} columns (variables). The similarity \eqn{G_{jk}}{Gjk} between objects \eqn{j} and \eqn{k} is computed as

    \deqn{G_{jk} = \frac{\sum_{i=1}^{n} w_{ijk} s_{ijk}}{\sum_{i=1}^{n} w_{ijk}}}{%
    Gjk = sum(Wijk * Sijk) / sum(Wijk)},

where \eqn{w_{ijk}}{Wijk} is the weight of variable \eqn{i} for the \eqn{j}-\eqn{k} pair, and \eqn{s_{ijk}}{Sijk} is the partial similarity of variable \eqn{i} for the \eqn{j}-\eqn{k} pair,

and where \eqn{w_{ijk} = 0}{Wijk = 0} if objects \eqn{j} and \eqn{k} cannot be compared because \eqn{x_{ij}}{Xij} or \eqn{x_{ik}}{Xik} is unknown (i.e. \code{\link{NA}}).


For binary variables, \eqn{s_{ijk} = 0}{Sijk = 0} if \eqn{x_{ij} \neq x_{ik}}{Xij is not equal to Xik}, and  \eqn{s_{ijk} = 1}{Sijk = 1} if \eqn{x_{ij} = x_{ik} = 1}{Xij = Xik = 1} or if \eqn{x_{ij} = x_{ik} = 0}{Xij = Xik = 0}.

For asymmetric binary variables, same as above except that \eqn{w_{ijk} = 0}{Wijk = 0} if \eqn{x_{ij} = x_{ik} = 0}{Xij = Xik = 0}.

For nominal variables, \eqn{s_{ijk} = 0}{Sijk = 0} if \eqn{x_{ij} \neq x_{ik}}{Xij is not equal to Xik} and \eqn{s_{ijk} = 1}{Sijk = 1} if \eqn{x_{ij} = x_{ik}}{Xij = Xik}.

For continuous variables,

    \deqn{s_{ijk} = 1 - \frac{|x_{ij} - x_{ik}|} {x_{i.max} - x_{i.min}} }{%
    Sijk = 1 - [ |Xij - Xik| / (Xi.max - Xi.min) ]}

where \eqn{x_{i.max}}{Xi.max} and \eqn{x_{i.min}}{Xi.min} are the maximum and minimum values of variable \eqn{i}, respectively.


For ordinal variables, when \code{ord = "podani"} or \code{ord = "metric"}, all \eqn{x_{ij}}{Xij} are replaced by their ranks \eqn{r_{ij}}{Rij} determined over all objects (such that ties are also considered), and then


if \code{ord = "podani"}

    \eqn{s_{ijk} = 1}{Sijk = 1} if \eqn{r_{ij} = r_{ik}}{Rij = Rik}, otherwise

    \deqn{ s_{ijk} = 1 - \frac{|r_{ij} - r_{ik}| - (T_{ij} - 1)/2 - (T_{ik} - 1)/2 }{r_{i.max} - r_{i.min} - (T_{i.max} - 1)/2 - (T_{i.min}-1)/2 }}{%
    Sijk = 1 - [ |Rij - Rik| - (Tij - 1) / 2 - (Tik - 1) / 2  / Ri.max - Ri.min - (Ti.max - 1) / 2 - (Ti.min - 1) / 2 ]}

where \eqn{T_{ij}}{Tij} is the number of objects which have the same rank score for variable \eqn{i} as object \eqn{j} (including \eqn{j} itself), \eqn{T_{ik}}{Tik} is the number of objects which have the same rank score for variable \eqn{i} as object \eqn{k} (including \eqn{k} itself), \eqn{r_{i.max}}{Ri.max} and \eqn{r_{i.min}}{Ri.min} are the maximum and minimum ranks for variable \eqn{i}, respectively, \eqn{T_{i,max}}{Ti.max} is the number of objects with the maximum rank, and \eqn{T_{i.min}}{Ti.min} is the number of objects with the minimum rank.


if \code{ord = "metric"}

    \deqn{s_{ijk} = 1 - \frac{|r_{ij} - r_{ik}|}{r_{i.max} - r_{i.min}} }{%
    Sijk = 1 - [ |Rij - Rik| / (Ri.max - Ri.min) ]}

When \code{ord = "classic"}, ordinal variables are simply treated as continuous variables.
}

\value{
an object of class \code{\link{dist}} with the following attributes: \code{Labels}, \code{Types} (the variable types, where 'C' is continuous/numeric, 'O' is ordinal, 'B' is symmetric binary, 'A' is asymmetric binary, and 'N' is nominal), \code{Size}, \code{Metric}.
}
\references{
Gower, J. C. (1971) A general coefficient of similarity and some of its properties. \emph{Biometrics} \bold{27}:857-871.

Legendre, P. and L. Legendre (1998) \emph{Numerical Ecology}. 2nd English edition. Amsterdam: Elsevier.

Podani, J. (1999) Extending Gower's general coefficient of similarity to ordinal characters. \emph{Taxon} \bold{48}:331-340.
}
\author{Etienne \enc{Laliberté}{Laliberte} \email{etiennelaliberte@gmail.com} \url{http://www.elaliberte.info}, with some help from Philippe Casgrain for the C interface.}
\seealso{\code{\link[cluster]{daisy}} is similar but less flexible, since it does not include variable weights and does not treat ordinal variables as described by Podani (1999). Using \code{ord = "classic"} reproduces the behaviour of \code{\link[cluster]{daisy}}.}
\examples{
ex1 <- gowdis(dummy$trait)
ex1

# check attributes
attributes(ex1)

# to include weights
w <- c(4,3,5,1,2,8,3,6)
ex2 <- gowdis(dummy$trait, w)
ex2

# variable 7 as asymmetric binary
ex3 <- gowdis(dummy$trait, asym.bin = 7)
ex3

# example with trait data from New Zealand vascular plant species
ex4 <- gowdis(tussock$trait)
}
\keyword{multivariate}
