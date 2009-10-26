\encoding{latin1}
\name{tussock}
\alias{tussock}
\concept{New Zealand}
\concept{functional composition}
\concept{functional trait}
\concept{functional diversity}
\docType{data}
\title{Functional Composition of Short-Tussock Grasslands}
\description{
\code{tussock} contains data on 16 functional traits measured on 53 vascular plant species from New Zealand short-tussock grasslands. It also contains the relative abundances (percent cover) of these 53 species from 30 8x50-m plots.
}
\usage{tussock}
\format{
\code{tussock} is a list of 2 components:
\describe{
\item{trait}{data frame of 16 functional traits measured on 53 plant species: growth form (sensu Cornelissen et al. 2003), reproductive plant height (\eqn{\mathrm{m}}{m}), leaf dry matter content (\eqn{\mathrm{mg} \; \mathrm{g}^{-1}}{mg g^-1}), leaf nitrogen concentration (\eqn{\mathrm{mg} \; \mathrm{g}^{-1}}{mg g^-1}), leaf phosphorous concentration (\eqn{\mathrm{mg} \; \mathrm{g}^{-1}}{mg g^-1}), leaf sulphur concentration (\eqn{\mathrm{mg} \; \mathrm{g}^{-1}}{mg g^-1}), specific leaf area (\eqn{\mathrm{m}^2 \; \mathrm{kg}^{-1}}{m^2 kg^-1}), nutrient uptake strategy (sensu Cornelissen et al. 2003), Raunkiaer life form, clonality, leaf size (\eqn{\mathrm{mm}^2}{mm^2}), primary dispersal mode, seed mass (\eqn{\mathrm{mg}}{mg}), resprouting capacity, pollination syndrome, and lifespan (an ordinal variable stored as \code{\link{ordered}}).}

\item{abun}{matrix containing the relative abundances (percent cover) of the 53 species in 30 plots}
}
}
\details{
The functional traits were measured using standardized methodologies (Cornelissen et al. 2003). Each of the 30 experimental plots from which species cover was estimated is 8x50 m. Relative abundances of all vascular plant species were estimated in November 2007. To do so, 20 1x1-m quadrats per plot were randomly positioned along two longitudinal transects and cover of each species was estimated using a modified Braun-Blanquet scale. This data was pooled at the plot scale to yield the percent cover data.
}
\source{
Etienne \enc{Laliberté}{Laliberte} \email{etiennelaliberte@gmail.com}

\url{http://www.elaliberte.info}
}
\references{
Cornelissen, J. H. C., S. Lavorel, E. Garnier, S. Diaz, N. Buchmann, D. E. Gurvich, P. B. Reich, H. ter Steege, H. D. Morgan, M. G. A. van der Heijden, J. G. Pausas and H. Poorter. (2003) A handbook of protocols for standardised and easy measurement of plant functional traits worldwide. \emph{Australian Journal of Botany} \bold{51}:335-380.


\enc{Laliberté}{Laliberte}, E., Norton, D. A. and D. Scott. (2008) Impacts of rangeland development on plant functional diversity, ecosystem processes and services, and resilience. \emph{Global Land Project (GLP) Newsletter} \bold{4}:4-6.

Scott, D. (1999) Sustainability of New Zealand high-country pastures under contrasting development inputs. 1. Site, and shoot nutrients. \emph{New Zealand Journal of Agricultural Research} \bold{42}:365-383.
}
\keyword{datasets}
