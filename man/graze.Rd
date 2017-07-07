\name{graze}
\alias{graze}
\docType{data}
\title{ Site information and grazed vegetation data. }
\description{
  This data frame contains site location, landscape context and dominant plant species abundances for 25 of the plant species found in 50 grazed pastures in the northeastern United States. Elements are the mean values for canopy cover for ten 0.5 x 2 m quadrats.
}
\usage{data(graze)}
\format{
  A data frame with 50 observations on the following 25 variables.
  \describe{
    \item{\code{sitelocation}}{Site location along a geographic gradient.}
    \item{\code{forestpct}}{Percentage forest cover within a 1-km radius.}
    \item{\code{ACMI2}}{Percentage canopy cover.}
    \item{\code{ANOD}}{Percentage canopy cover.}
    \item{\code{ASSY}}{Percentage canopy cover.}
    \item{\code{BRIN2}}{Percentage canopy cover.}
    \item{\code{CIAR4}}{Percentage canopy cover.}
    \item{\code{CIIN}}{Percentage canopy cover.}
    \item{\code{CIVU}}{Percentage canopy cover.}
    \item{\code{DAGL}}{Percentage canopy cover.}
    \item{\code{ELRE4}}{Percentage canopy cover.}
    \item{\code{GAMO}}{Percentage canopy cover.}
    \item{\code{LOAR10}}{Percentage canopy cover.}
    \item{\code{LOCO6}}{Percentage canopy cover.}
    \item{\code{LOPE}}{Percentage canopy cover.}
    \item{\code{OXST}}{Percentage canopy cover.}
    \item{\code{PLMA2}}{Percentage canopy cover.}
    \item{\code{POPR}}{Percentage canopy cover.}
    \item{\code{PRVU}}{Percentage canopy cover.}
    \item{\code{RAAC3}}{Percentage canopy cover.}
    \item{\code{RUCR}}{Percentage canopy cover.}
    \item{\code{SORU2}}{Percentage canopy cover.}
    \item{\code{STGR}}{Percentage canopy cover.}
    \item{\code{TAOF}}{Percentage canopy cover.}
    \item{\code{TRPR2}}{Percentage canopy cover.}
    \item{\code{TRRE3}}{Percentage canopy cover.}
    \item{\code{VEOF2}}{Percentage canopy cover.}
  }
}
\details{
    Site locations fall along a southwest-northeast transect through the northeastern United States. This is a synthetic gradient calculated from latitude and longitude.
  Forest landcover is taken from the USGS 1992 National Land Cover Dataset. All forest classes were combined, and the percentage within 1 km of the sample sites was calculated using a GIS.
}
\source{
    Details of these data are available in Tracy and Sanderson (2000) and Goslee and Sanderson (2010).
    The 1992 NLCD data can be obtained from http://www.mrlc.gov/.
    Species codes are from http://plants.usda.gov (2010). 
}
\references{
  Tracy, B.F. and M.A. Sanderson. 2000. Patterns of plant species richness in pasture lands of the northeast United States. Plant Ecology 149:169-180.

  Goslee, S.C., Sanderson, M.A. 2010. Landscape Context and Plant Community Composition in Grazed Agricultural Systems. Landscape Ecology 25:1029-1039.
}
\examples{
data(graze)
}
\keyword{datasets}
