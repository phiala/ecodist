\name{graze}
\alias{graze}
\docType{data}
\title{ Site information and grazed vegetation data. }
\description{
  These data contain site location, landscape context and dominant plant species abundances for 50 grazed pastures in the northeastern United States. Elements are the mean values for canopy cover for ten 0.5 x 2 m quadrats.
}
\usage{data(graze)}
\format{
  A data frame with 50 observations on the following 67 variables.
  \describe{
    \item{\code{sitelocation}}{ Site location along a geographic gradient. }
    \item{\code{forestpct}}{Percentage forest cover within a 1-km radius.}
    \item{\code{ACMI2}}{percentage canopy cover}
    \item{\code{AGGI2}}{percentage canopy cover}
    \item{\code{AGROS2}}{percentage canopy cover}
    \item{\code{ANAR6}}{percentage canopy cover}
    \item{\code{ANOD}}{percentage canopy cover}
    \item{\code{ARMI2}}{percentage canopy cover}
    \item{\code{ASSY}}{percentage canopy cover}
    \item{\code{BRCA5}}{percentage canopy cover}
    \item{\code{BRIN2}}{percentage canopy cover}
    \item{\code{CEFO2}}{percentage canopy cover}
    \item{\code{CIAR4}}{percentage canopy cover}
    \item{\code{CIIN}}{percentage canopy cover}
    \item{\code{CIVU}}{percentage canopy cover}
    \item{\code{DACA6}}{percentage canopy cover}
    \item{\code{DAGL}}{percentage canopy cover}
    \item{\code{DRYOP}}{percentage canopy cover}
    \item{\code{ELRE4}}{percentage canopy cover}
    \item{\code{ERAN}}{percentage canopy cover}
    \item{\code{EUGRG}}{percentage canopy cover}
    \item{\code{FERU2}}{percentage canopy cover}
    \item{\code{FRVI}}{percentage canopy cover}
    \item{\code{GAMO}}{percentage canopy cover}
    \item{\code{HIERA}}{percentage canopy cover}
    \item{\code{HOLA}}{percentage canopy cover}
    \item{\code{HYPU}}{percentage canopy cover}
    \item{\code{HYRA3}}{percentage canopy cover}
    \item{\code{JUEF}}{percentage canopy cover}
    \item{\code{LEAU2}}{percentage canopy cover}
    \item{\code{LEVU}}{percentage canopy cover}
    \item{\code{LOAR10}}{percentage canopy cover}
    \item{\code{LOCO6}}{percentage canopy cover}
    \item{\code{LOPE}}{percentage canopy cover}
    \item{\code{MOVE}}{percentage canopy cover}
    \item{\code{OXST}}{percentage canopy cover}
    \item{\code{PEGL2}}{percentage canopy cover}
    \item{\code{PHAR3}}{percentage canopy cover}
    \item{\code{PHPR3}}{percentage canopy cover}
    \item{\code{PLLA}}{percentage canopy cover}
    \item{\code{PLMA2}}{percentage canopy cover}
    \item{\code{POAV}}{percentage canopy cover}
    \item{\code{POCA17}}{percentage canopy cover}
    \item{\code{POLA4}}{percentage canopy cover}
    \item{\code{POLYG}}{percentage canopy cover}
    \item{\code{PONO3}}{percentage canopy cover}
    \item{\code{POPR}}{percentage canopy cover}
    \item{\code{POSI2}}{percentage canopy cover}
    \item{\code{POTEN}}{percentage canopy cover}
    \item{\code{PRVU}}{percentage canopy cover}
    \item{\code{RAAC3}}{percentage canopy cover}
    \item{\code{RUBUS}}{percentage canopy cover}
    \item{\code{RUCR}}{percentage canopy cover}
    \item{\code{SIOF}}{percentage canopy cover}
    \item{\code{SOLID}}{percentage canopy cover}
    \item{\code{SORU2}}{percentage canopy cover}
    \item{\code{STGR}}{percentage canopy cover}
    \item{\code{STME2}}{percentage canopy cover}
    \item{\code{SYLA4}}{percentage canopy cover}
    \item{\code{TAOF}}{percentage canopy cover}
    \item{\code{THVU}}{percentage canopy cover}
    \item{\code{TRAU2}}{percentage canopy cover}
    \item{\code{TRHY}}{percentage canopy cover}
    \item{\code{TRPR2}}{percentage canopy cover}
    \item{\code{TRRE3}}{percentage canopy cover}
    \item{\code{VEOF2}}{percentage canopy cover}
    \item{\code{VICIA}}{percentage canopy cover}
  }
}
\details{
    Site locations fall along a southwest-northeast transect through the northeastern United States. This is a synthetic gradient calculated from latitude and longitude.
  Forest landcover is taken from the USGS 1992 National Land Cover Dataset. All forest classes were combined, and the percentage within 1 km of the sample sites was calculated using a GIS.
  Species codes follow the USDA Plants database. Species with a maximum cover <= 1 have been omitted.
}
\source{
    Details of these data are available in Tracy and Sanderson (2000) and Goslee, Sanderson and Tracy (ms in review).
    The 1992 NLCD data can be obtained from http://www.mrlc.gov/.
    Species codes are from http://plants.usda.gov. 
}
\references{
  Tracy, B.F. and M.A. Sanderson. 2000. Patterns of plant species richness in pasture lands of the northeast United States. Plant Ecology 149:169-180.
}
\examples{
data(graze)
}
\keyword{datasets}
