\name{graze}
\alias{graze}
\docType{data}
\title{ Grazed vegetation data. }
\description{
  These data contain the plant community compositions of 50 grazed pastures in the northeastern United States. Elements are the mean values for canopy cover for ten 0.5 x 2 m quadrats.
}
\usage{data(graze)}
\format{
  A data frame with 50 observations on 115 species variables.
}
\details{
  Species names follow the USDA Plants database.
}
\source{
    Details of these data are available in Tracy and Sanderson (2000) and Goslee, Sanderson and Tracy (ms in review).
  }
\references{
  Species codes are from http://plants.usda.gov. 
  Tracy, B.F. and M.A. Sanderson. 2000. Patterns of plant species richness in pasture lands of the northeast United States. Plant Ecology 149:169-180.
}
\examples{
data(graze)
## maybe str(graze) ; plot(graze) ...
}
\keyword{datasets}
