\name{forestpct}
\alias{forestpct}
\docType{data}
\title{ Forest landcover data}
\description{
  Percentage forest cover within a 1-km radius of 50 sampling locations.
}
\usage{data(forestpct)}
\format{
  The format is:
 num [1:50] 63.9 71.3 72.5 77.1 18.4 ...
}
\details{
  Forest landcover is taken from the USGS 1992 National Land Cover Dataset. All forest classes were combined, and the percentage within 1 km of the sample sites was calculated using a GIS.
}
\source{
  The 1992 NLCD data can be obtained from http://www.mrlc.gov/.
}
\examples{
data(forestpct)
## maybe str(forestpct) ; plot(forestpct) ...
}
\keyword{datasets}
