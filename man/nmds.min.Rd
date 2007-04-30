\name{nmds.min}
\alias{nmds.min}
\title{ Find minimum stress configuration }
\description{
  Finds minimum stress configuration from output of nmds()
}
\usage{
nmds.min(x, dims = 0)
}
\arguments{
  \item{x}{ output from nmds() }
  \item{dims}{ desired dimensionality of result. If dims == 0 (the default) then the overall minimum stress configuration is chosen. }
}
\value{
  Configuration of minimum stress ordination (dataframe of
coordinates).
}
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{nmds}} }
\examples{

\dontrun{
# Example of multivariate analysis using built-in iris dataset
data(iris)
iris.md <- distance(iris[,1:4], "mahal")

# Minimum-stress 2-dimensional nonmetric multidimensional scaling configuration
iris.nmds <- nmds(iris.md)
iris.nmin <- nmds.min(iris.nmds)

# Plot NMDS result with symbols denoting species
plot(iris.nmin, pch=as.numeric(iris[,5]))
}

# A full example is available in the NMDS 
# section of the main help file for \link{ecodist}.
}
\keyword{ multivariate }
