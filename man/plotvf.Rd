\name{plotvf}
\alias{plotvf}
\title{ Plots fitted vectors onto an ordination diagram }
\description{
Add vector fitting arrows to an existing ordination plot.
}
\usage{
plotvf(x, pval = 1, cex = 0.8, ascale = 0.9, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ output from vf() }
  \item{pval}{ optional, critical p-value for choosing variables to plot }
  \item{cex}{ text size }
  \item{ascale}{ optional, proportion of plot area to use when calculating arrow length }
  \item{\dots}{ optional, other graphics parameters }
}
\value{
  Adds arrows to an existing ordination plot. Only arrows with a p-value less than pval are added. By default, all variables are shown.
  }
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{vf}}}
\examples{

\dontrun{

# Fit vectors for the main variables to an NMDS configuration
iris.vf <- vf(iris.nmin, iris[,1:4])
plotvf(iris.vf, pval = 0.05, col="blue")
}

# A full example is available in the NMDS 
# section of the main help file for \link{ecodist}.

}
\keyword{ multivariate }
\keyword{ hplot }
