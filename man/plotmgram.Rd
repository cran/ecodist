\name{plotmgram}
\alias{plotmgram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot a Mantel correlogram }
\description{
  Plot a Mantel correlogram with solid symbols for significant values.
}
\usage{
plotmgram(amgram, pval = 0.05, xlab = "Distance", ylab = "Mantel r", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amgram}{ output from mgram() }
  \item{pval}{ cut-off level for statistical significance. }
  \item{xlab}{ x-axis label. }
  \item{ylab}{ y-axis label. }
  \item{\dots}{ optional, any additional graphics parameters. }
}

\value{
  draws a plot (graphics device must be active).
}
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{mgram}} }
\examples{
\dontrun{
z.mgram <- mgram(z.d, space.d)
plotmgram(z.mgram)
}

# A full example is available in the Mantel correlogram 
# section of the main help file for \link{ecodist}.
}
\keyword{ multivariate }
\keyword{ hplot }
