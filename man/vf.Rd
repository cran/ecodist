\name{vf}
\alias{vf}
\title{ Vector fitting }
\description{
  Fits environmental variables to an ordination configuration.
}
\usage{
vf(ord, vars, nperm = 100)
}
\arguments{
  \item{ord}{ matrix containing an ordination result with axes as columns. }
  \item{vars}{ matrix with ancillary variables as columns. }
  \item{nperm}{ number of permutation for the significance test. If nperm = 0, the test will be omitted. }
}
\details{
 Vector fitting finds the maximum correlation of the individual variables with aconfiguration of samples in ordination space. 
}
\value{
  a matrix with the first n columns containing the scores for every variable in each of the n dimensions of the ordination space. r is the maximum correlation of the variable with the ordination space, and pval is the result of the permutation test.
  }
\references{ Jongman, R.H.G., C.J.F. ter Braak and O.F.R. van Tongeren. 1995. Data analysis in community and landscape ecology. Cambridge University Press, New York. }
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{plotvf}} }
\examples{

\dontrun{
# Fit vectors for the main variables to an NMDS configuration
iris.vf <- vf(iris.nmin, iris[,1:4])
plotvf(iris.vf, col="blue")
}

# A full example is available in the NMDS 
# section of the main help file for \link{ecodist}.

}
\keyword{ multivariate }
