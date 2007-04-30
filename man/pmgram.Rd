\name{pmgram}
\alias{pmgram}
\title{ Partial Mantel correlogram }
\description{
 This function calculates simple and partial multivariate correlograms. 
}
\usage{
pmgram(data, space, partial, breaks, nclass, stepsize, resids = FALSE, nperm = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ lower-triangular dissimilarity matrix. }
  \item{space}{ lower-triangular matrix of geographic distances. }
  \item{partial}{ optional, lower-triangular dissimilarity matrix of ancillary data. }
  \item{breaks}{ locations of class breaks. If specified, overrides nclass and stepsize. }
  \item{nclass}{ number of distance classes. If not specified, Sturge's rule will be used
to determine an appropriate number of classes. }
  \item{stepsize}{ width of each distance class. If not specified, nclass and the range of space.d will be used to calculate an appropriate default. }
  \item{resids}{ if resids=TRUE, will return the residuals for each distance class. }
  \item{nperm}{ number of permutations to use. If set to 0, the permutation test will be omitted. }
}
\details{
This function does four different analyses:
	If data has 1 column and partial is missing, 
	calculates a multivariate correlogram for data.

If data has 2 columns and partial is missing,
	calculates Mantel r between the two columns for each distance class.

If data has 1 column and partial exists,
	does a multivariate correlogram for the residuals,
	taking residuals over whole extent.

If data has 2 columns and partial exists,
	does a partial multivariate correlogram, 
	calculating partial for each distance class separately.

}
\value{
  Returns a matrix with one row for each distance class and 4 columns:
  \item{lag }{midpoint of the distance class.}
  \item{ngroup }{number of distances in that class.}
  \item{mantelr }{Mantel r value.}
  \item{pval }{two-sided p-value.}
}
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{mgram}}, \code{\link{mantel}} }
\examples{

\dontrun{
# take partial correlogram of z2 on the residuals of z1 ~ space.d
z.pmgram <- pmgram(z1.d, space.d, z2.d, nperm=10)
plotmgram(z.pmgram, pval=0.1)
}

# A full example is available in the partial Mantel 
# correlogram section of the main help file for \link{ecodist}.

}
\keyword{ multivariate }
