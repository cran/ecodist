\name{mgram}
\alias{mgram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mantel correlogram }
\description{
	Calculates simple and partial Mantel correlograms. 
}
\usage{
mgram(species.d, space.d, breaks, nclass, stepsize, nperm = 1000, mrank = FALSE, nboot = 500, pboot = 0.9, cboot = 0.95, alternative = "two.sided", trace = FALSE)
}
\arguments{
  \item{species.d}{ lower-triangular dissimilarity matrix. }
  \item{space.d}{ lower-triangular matrix of geographic distances. }
  \item{breaks}{ locations of class breaks. If specified, overrides nclass and stepsize. }
  \item{nclass}{ number of distance classes. If not specified, Sturge's rule will be used
to determine an appropriate number of classes. }
  \item{stepsize}{ width of each distance class. If not specified, nclass and the range of space.d will be used to calculate an appropriate default. }
  \item{nperm}{ number of permutations to use. If set to 0, the permutation test will be omitted. }
  \item{mrank}{ if this is set to F (the default option), Pearson correlations will be used. If
 set to T, the Spearman correlation (correlation ranked distances) will be used. }
  \item{nboot}{ number of iterations to use for the bootstrapped confidence limits. If set to 0,
 the bootstrapping will be omitted. }
  \item{pboot}{ the level at which to resample the data for the bootstrapping procedure. }
  \item{cboot}{ the level of the confidence limits to estimate. }
  \item{alternative}{ default is "two.sided", and returns p-values for H0: rM = 0. The alternative is "one.sided", which returns p-values for H0: rM <= 0. }
  \item{trace}{ if TRUE, returns progress indicators. }
}
\details{
This function calculates Mantel correlograms. The Mantel correlogram is essentially a multivariate autocorrelation function.
The Mantel r represents the dissimilarity in variable composition  (often 
species composition) at a particular lag distance.
}
\value{
	Returns a matrix with one row for each distance class and 6 columns:
  \item{lag }{midpoint of the distance class.}
  \item{ngroup }{number of distances in that class.}
  \item{mantelr }{Mantel r value.}
  \item{pval }{p-value for the test chosen.}
  \item{llim }{lower bound of confidence limit for mantelr.}
  \item{ulim }{upper bound of confidence limit for mantelr.}
}
\references{ Legendre, P. and M. Fortin. 1989. Spatial pattern and ecological analysis. 
Vegetatio 80:107-138. }
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\seealso{ \code{\link{mantel}}, \code{\link{plotmgram}}, \code{\link{pmgram}} }
\examples{
\dontrun{
z.mgram <- mgram(z.d, space.d)
plotmgram(z.mgram)
}

# A full example is available in the Mantel correlogram 
# section of the main help file for \link{ecodist}.

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
