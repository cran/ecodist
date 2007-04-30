\name{corgen}
\alias{corgen}
\title{ Generate correlated data }
\description{
  Generate correlated data vectors of a given length. }
\usage{
corgen(len, r, population = FALSE, epsilon = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{len}{ Length of data vectors. }
  \item{r}{ Correlation between data vectors. }
  \item{population}{ TRUE for vectors drawn from two populations with correlation r, otherwise r is the sample correlation. }
  \item{epsilon}{ If epsilon = 0, it has no effect, otherwise the sampling process is repeated until the sample correlation is within epsilon of r. This option allows the production of exactly-correlated data, within the limits of epsilon. }
}
\value{
  If it is a LIST, use
  \item{x }{First data vector. }
  \item{y }{Second data vector. }
  ...
}
\details{
    This function is more sophisticated that the corgen function available in the stats package.
}
\author{ Sarah Goslee }

}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

# create two random variables of length 100 with correlation 
# of 0.10 +/- 0.00001
corgen(100, .1, epsilon=0.00001)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
