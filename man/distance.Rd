\name{distance}
\alias{distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate dissimilarity/distance metrics }
\description{
  This function calculates a variety of dissimilarity or distance metrics. Although it duplicates the functionality of dist() and bcdist(), it is written entirely in R, and can easily be extended to include additional metrics.
}
\usage{
distance(x, method = "euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ matrix or data frame with rows as samples and columns as variables (such as 
species). Distances will be calculated for each pair of rows. }
  \item{method}{ Currently 5 dissimilarity metrics can be calculated: 
euclidean
bray-curtis
manhattan
mahalanobis
jaccard
Partial matching will work for selecting a method.
  }
}

\value{
  Returns a lower-triangular distance matrix as an object of class "dist".
  
}
\author{ Sarah Goslee <Sarah.Goslee@ars.usda.gov> }

\seealso{ \code{\link{dist}} }
\examples{
\dontrun{
data(iris)
iris.md <- distance(iris[,1:4], "mahal")
}

# A full example is available in the distance 
# section of the main help file for \link{ecodist}.
}
\keyword{ multivariate }
