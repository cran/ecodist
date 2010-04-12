\name{crosstab}
\alias{crosstab}
\title{ Data formatting }
\description{
  Converts field data of the form site, species, observation into a matrix.
}
\usage{
crosstab(rowlab, collab, values, type = "sum", allrows, allcols)
}
\arguments{
  \item{rowlab}{ row labels, e.g. site names. }
  \item{collab}{ column labels, e.g. species names. }
  \item{values}{ data values. }
  \item{type}{ function to use to combine data, one of "sum" (default), "min", "max", "mean". }
  \item{allrows}{ optional, list of desired row names that may not appear in rowlab. }
  \item{allcols}{ optional, list of desired column names that may not appear in collab. }
}

\details{
Field data are often recorded as a separate row for each site-species
combination. This function turns them into a matrix for further analysis based on unique row and column labels.

The three vectors should all be the same length (including duplicates).

If allrows or allcols exists, rows and/or columns of zeros are inserted for any elements of allrows/allcols not present in rowlab/collab.
}
\value{
matrix with rowlab as row headings, collab as columns, and
values as the data.
}
\author{ Sarah Goslee, Sarah.Goslee@ars.usda.gov }

\examples{
# Make a random example
plotnames <- sort(rep(1:5, 3))
speciesnames <- rep(c("A", "B", "C"), 5)
freqdata <- runif(15)

data.m <- crosstab(plotnames, speciesnames, freqdata)
}
\keyword{ manip }

