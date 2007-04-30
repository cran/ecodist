\name{permresults}
\alias{permresults}
\docType{data}
\title{ Mantel test results for number of permutations, correlation and dataset size.}
\description{
  These results explore the variability in randomization tests of significance and in bootstrap confidence limits for different sample sizes, strengths of relationship and numbers of permutations.
}
\usage{data(permresults)}
\format{
  A data frame with 15000 observations on the following 17 variables.
  \describe{
    \item{\code{n}}{The number of samples.}
    \item{\code{r}}{The Mantel r of the dissimilarity vectors.}
    \item{\code{p500}}{The p-values generated using 500 permutations.}
    \item{\code{l100}}{The lower confidence limit generated using 100 subsamples.}
    \item{\code{u100}}{The upper confidence limit generated using 100 subsamples.}
    \item{\code{p1000}}{p-values with 1,000 permutations}
    \item{\code{l500}}{lower CI with 500 permutations}
    \item{\code{u500}}{upper CI with 500 permutations}
    \item{\code{p10000}}{p-values with 10,000 permutations}
    \item{\code{l1000}}{lower CI with 1,000 permutations}
    \item{\code{u1000}}{upper CI with 1,000 permutations}
    \item{\code{p100000}}{p-values with 100,000 permutations}
    \item{\code{l10000}}{lower CI with 10,000 permutations}
    \item{\code{u10000}}{upper CI with 10,000 permutations}
    \item{\code{p1000000}}{p-values with 1,000,000 permutations}
    \item{\code{l100000}}{lower CI with 100,000 permutations}
    \item{\code{u100000}}{upper CI with 100,000 permutations}
  }
}
\examples{
data(permresults)
}
\keyword{datasets}
