\name{addord}
\alias{addord}
\title{
Fit new points to an existing NMDS configuration.
}
\description{
Uses a brute force algorithm to find the location for each new point that minimizes overall stress.
}
\usage{
addord(origconf, fulldat, fulldist, isTrain, bfstep = 10, maxit = 50, epsilon = 1e-12)
}

\arguments{
  \item{origconf}{
    The original ordination configuration.
}
  \item{fulldat}{
    The dataset containing original and new points.
}
  \item{fulldist}{
    A dissimilarity matrix calculated on \code{fulldat}.
}
  \item{isTrain}{
    A boolean vector of length \code{nrow(fulldat)} indicating which rows were training data used in determining \code{origconf} (TRUE), or are new points (FALSE).
}
  \item{bfstep}{
    A tuning parameter for the brute force algorithm describing the size of grid to use.
}
  \item{maxit}{
    The maximum number of iterations to use.
}
  \item{epsilon}{
	Tolerance value for convergence.
}
}
\details{
   A region comprising the original ordination configuration plus one standard deviation is divided into a grid of \code{bfstep} rows and columns. For a new point, the grid cell with the lowest stress is identified. That cell is divided into a finer grid, and the lowest-stress cell identified. This process is repeated up to \code{maxit} times, or until stress changes less than \code{epsilon}.
}
\value{
  \item{fullfitconf }{The new ordination configuration containing training and new points.}
  \item{stress }{The stress value for each point.}
  \item{isTrain }{The boolean vector indicating training set membership, for reference.}
}
\author{ Sarah Goslee }

\examples{
data(iris)
iris.d <- dist(iris[,1:4])

### nmds() is timeconsuming, so this was generated
### in advance and saved.
### set.seed(1234)
### iris.nmds <- nmds(iris.d, nits=20, mindim=1, maxdim=4)
### save(iris.nmds, file="ecodist/data/iris.nmds.rda")
data(iris.nmds)

# examine fit by number of dimensions
plot(iris.nmds)

# choose the best two-dimensional solution to work with
iris.nmin <- min(iris.nmds, dims=2)

# rotate the configuration to maximize variance
iris.rot <- princomp(iris.nmin)$scores

# rotation preserves distance apart in ordination space
cor(dist(iris.nmin), dist(iris.rot))

# fit the data to the ordination as vectors
### vf() is timeconsuming, so this was generated
### in advance and saved.
### set.seed(1234)
### iris.vf <- vf(iris.nmin, iris[,1:4], nperm=1000)
### save(iris.vf, file="ecodist/data/iris.vf.rda")
data(iris.vf)

# repeat for the rotated ordination
### vf() is timeconsuming, so this was generated
### in advance and saved.
### set.seed(1234)
### iris.vfrot <- vf(iris.rot, iris[,1:4], nperm=1000)
### save(iris.vfrot, file="ecodist/data/iris.vfrot.rda")
data(iris.vfrot)

par(mfrow=c(1,2))
plot(iris.nmin, col=as.numeric(iris$Species), pch=as.numeric(iris$Species), main="NMDS")
plot(iris.vf)
plot(iris.rot, col=as.numeric(iris$Species), pch=as.numeric(iris$Species), main="Rotated NMDS")
plot(iris.vfrot)

####### addord example

# generate new data points to add to the ordination
# this might be new samples, or a second dataset

iris.new <- structure(list(Sepal.Length = c(4.6, 4.9, 5.4, 5.2, 6, 6.5, 6, 
6.8, 7.3), Sepal.Width = c(3.2, 3.5, 3.6, 2.3, 2.8, 3, 2.7, 3.1, 
3.2), Petal.Length = c(1.2, 1.5, 1.5, 3.5, 4.1, 4.2, 4.8, 5, 
5.7), Petal.Width = c(0.26, 0.26, 0.26, 1.2, 1.3, 1.4, 1.8, 2, 
2), Species = structure(c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L), .Label = c("setosa", 
"versicolor", "virginica"), class = "factor")), .Names = c("Sepal.Length", 
"Sepal.Width", "Petal.Length", "Petal.Width", "Species"), class = "data.frame", row.names = c(NA, 
-9L))

# provide a dist object containing original and new data
# provide a logical vector indicating which samples were used to
# construct the original configuration

iris.full <- rbind(iris, iris.new)
all.d <- dist(iris.full[,1:4])
is.orig <- c(rep(TRUE, nrow(iris)), rep(FALSE, nrow(iris.new)))

### addord() is timeconsuming, so this was generated
### in advance and saved.
### set.seed(1234)
### iris.fit <- addord(iris.nmin, iris.full[,1:4], all.d, is.orig, maxit=100)
### save(iris.fit, file="ecodist/data/iris.fit.rda")
data(iris.fit)

plot(iris.fit$conf, col=iris.full$Species, pch=c(18, 4)[is.orig + 1], xlab="NMDS 1", ylab="NMDS 2")
title("Demo: adding points to an ordination")
legend("bottomleft", c("Training set", "Added point"), pch=c(4, 18))
legend("topright", levels(iris$Species), fill=1:3)

}

\keyword{ multivariate }
