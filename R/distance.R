distance <- function(x, method="euclidean")
{
# calculates similarity and dissimilarity coefficients
# as described in Legendre and Legendre 1998
# returns lower-triangle
###
# Sarah Goslee
# 2 March 2006
# revised 31 March 2008
### 
# uses clever matrix math to calculate the pieces needed
# by dissimilarity matrices, to make it easy to add new
# indices.
###
# to add a new metric:
# add it to the commented list below
# add it to the end of the METHODS <- c(...) list
# add the code at the appropriate point at the bottom of
# the function

# 1: euclidean
# 2: bray-curtis
# 3: manhattan
# 4: mahalanobis
# 5: jaccard
# 6: simple difference
# 7: sorensen

pairedsum <- function(x)
{
	### paired sums
   ### returns an N by N by P matrix containing each
   ### combination of 
    N <- nrow(x)
	 P <- ncol(x)
	 A <- numeric(N * N * P)
    A <- .C("psum",
		as.double(as.vector(t(x))),
		as.integer(N),
		as.integer(P),
		A = as.double(A),
		PACKAGE = "ecodist")$A
	 
	 A <- array(A, dim=c(N, N, P))
    A
}

paireddiff <- function(x)
{
### paired differences
    N <- nrow(x)
	 P <- ncol(x)
	 A <- numeric(N * N * P)
    A <- .C("pdiff",
		as.double(as.vector(t(x))),
		as.integer(N),
		as.integer(P),
		A = as.double(A),
		PACKAGE = "ecodist")$A
	 
	 A <- array(A, dim=c(N, N, P))
    A
}

jointpresence <- function(x)
{
### joint count of presences
    N <- nrow(x)
	 P <- ncol(x)
	 A <- numeric(N * N * P)
    A <- .C("jpres",
		as.double(as.vector(t(x))),
		as.integer(N),
		as.integer(P),
		A = as.double(A),
		PACKAGE = "ecodist")$A
	 
	 A <- array(A, dim=c(N, N, P))
    A
}

jointabsence <- function(x)
{
### joint count of absences
    N <- nrow(x)
	 P <- ncol(x)
	 A <- numeric(N * N * P)
    A <- .C("jabs",
		as.double(as.vector(t(x))),
		as.integer(N),
		as.integer(P),
		A = as.double(A),
		PACKAGE = "ecodist")$A
	 
	 A <- array(A, dim=c(N, N, P))
    A
}

firstonly <- function(x)
{
### present only in first sample
    N <- nrow(x)
	 P <- ncol(x)
	 A <- numeric(N * N * P)
    A <- .C("jfirst",
		as.double(as.vector(t(x))),
		as.integer(N),
		as.integer(P),
		A = as.double(A),
		PACKAGE = "ecodist")$A
	 
	 A <- array(A, dim=c(N, N, P))
    A
}

secondonly <- function(x)
{
### present only in second sample
    N <- nrow(x)
	 P <- ncol(x)
	 A <- numeric(N * N * P)
    A <- .C("jsec",
		as.double(as.vector(t(x))),
		as.integer(N),
		as.integer(P),
		A = as.double(A),
		PACKAGE = "ecodist")$A
	 
	 A <- array(A, dim=c(N, N, P))
    A
}

x <- as.matrix(x)

## code borrowed from dist()
    METHODS <- c("euclidean", "bray-curtis", "manhattan", "mahalanobis", "jaccard", "difference", "sorensen")

    method <- pmatch(method, METHODS)
    if (is.na(method)) 
        stop("invalid distance method")
    if (method == -1) 
        stop("ambiguous distance method")
    N <- nrow(x)
	 P <- ncol(x)


if(method == 1) 
{
# Euclidean distance
   A <- paireddiff(x)
	D <- sqrt(apply(A, 1:2, function(x)sum(x * x)))
}

if(method == 2)
{
# Bray-Curtis distance
    A <- paireddiff(x)
	 A <- apply(A, 1:2, function(x)sum(abs(x)))
	 B <- pairedsum(x)
	 B <- apply(B, 1:2, sum)
    D <- A / B
}

if(method == 3)
{
# unstandardized manhattan distance
A <- paireddiff(x)
D <- apply(A, 1:2, function(x)sum(abs(x)))
}

if(method == 4)
{
# pairwise Mahalanobis distance
# same as mahal()
	icov <- solve(cov(x))
	A <- paireddiff(x)
	A1 <- apply(A, 1, function(z)(z %*% icov %*% t(z)))
	D <- A1[seq(1, N*N, by=(N+1)), ]
}

if(method == 5)
{
# Jaccard distance
    A <- jointpresence(x)
	 A <- apply(A, 1:2, sum)
	 B <- firstonly(x)
	 B <- apply(B, 1:2, sum)
	 C <- secondonly(x)
	 C <- apply(C, 1:2, sum)
	 D <- 1 - A / (A + B + C)
}

if(method == 6)
{
# simple difference, NOT symmetric
D <- paireddiff(x)[,,1, drop=TRUE]
}

if(method == 6)
{
# Sorensen distance
    A <- jointpresence(x)
	 A <- apply(A, 1:2, sum)
	 B <- firstonly(x)
	 B <- apply(B, 1:2, sum)
	 C <- secondonly(x)
	 C <- apply(C, 1:2, sum)
	 D <- 1 - (2*A) / (2*A + B + C)
}
	 

## Make the results lower triangular	 
    D <- D[col(D) < row(D)]

## give the results attributes similar to dist()
    attr(D, "Size") <- N
	 attr(D, "Labels") <- rownames(x)
	 attr(D, "Diag") <- FALSE
	 attr(D, "Upper") <- FALSE
    attr(D, "method") <- METHODS[method]
    class(D) <- "dist"
    D
}
