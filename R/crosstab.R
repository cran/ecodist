crosstab <- function(rowlab, collab, values, type="sum", allrows, allcols)

{

# Converts field data in the form:
# site, species, observation
# into a site by species table
#
# rowlab can be a matrix with three columns
# or each element can be specified individually.
#
# By default, takes the sum of the data
# Can also use "mean", "max" or "min"
# Sarah C. Goslee
# 24 Sept 2003
#
###
#
# if allrows or allcols exists, will expand the matrix to
# include all of those elements into the rows or colums of
# the results matrix
# scg 20 Apr 2006

rowlab <- as.vector(rowlab)
collab <- as.vector(collab)

results <- switch(type,
	mean = tapply(values, list(rowlab, collab), mean),
	max = tapply(values, list(rowlab, collab), max),
	min = tapply(values, list(rowlab, collab), min),
	sum = tapply(values, list(rowlab, collab), sum),
)

if(!missing(allrows)) {
allrows <- as.vector(allrows)
allrows <- c(rownames(results), allrows)
allrows <- sort(unique(allrows))
newrows <- allrows[!(allrows %in% rownames(results))]
temp <- matrix(NA, ncol=ncol(results), nrow=length(newrows))
colnames(temp) <- colnames(results)
rownames(temp) <- newrows
results <- rbind(results, temp)
results <-results[order(rownames(results)),]
}

if(!missing(allcols)) {
allcols <- as.vector(allcols)
allcols <- c(colnames(results), allcols)
allcols <- sort(unique(allcols))
newcols <- allcols[!(allcols %in% colnames(results))]
temp <- matrix(NA, nrow=nrow(results), ncol=length(newcols))
rownames(temp) <- rownames(results)
colnames(temp) <- newcols
results <- cbind(results, temp)
results <-results[,order(colnames(results))]
}

results[is.na(results)] <- 0

results

}

