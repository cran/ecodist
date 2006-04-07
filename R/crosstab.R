crosstab <- function(rowlab, collab, values, type="sum")

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

results <- switch(type,
	mean = tapply(values, list(rowlab, collab), mean),
	max = tapply(values, list(rowlab, collab), max),
	min = tapply(values, list(rowlab, collab), min),
	sum = tapply(values, list(rowlab, collab), sum),
)

results[is.na(results)] <- 0

results

}

