lower <- function(m)
{
# Takes the lower triangle of a matrix
# Does NOT check for symmetric matrix

if(ncol(m) != nrow(m))
	stop("Matrix not square.")

	m[col(m) < row(m)]
}
