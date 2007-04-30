nmds.min <- function(x, dims=0)

{
# returns the minimum-stress configuration from nmds output.
# (Results from nmds)
# if dims==0, returns the overall lowest-stress configuration 
# Otherwise, returns the lowest-stress configuration of dimensionality dims

    if(dims == 0) {
        x.min <- x$conf[x$stress == min(x$stress)]
    }
    else {
        x.dims <- sapply(x$conf, ncol)
        x$conf <- x$conf[x.dims == dims]
        x$stress <- x$stress[x.dims == dims]
        x.min <- x$conf[x$stress == min(x$stress)]
    }
    x.min <- x.min[[1]]
    data.frame(x.min)
}
