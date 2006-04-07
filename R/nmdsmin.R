nmds.min <- function(x)

{
# returns the minimum-stress configuration from nmds output.
# (Results from nmds)

x.min <- x$conf[x$stress == min(x$stress)]
x.min <- x.min[[1]]
data.frame(x.min)
}
