plotvf <- function(x, cex=.8, ascale=.9, ...)
{
# add vector fitting arrows to an existing ordination plot
# x is output from vf
# other arguments are from arrows command and graphics parameters.

plotlim <- par()$usr
plotlim <- min((plotlim[2]-plotlim[1]), (plotlim[4]-plotlim[3]))
ascale <- ascale * (plotlim/2)


    for (i in 1:dim(x)[[1]]) {
        arrows(0, 0, x[i, 1] * x[i,3] * ascale, x[i, 2] * x[i,3] * ascale, ...)
    }
    text(x[, 1] * x[,3] * (ascale+.1), x[, 2] * x[,3]*(ascale+.1), dimnames(x)[[1]], cex = cex, ...)
}

