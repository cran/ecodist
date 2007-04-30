plotmgram <- function(amgram, pval = 0.05, xlab = "Distance", ylab = 
                "Mantel r", ...)
{
# amgram is the output from mgram
# pval is the p-value to be considered signficant
# ... are additional graphics parameters

        if(is.list(amgram)) amgram <- amgram$answer.m

        pval.v <- amgram[, 4]
        plot(amgram[, 1], amgram[, 3], type = "l", xlab = xlab, ylab = 
                ylab, ...)
        points(amgram[pval.v <= pval, 1], amgram[pval.v <= pval, 3], pch = 16)
        points(amgram[pval.v > pval, 1], amgram[pval.v > pval, 3], pch = 1)
        invisible()
}
