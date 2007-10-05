corgen <- function (len, x, r, population = FALSE, epsilon = 0) 
{
   # generate data vectors of length len with a given correlation value r
   # if population = TRUE, samples are drawn from populations with
   # correlation r, otherwise the correlation of the samples is r.
   #
   # if epsilon = 0, then cor(x, y) is approximately r
   # if epsilon is anything else, then cor(x, y) is within 
   # epsilon of r

    rsign <- sign(r)
    if(rsign == 0) rsign <- 1
    r <- abs(r)
    if(missing(x)) {
        if(missing(len)) {
            stop("Must specify x or len.\n")
        }
        else {
            x <- scale(rnorm(len))
            x.orig <- x
        }
    }
    else {
        len <- length(x)
        x.orig <- x
        x <- scale(x)
    }
    y <- scale(rnorm(len))
    if (!population) {
        # generate uncorrelated data
        xy <- princomp(cbind(x, y))$scores
   
        # impose a correlation
        x <- xy[, 1]
        x.orig <- x
        y <- xy[, 2]
    }
    a <- r/sqrt(1 - r^2)
    y <- x * a + y
    if (epsilon > 0) {
        while (abs(cor(x, y) - r) > epsilon) {

            y <- scale(rnorm(len))
            if (!population) {
                xy <- princomp(cbind(x, y))$scores
                x <- xy[, 1]
                x.orig <- x
                y <- xy[, 2]
            }
            a <- r/sqrt(1 - r^2)
            y <- x * a + y

        }
    }
    list(x = x.orig, y = y*rsign)
}

