corgen <- function(len, r, population = FALSE, epsilon = 0)
{
   # generate data vectors of length len with a given correlation value r
   # if population = TRUE, samples are drawn from populations with
   # correlation r, otherwise the correlation of the samples is r.
   #
   # if epsilon = 0, then cor(x, y) is approximately r
   # if epsilon is anything else, then cor(x, y) is within 
   # epsilon of r

   x <- scale(rnorm(len))
   y <- scale(rnorm(len))

   if(!population) {
      # generate uncorrelated data
      xy <- princomp(cbind(x, y))$scores
   }
   
   # impose a correlation
   x <- xy[,1]
   y <- xy[,2]
   a <- r / sqrt(1 - r ^ 2)
   y <- x * a + y

   if(epsilon > 0) {
      while(abs(cor(x, y) - r) > epsilon) {
         x <- scale(rnorm(len))
         y <- scale(rnorm(len))
         # generate uncorrelated data
         xy <- princomp(cbind(x, y))$scores
         
         # impose a correlation
         x <- xy[,1]
         y <- xy[,2]
         a <- r / sqrt(1 - r ^ 2)
         y <- x * a + y
      }
   }

   list(x=x, y=y)
}

