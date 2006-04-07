#include <S.h>
#include <math.h>

#define RANDIN  seed_in((long *)NULL)
#define RANDOUT seed_out((long *)NULL)
#define UNIF unif_rand()
#define S_EVALUATOR

void bootstrap(double *x, double *y, long *n, long *xlen, long *nboot, double *pboot, double *bootcor, long *rarray, long *rmat, double *xdif, double *ydif)

{

long i, j, k, l;
double r;
double nsamp;
double xmean, ymean;
double xsum;
double xxsum, yysum;

S_EVALUATOR

/* Set random seed using Splus function */

RANDIN;


for(i = 0; i < *nboot; i++) {

/* Set up rarray. */

   for(j = 0; j < *n; j++) {
      r = UNIF;
      if(r > *pboot)
         rarray[j] = 0;
      else rarray[j] = 1;
   }

/* Turn rarray into a lower-triangular sampling matrix. */
/* 1 means include, 0 means omit. */

   l = 0;
   for(j = 1; j < *n; j++) {
      for(k = 0; k < j; k++) {
         if(rarray[j] == 0 || rarray[k] == 0)
            rmat[l] = 0;
         else rmat[l] = 1;
         l++;
      }
   }


   nsamp = 0;
   for(j = 0; j < *xlen; j++) {
      nsamp += rmat[j];
   }


/* Calculate means for x and y. */

   xmean = 0;
   ymean = 0;
   for(j = 0; j < *xlen; j++) {
      if(rmat[j] == 1) {
         xmean += x[j];
         ymean += y[j];
      }
   }
   xmean = xmean/nsamp;
   ymean = ymean/nsamp;

/* Calculate deviations for x and y. */

   for(j = 0; j < *xlen; j++) {
      if(rmat[j] == 1) {
         xdif[j] = x[j] - xmean;
         ydif[j] = y[j] - ymean;
      }
      else {
         xdif[j] = 0;
         ydif[j] = 0;
      }
   }

   xsum = 0;
   xxsum = 0; 
   yysum = 0;

   for(j = 0; j < *xlen; j++) {
      if(rmat[j] == 1) {
         xsum += (xdif[j] * ydif[j]);
         xxsum += (xdif[j] * xdif[j]);
         yysum += (ydif[j] * ydif[j]);
      }
   }

   bootcor[i] = (xsum) / sqrt(xxsum * yysum);

}

/* Reset random seed using an Splus function. */

RANDOUT;

}


void permute(double *x, double *y, long *n, long *xlen, long *nperm, double *zstats, double *tmat, long *rarray)

{

long i, k, l, m;
double cumsum;
long temp;

S_EVALUATOR

/* Set random seed using Splus function */

RANDIN;

/* Calculate first z-statistic (unpermuted data). */

cumsum = 0;

for(k = 0; k < *xlen; k++) {
   cumsum += x[k] * y[k];
}

zstats[0] = cumsum / *xlen;

/* Start permutation routine */

for(i = 1; i < *nperm; i++) {

/* Set up rarray. */

   for(k = 0; k < *n; k++) {
      rarray[k] = k;
   }

/* Convert x to a full matrix. */

   m = 0;
   for(k = 1; k < *n; k++) {
      for(l = 0; l < k; l++) {
         tmat[k * *n + l] = x[m];
         tmat[l * *n + k] = x[m];
         m++;
      }
   }

/* Randomize rarray using an Splus function. */

   for(k = 0; k < (*n - 1); k++) {
      l = *n - k - 1;
      m = (long)((float)l * UNIF);
      if(m > l) m = l;
      temp = rarray[l];
      rarray[l] = rarray[m];
      rarray[m] = temp;
   }

/* Reorder x. */

   m = 0;
   for(k = 1; k < *n; k++) {
      for(l = 0; l < k; l++) {
         x[m] = tmat[rarray[k] * *n + rarray[l]];
         m++;
      }
   }


/* Calculate new sum of products. */

   cumsum = 0;

   for(k = 0; k < *xlen; k++) {
         cumsum += x[k] * y[k];
   
   }

   zstats[i] = cumsum / *xlen;

}

/* Reset random seed using an Splus function. */

RANDOUT;

}



void permpart(double *hmat, double *bmat, double *omat, double *y, double *xcor, double *ycor, long *n, long *ncol, long *xlen, long *nperm, double *zstats, double *tmat, long *rarray)

{

long i, k, l, m;
double cumsum;
double bsum;
double w1, w2;
long temp;

S_EVALUATOR

/* Set random seed using Splus function */

RANDIN;

/* Calculate first z-statistic (unpermuted data). */

cumsum = 0;

for(k = 0; k < *xlen; k++) {
   cumsum += xcor[k] * ycor[k];
}

zstats[0] = cumsum / *xlen;


/* Start permutation routine */

for(i = 1; i < *nperm; i++) {

/* Set up rarray. */

   for(k = 0; k < *n; k++) {
      rarray[k] = k;
   }


/* Convert y to a full matrix. */

   m = 0;
   for(k = 1; k < *n; k++) {
      for(l = 0; l < k; l++) {
         tmat[k * *n + l] = y[m];
         tmat[l * *n + k] = y[m];
         m++;
      }
   }

/* Randomize rarray using an Splus function. */

   for(k = 0; k < (*n - 1); k++) {
      l = *n - k - 1;
      m = (long)((float)l * UNIF);
      if(m > l) m = l;
      temp = rarray[l];
      rarray[l] = rarray[m];
      rarray[m] = temp;
   }


/* Reorder y. */

   m = 0;
   for(k = 1; k < *n; k++) {
      for(l = 0; l < k; l++) {
         y[m] = tmat[rarray[k] * *n + rarray[l]];
         m++;
      }
   }


/* Calculate residuals for y */

/* Calculate bmat */

for(k = 0; k < *ncol; k++) {
   bmat[k] = 0;
}

for(k = 0; k < *ncol; k++) {
   for(l = 0; l < *xlen; l++) {
      bmat[k] = bmat[k] + hmat[l * *ncol + k] * y[l];
   }
}

/* Calculate ycor (residuals) */

for(k = 0; k < *xlen; k++) {
   ycor[k] = 0;
}

for(k = 0; k < *xlen; k++) {
   bsum = 0;
   for(l = 0; l < *ncol; l++) {
      bsum = bsum + bmat[l] * omat[l * *xlen + k];
   }
   ycor[k] = y[k] - bsum;
}


/* Standardize residuals so z = r */

w1 = 0;
w2 = 0;

for(k = 0; k < *xlen; k++) {
   w1 = w1 + ycor[k];
   w2 = w2 + ycor[k] * ycor[k];
}
w1 = w1 / *xlen;
w2 = sqrt(w2 / *xlen - w1 * w1);
for(k = 0; k < *xlen; k++) {
   ycor[k] = (ycor[k] - w1) / w2;
}

/* Calculate new sum of products. */

   cumsum = 0;

   for(k = 0; k < *xlen; k++) {
      cumsum += xcor[k] * ycor[k];
   }

   zstats[i] = cumsum / *xlen;

}

/* Reset random seed using an Splus function. */

RANDOUT;

}


void xbootstrap(double *x, double *y, long *n, long *xlen, long *nboot, double *pboot, double *bootcor, long *rarray, long *rmat, double *xdif, double *ydif)

{

long i, j, k;
double r;
double nsamp;
double xmean, ymean;
double xsum;
double xxsum, yysum;

S_EVALUATOR

/* Set random seed using Splus function */

RANDIN;


for(i = 0; i < *nboot; i++) {

/* Set up rarray. */

   for(j = 0; j < *n; j++) {
      r = UNIF;
      if(r > *pboot)
         rarray[j] = 0;
      else rarray[j] = 1;
   }

/* Turn rarray into a square sampling matrix. */
/* 1 means include, 0 means omit. */

   for(j = 0; j < *xlen; j++) {
      rmat[j] = 1;
   }

   for(j = 0; j < *n; j++) {
      for(k = 0; k <= j; k++) {
         if(rarray[j] == 0 || rarray[k] == 0) {
            rmat[j * *n + k] = 0;
            rmat[k * *n + j] = 0;
         }
      }
   }

   nsamp = 0;
   for(j = 0; j < *xlen; j++) {
      nsamp += rmat[j];
   }


/* Calculate means for x and y. */

   xmean = 0;
   ymean = 0;
   for(j = 0; j < *xlen; j++) {
      if(rmat[j] == 1) {
         xmean += x[j];
         ymean += y[j];
      }
   }
   xmean = xmean/nsamp;
   ymean = ymean/nsamp;

/* Calculate deviations for x and y. */

   for(j = 0; j < *xlen; j++) {
      if(rmat[j] == 1) {
         xdif[j] = x[j] - xmean;
         ydif[j] = y[j] - ymean;
      }
      else {
         xdif[j] = 0;
         ydif[j] = 0;
      }
   }


   xsum = 0;
   xxsum = 0; 
   yysum = 0;

   for(j = 0; j < *xlen; j++) {
      if(rmat[j] == 1) {
         xsum += (xdif[j] * ydif[j]);
         xxsum += (xdif[j] * xdif[j]);
         yysum += (ydif[j] * ydif[j]);
      }
   }

   bootcor[i] = (xsum) / sqrt(xxsum * yysum);

}

/* Reset random seed using an Splus function. */

RANDOUT;

}



void xpermute(double *x, double *y, long *n, long *xlen, long *nperm, double *zstats, double *tmat, long *rarray)

{

long i, k, l, m;
double cumsum;
long temp;

S_EVALUATOR

/* Set random seed using Splus function */

RANDIN;

/* Calculate first z-statistic (unpermuted data). */

cumsum = 0;

for(k = 0; k < *xlen; k++) {
   cumsum += x[k] * y[k];
}

zstats[0] = cumsum;


/* Start permutation routine */

for(i = 1; i < *nperm; i++) {

/* Set up rarray. */

   for(k = 0; k < *n; k++) {
      rarray[k] = k;
   }

/* Randomize rarray using an Splus function. */

   for(k = 0; k < (*n - 1); k++) {
      l = *n - k - 1;
      m = (long)((float)l * UNIF);
      if(m > l) m = l;
      temp = rarray[l];
      rarray[l] = rarray[m];
      rarray[m] = temp;
   }

/* Reorder x. */

   for(k = 0; k < *n; k++) {
      for(l = 0; l <= k; l++) {
         x[k * *n + l] = tmat[rarray[k] * *n + rarray[l]];
         x[l * *n + k] = tmat[rarray[l] * *n + rarray[k]];
      }
   }


/* Calculate new sum of products. */

   cumsum = 0;

   for(k = 0; k < *xlen; k++) {
         cumsum += x[k] * y[k];
   
   }

   zstats[i] = cumsum;

}

/* Reset random seed using an Splus function. */

RANDOUT;

}




void xpermpart(double *hmat, double *y, double *xcor, double *ycor, long *n, long *xlen, long *nperm, double *zstats, double *tmat, long *rarray)

{

long i, k, l, m;
double cumsum;
long temp;

S_EVALUATOR

/* Set random seed using Splus function */

RANDIN;

/* Calculate residuals for y */

for(k = 0; k < *xlen; k++) {
   ycor[k] = 0;
}


for(k = 0; k < *xlen; k++) {
   for(l = 0; l < *xlen; l++) {
      ycor[k] = ycor[k] + hmat[k * *xlen + l] * y[l];
   }
}


/* Calculate first z-statistic (unpermuted data). */

cumsum = 0;

for(k = 0; k < *xlen; k++) {
   cumsum += xcor[k] * ycor[k];
}

zstats[0] = cumsum;


/* Start permutation routine */

for(i = 1; i < *nperm; i++) {

/* Set up rarray. */

   for(k = 0; k < *n; k++) {
      rarray[k] = k;
   }

/* Randomize rarray using an Splus function. */

   for(k = 0; k < (*n - 1); k++) {
      l = *n - k - 1;
      m = (long)((float)l * UNIF);
      if(m > l) m = l;
      temp = rarray[l];
      rarray[l] = rarray[m];
      rarray[m] = temp;
   }

/* Reorder y. */

   for(k = 0; k < *n; k++) {
      for(l = 0; l <= k; l++) {
         y[k * *n + l] = tmat[rarray[k] * *n + rarray[l]];
         y[l * *n + k] = tmat[rarray[l] * *n + rarray[k]];
      }
   }


/* Calculate residuals for y */

for(k = 0; k < *xlen; k++) {
   ycor[k] = 0;
}

for(k = 0; k < *xlen; k++) {
   for(l = 0; l < *xlen; l++) {
      ycor[k] = ycor[k] + hmat[k * *xlen + l] * y[l];
   }
}


/* Calculate new sum of products. */

   cumsum = 0;

   for(k = 0; k < *xlen; k++) {
         cumsum += xcor[k] * ycor[k];
   
   }

   zstats[i] = cumsum;

}

/* Reset random seed using an Splus function. */

RANDOUT;

}



void bcdist(double *x, long *pnrow, long *pncol, double *dist)
{
long i, j, k, l;
long nrow, ncol;
double sumi, sumj;
double minsum;


l = 0;
nrow = *pnrow;
ncol = *pncol;

for(i = 0; i < (nrow - 1); i++) {
   for(j = (i + 1); j < (nrow); j++) {
      minsum = 0;
      sumi = 0;
      sumj = 0;
      for(k = 0; k < ncol; k++) {
         if(x[i * ncol + k] < x[j * ncol + k]) 
            minsum += x[i * ncol + k];
         else 
            minsum += x[j * ncol + k];
         sumi += x[i * ncol + k];
         sumj += x[j * ncol + k]; 
      }
   if((sumi + sumj) == 0) 
      dist[l] = 0;
   else
      dist[l] = (1 - (2 * minsum) / (sumi + sumj));
   l++;
   }
}
}


void weight(long *n, double *datadist, double *d1, double *d2, double *w)

{
long i;
double m1, m2;
double w1, w2;
double pi;

pi = 2 * acos(0);

for(i = 0; i < *n * *n; i++) {
   if(datadist[i] != 0) {
      if(d1[i] < datadist[i])
         m1 = d1[i] / datadist[i];
      else m1 = 1;

      if(d2[i] < datadist[i])
         m2 = d2[i] / datadist[i];
      else m2 = 1;
   }
   else {
      m1 = 0;
      m2 = 0;
   }

   w1 = 1 - (acos(m1) + acos(m2)) / pi;

   if(datadist[i] != 0) {
      m1 = d1[i] / datadist[i];
      if(m1 > 1)
         m1 = 1;

      m2 = d2[i] / datadist[i];
      if(m2 > 1)
         m2 = 1;
   }
   else {
      m1 = 0;
      m2 = 0;
   }

   w2 = 0.75 - (acos(m1) + acos(m2)) / (2 * pi);

   if((datadist[i] * datadist[i]) >= (d1[i] * d1[i] + d2[i] * d2[i]))
      w1 = 0;
   if((datadist[i] * datadist[i]) < (d1[i] * d1[i] + d2[i] * d2[i]))
      w2 = 0;

   w[i] = w1 + w2;
}
}

void newpermone(double *x, long *dclass, long *n, long *xlen, long *nperm, double *zstats, double *tmat, long *rarray)

{

long i, k, l, m;
double cumsum;
long temp;

S_EVALUATOR

/* Set random seed using Splus function */

RANDIN;

/* Calculate first z-statistic (unpermuted data). */

cumsum = 0;

for(k = 0; k < *xlen; k++) {
	if(dclass[k] == 0) {
   	cumsum += x[k];
	}
}

zstats[0] = cumsum;

/* Start permutation routine */

for(i = 1; i < *nperm; i++) {

/* Set up rarray. */

   for(k = 0; k < *n; k++) {
      rarray[k] = k;
   }

/* Convert x to a full matrix. */

   m = 0;
   for(k = 1; k < *n; k++) {
      for(l = 0; l < k; l++) {
         tmat[k * *n + l] = x[m];
         tmat[l * *n + k] = x[m];
         m++;
      }
   }

/* Randomize rarray using an Splus function. */

   for(k = 0; k < (*n - 1); k++) {
      l = *n - k - 1;
      m = (long)((float)l * UNIF);
      if(m > l) m = l;
      temp = rarray[l];
      rarray[l] = rarray[m];
      rarray[m] = temp;
   }

/* Reorder x. */

   m = 0;
   for(k = 1; k < *n; k++) {
      for(l = 0; l < k; l++) {
         x[m] = tmat[rarray[k] * *n + rarray[l]];
         m++;
      }
   }


/* Calculate new sum of products. */

   cumsum = 0;

   for(k = 0; k < *xlen; k++) {
      if(dclass[k] == 0) {
          cumsum += x[k];
      }
   }

   zstats[i] = cumsum;

}

/* Reset random seed using an Splus function. */

RANDOUT;

}




void newpermtwo(double *x, double *y, long *n, long *xlen, long *nperm, double *zstats, double *tmat, long *rarray)

{

long i, k, l, m;
double cumsum;
long temp;
float naval = -9999;

S_EVALUATOR

/* Set random seed using Splus function */

RANDIN;

/* Calculate first z-statistic (unpermuted data). */

cumsum = 0;

for(k = 0; k < *xlen; k++) {
	if(x[k] != naval) {
	   cumsum += x[k] * y[k];
	}
}

zstats[0] = cumsum;

/* Start permutation routine */

for(i = 1; i < *nperm; i++) {

/* Set up rarray. */

   for(k = 0; k < *n; k++) {
      rarray[k] = k;
   }

/* Convert x to a full matrix. */

   m = 0;
   for(k = 1; k < *n; k++) {
      for(l = 0; l < k; l++) {
         tmat[k * *n + l] = x[m];
         tmat[l * *n + k] = x[m];
         m++;
      }
   }

/* Randomize rarray using an Splus function. */

   for(k = 0; k < (*n - 1); k++) {
      l = *n - k - 1;
      m = (long)((float)l * UNIF);
      if(m > l) m = l;
      temp = rarray[l];
      rarray[l] = rarray[m];
      rarray[m] = temp;
   }

/* Reorder x. */

   m = 0;
   for(k = 1; k < *n; k++) {
      for(l = 0; l < k; l++) {
         x[m] = tmat[rarray[k] * *n + rarray[l]];
         m++;
      }
   }


/* Calculate new sum of products. */

   cumsum = 0;

   for(k = 0; k < *xlen; k++) {
      if(x[k] != naval) {
          cumsum += x[k] * y[k];
      }
   }
	
   zstats[i] = cumsum;

}

/* Reset random seed using an Splus function. */

RANDOUT;

}




void psum(double *x, long *pnrow, long *pncol, double *dist)
{
long row1, row2, col1;
long nrow, ncol;
long l;
double thisval, thatval;

l = 0;
nrow = *pnrow;
ncol = *pncol;

  for(col1 = 0; col1 < ncol; col1++) {
	for(row1 = 0; row1 < nrow; row1++) {
		thatval = x[row1 * ncol + col1];
		for(row2 = 0; row2 < nrow; row2++) {
			thisval = x[row2 * ncol + col1];
			dist[l] = thisval + thatval;
			l++;
		}
	}
}

}
    

void pdiff(double *x, long *pnrow, long *pncol, double *dist)
{
long row1, row2, col1;
long nrow, ncol;
long l;
double thisval, thatval;

l = 0;
nrow = *pnrow;
ncol = *pncol;

  for(col1 = 0; col1 < ncol; col1++) {
	for(row1 = 0; row1 < nrow; row1++) {
		thatval = x[row1 * ncol + col1];
		for(row2 = 0; row2 < nrow; row2++) {
			thisval = x[row2 * ncol + col1];
			dist[l] = thisval - thatval;
			l++;
		}
	}
}

}
    

void jpres(double *x, long *pnrow, long *pncol, double *dist)
{
long row1, row2, col1;
long nrow, ncol;
long l;
double thisval, thatval;

l = 0;
nrow = *pnrow;
ncol = *pncol;

  for(col1 = 0; col1 < ncol; col1++) {
	for(row1 = 0; row1 < nrow; row1++) {
		thatval = x[row1 * ncol + col1];
		for(row2 = 0; row2 < nrow; row2++) {
			thisval = x[row2 * ncol + col1];
			if(thisval > 0 & thatval > 0) {
				dist[l] = 1;
			}
			else {
				dist[l] = 0;
			}
			l++;
		}
	}
}

}
    

void jabs(double *x, long *pnrow, long *pncol, double *dist)
{
long row1, row2, col1;
long nrow, ncol;
long l;
double thisval, thatval;

l = 0;
nrow = *pnrow;
ncol = *pncol;

  for(col1 = 0; col1 < ncol; col1++) {
	for(row1 = 0; row1 < nrow; row1++) {
		thatval = x[row1 * ncol + col1];
		for(row2 = 0; row2 < nrow; row2++) {
			thisval = x[row2 * ncol + col1];
			if(thisval == 0 & thatval == 0) {
				dist[l] = 1;
			}
			else {
				dist[l] = 0;
			}
			l++;
		}
	}
}

}
    

void jfirst(double *x, long *pnrow, long *pncol, double *dist)
{
long row1, row2, col1;
long nrow, ncol;
long l;
double thisval, thatval;

l = 0;
nrow = *pnrow;
ncol = *pncol;

  for(col1 = 0; col1 < ncol; col1++) {
	for(row1 = 0; row1 < nrow; row1++) {
		thatval = x[row1 * ncol + col1];
		for(row2 = 0; row2 < nrow; row2++) {
			thisval = x[row2 * ncol + col1];
			if(thisval > 0 & thatval == 0) {
				dist[l] = 1;
			}
			else {
				dist[l] = 0;
			}
			l++;
		}
	}
}

}
    

void jsec(double *x, long *pnrow, long *pncol, double *dist)
{
long row1, row2, col1;
long nrow, ncol;
long l;
double thisval, thatval;

l = 0;
nrow = *pnrow;
ncol = *pncol;

  for(col1 = 0; col1 < ncol; col1++) {
	for(row1 = 0; row1 < nrow; row1++) {
		thatval = x[row1 * ncol + col1];
		for(row2 = 0; row2 < nrow; row2++) {
			thisval = x[row2 * ncol + col1];
			if(thisval == 0 & thatval > 0) {
				dist[l] = 1;
			}
			else {
				dist[l] = 0;
			}
			l++;
		}
	}
}

}

