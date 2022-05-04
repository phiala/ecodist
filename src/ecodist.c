#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h> /* for dgemm */

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

void bootstrap(double *x, double *y, int *n, int *xlen, int *nboot, double *pboot, double *bootcor, int *rarray, int *rmat, double *xdif, double *ydif)

{

int i, j, k, l;
double r;
double nsamp;
double xmean, ymean;
double xsum;
double xxsum, yysum;


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


void permute(double *x, double *y, int *n, int *xlen, int *nperm, double *zstats, double *tmat, int *rarray)

{

int i, k, l, m;
double cumsum;
int temp;


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
      m = (int)((float)l * UNIF);
      if(m > l) m = l;
      temp = rarray[l];
      rarray[l] = rarray[m];
      rarray[m] = temp;
   }


/* Reorder x and take lower triangle. */

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



void permpart(double *hmat, double *bmat, double *omat, double *y, double *xcor, double *ycor, int *n, int *ncol, int *xlen, int *nperm, double *zstats, double *tmat, int *rarray)

{

int i, k, l, m;
double cumsum;
double bsum;
double w1, w2;
int temp;


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
      m = (int)((float)l * UNIF);
      if(m > l) m = l;
      temp = rarray[l];
      rarray[l] = rarray[m];
      rarray[m] = temp;
   }


/* Reorder y and take lower triangle. */

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



void bcdistc(double *x, int *pnrow, int *pncol, double *dist)
{
int i, j, k, l;
int nrow, ncol;
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

void newpermone(double *x, int *dclass, int *n, int *xlen, int *nperm, double *zstats, double *tmat, int *rarray)

{

int i, k, l, m;
double cumsum;
int temp;


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
      m = (int)((float)l * UNIF);
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




void newpermtwo(double *x, double *y, int *n, int *xlen, int *nperm, double *zstats, double *tmat, int *rarray)

{

int i, k, l, m;
double cumsum;
int temp;
float naval = -9999;


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
      m = (int)((float)l * UNIF);
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




void psum(double *x, int *pnrow, int *pncol, double *dist)
{
int row1, row2, col1;
int nrow, ncol;
int l;
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
    

void pdiff(double *x, int *pnrow, int *pncol, double *dist)
{
int row1, row2, col1;
int nrow, ncol;
int l;
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
    

void jpres(double *x, int *pnrow, int *pncol, double *dist)
{
int row1, row2, col1;
int nrow, ncol;
int l;
double thisval, thatval;

l = 0;
nrow = *pnrow;
ncol = *pncol;

  for(col1 = 0; col1 < ncol; col1++) {
	for(row1 = 0; row1 < nrow; row1++) {
		thatval = x[row1 * ncol + col1];
		for(row2 = 0; row2 < nrow; row2++) {
			thisval = x[row2 * ncol + col1];
			if((thisval > 0) & (thatval > 0)) {
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
    

void jabs(double *x, int *pnrow, int *pncol, double *dist)
{
int row1, row2, col1;
int nrow, ncol;
int l;
double thisval, thatval;

l = 0;
nrow = *pnrow;
ncol = *pncol;

  for(col1 = 0; col1 < ncol; col1++) {
	for(row1 = 0; row1 < nrow; row1++) {
		thatval = x[row1 * ncol + col1];
		for(row2 = 0; row2 < nrow; row2++) {
			thisval = x[row2 * ncol + col1];
			if((thisval == 0) & (thatval == 0)) {
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
    

void jfirst(double *x, int *pnrow, int *pncol, double *dist)
{
int row1, row2, col1;
int nrow, ncol;
int l;
double thisval, thatval;

l = 0;
nrow = *pnrow;
ncol = *pncol;

  for(col1 = 0; col1 < ncol; col1++) {
	for(row1 = 0; row1 < nrow; row1++) {
		thatval = x[row1 * ncol + col1];
		for(row2 = 0; row2 < nrow; row2++) {
			thisval = x[row2 * ncol + col1];
			if((thisval > 0) & (thatval == 0)) {
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
    

void jsec(double *x, int *pnrow, int *pncol, double *dist)
{
int row1, row2, col1;
int nrow, ncol;
int l;
double thisval, thatval;

l = 0;
nrow = *pnrow;
ncol = *pncol;

  for(col1 = 0; col1 < ncol; col1++) {
	for(row1 = 0; row1 < nrow; row1++) {
		thatval = x[row1 * ncol + col1];
		for(row2 = 0; row2 < nrow; row2++) {
			thisval = x[row2 * ncol + col1];
			if((thisval == 0) & (thatval > 0)) {
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


void mrmperm(double *x, double *y, int *p, int *nd, int *n, int *nperm, double *r2all, double *ball, double *fall, double *tmat, int *rarray, double *XX, double *XY, double *YY, double *b)

{

int i, k, l;
int m;
int temp;
double SSE=0.0, SSTO=0.0, SSR=0.0;
double r2=0, f=0;
double btemp=0.0;
int bcount = 0;
char *transt = "T", *transn = "N";
double one = 1.0, zero = 0.0;
int onei = 1;


/* Set random seed using R function */

RANDIN;

/* Start permutation routine */
for(i = 0; i < *nperm; i++) {

/* first do the unpermuted values */

/*	F77_CALL(dgemm)(transa, transb, &ncx, &ncy, &nrx, &one,
			x, &nrx, y, &nry, &zero, z, &ncx FCONE FCONE); */

/* take crossproduct t(X) %*% Y - WORKS */
    F77_CALL(dgemm)(transt, transn, 
            p, &onei, nd, 
            &one, x, nd, y, nd, 
            &zero, XY, p FCONE FCONE);

/* take crossproduct t(Y) %*% (Y) - WORKS */
    F77_CALL(dgemm)(transt, transn, 
            &onei, &onei, nd, 
            &one, y, nd, y, nd, 
            &zero, YY, &onei FCONE FCONE);

/* calculate regression coefficients XX %*% XY - WORKS */
    F77_CALL(dgemm)(transn, transn, 
            p, &onei, p, 
            &one, XX, p, XY, p, 
            &zero, b, p FCONE FCONE);

/* calculate regression components - WORKS */
    F77_CALL(dgemm)(transt, transn, 
            &onei, &onei, p, 
            &one, b, p, XY, p, 
            &zero, &btemp, &onei FCONE FCONE);

/* SSE - WORKS */    
    SSE = YY[0] - btemp;

/* SSTO - WORKS */  
    SSTO = 0;
    for(k = 0; k < *nd; k++) {
        SSTO = SSTO + y[k];
    }
    SSTO = YY[0] - (SSTO * SSTO) / *nd;

    SSR = SSTO - SSE;

    /* calculate R2 - WORKS */
    r2 = 1 - SSE / SSTO;

    /* calculate F* - WORKS */
    f = (SSR / (*p - 1)) / (SSE / (*nd - *p));

    r2all[i] = r2;
    fall[i] = f;

    /* calculate pseudo-t for regression coefficients - WORKS*/
    /* b / sqrt(1 - R2) */
    for(k=0; k<*p; k++) {
        ball[bcount] = b[k] / sqrt(1 - r2);
        bcount++;
    }


/* permute Y */
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
      m = (int)((float)l * UNIF);
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

}

/* Reset random seed using an Splus function. */

RANDOUT;

}


void xpermute(double *x, double *y, int *nrow, int *ncol, int *xlen, int *nperm, double *zstats, double *newx, int *rarray, int *carray)

{

int i, k, l, m;
double cumsum;
int temp;
int newk, newl;


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

cumsum = 0;

/* Set up rarray. */

   for(k = 0; k < *nrow; k++) {
      rarray[k] = k;
   }


/* Set up carray. */

   for(k = 0; k < *ncol; k++) {
      carray[k] = k;
   }

/* Randomize rarray using an Splus function. */

   for(k = 0; k < (*nrow - 1); k++) {
      l = *nrow - k - 1;
      m = (long)((float)l * UNIF);
      if(m > l) m = l;
      temp = rarray[l];
      rarray[l] = rarray[m];
      rarray[m] = temp;
   }


/* Randomize carray using an Splus function. */

   for(k = 0; k < (*ncol - 1); k++) {
      l = *ncol - k - 1;
      m = (long)((float)l * UNIF);
      if(m > l) m = l;
      temp = carray[l];
      carray[l] = carray[m];
      carray[m] = temp;
   }

/* Reorder x. */

	/* loop thru the rows
	 * swapping each value with its replacement */

   for(k = 0; k < *nrow; k++) {
   }


	for(l = 0; l < *nrow; l++) {
		newl = rarray[l];
		if(newl != l) {
			for(k = 0; k < *ncol; k++) {
				newx[k*(*nrow) + l] = x[k*(*nrow) + newl];
			}
		}
	}

	/* now x has the original info and newx has swapped rows */
	/* go thru x and set x identical to newx before swapping columns */


	for(k = 0; k < *ncol; k++) {
		for(l = 0; l < *nrow; l++) {
			x[k*(*nrow) + l] = newx[k*(*nrow) + l];
		}
	}


	/* loop thru the columns
	 * swapping each value with its replacement */

   for(k = 0; k < *ncol; k++) {
   }

	for(k = 0; k < *ncol; k++) {
		newk = carray[k];
		if(newk != k) {
			for(l = 0; l < *nrow; l++) {
				newx[k*(*nrow) + l] = x[newk*(*nrow) + l];
			}
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



void xpermpart(double *hmat, double *y, double *xcor, double *ycor, int *nrow, int *ncol, int *xlen, int *nperm, double *zstats, double *newy, int *rarray, int *carray)

{

int i, k, l, m;
double cumsum;
int temp;
int newk, newl;


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

   for(k = 0; k < *nrow; k++) {
      rarray[k] = k;
   }


/* Set up carray. */

   for(k = 0; k < *ncol; k++) {
      carray[k] = k;
   }

/* Randomize rarray using an Splus function. */

   for(k = 0; k < (*nrow - 1); k++) {
      l = *nrow - k - 1;
      m = (long)((float)l * UNIF);
      if(m > l) m = l;
      temp = rarray[l];
      rarray[l] = rarray[m];
      rarray[m] = temp;
   }


/* Randomize carray using an Splus function. */

   for(k = 0; k < (*ncol - 1); k++) {
      l = *ncol - k - 1;
      m = (long)((float)l * UNIF);
      if(m > l) m = l;
      temp = carray[l];
      carray[l] = carray[m];
      carray[m] = temp;
   }



/* Reorder y. */


	/* loop thru the rows
	 * swapping each value with its replacement */

   for(k = 0; k < *nrow; k++) {
   }


	for(l = 0; l < *nrow; l++) {
		newl = rarray[l];
		if(newl != l) {
			for(k = 0; k < *ncol; k++) {
				newy[k*(*nrow) + l] = y[k*(*nrow) + newl];
			}
		}
	}

	/* now y has the original info and newy has swapped rows */
	/* go thru y and set y identical to newy before swapping columns */


	for(k = 0; k < *ncol; k++) {
		for(l = 0; l < *nrow; l++) {
			y[k*(*nrow) + l] = newy[k*(*nrow) + l];
		}
	}


	/* loop thru the columns
	 * swapping each value with its replacement */

   for(k = 0; k < *ncol; k++) {
   }

	for(k = 0; k < *ncol; k++) {
		newk = carray[k];
		if(newk != k) {
			for(l = 0; l < *nrow; l++) {
				newy[k*(*nrow) + l] = y[newk*(*nrow) + l];
			}
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



