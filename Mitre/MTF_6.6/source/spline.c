/**********************************************************************
                                NOTICE
                                ------

                This "sinemtf" software was developed by 
                        The MITRE Corporation
and was produced for the U.S. Government under Contract numbers 
J-FBI-12-128, J-FBI-86-028, J-FBI-93-039, DAAB07-99-C-C201, DAAB07-00-C-C201, 
DAAB07-01-C-C201, DAAB07-01-C-N200, DAAB07-02-C-N200, W15P7T-04-C-D001,and 
W15P7T-04-C-D199 and is subject to the Rights in Data-General clause 
FAR 52.227-l4, Alt. IV (DEC 2007). 

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
-Redistribution of source code must retain the above copyright notice, this
list of conditions, and the following disclaimer.
-Redistribution in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**********************************************************************/
/******************************************************************************/
/*                                                                            */
/* Revisions to spline.c                                                      */
/*                                                                            */
/* Original generated  mal 8/9/01                                             */
/*                                                                            */
/* V5.0 Spline modified to allow external reference to some routines.         */
/*      If spline fit has errors, revert to s-curve fit instead of piecewise  */
/*                                                      mal, 11/03            */
/*                                                                            */
/* V6.0   Fixed bug caused when patch std_deviation = 0                       */
/*                                                      mal, 06/05            */
/* v6.4.03	Removed 3 no-op lines that caused compiler warnings			   */
/*                                                      mal, 04/14            */
/******************************************************************************/

/******************************************************************************/
/* This spline code is based on freeware code from two separate sources.      */
/* Pspline.f  which can be found via                                          */
/* http://cran.r-project.org/src/contrib/PACKAGES.html#pspline                */
/*                         and                                                */
/*    The smoothSpline routine and functions it calls within the              */
/* TSD (Time Series Analysis, Signal Processing, and Dynamics) Code Library   */
/*  which can be found via http://www.qmw.ac.uk/~ugte133/tsdlibs/title.htm    */
/******************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TIFF char   /* To allow mtf.h processing without hiccup */
#include "mtf.h"
#include "xterndcl.h"

/* These constants are somewhat heuristic and may need to be 
   adjusted as more experience is gained. 
*/
#define STD_FRAC  1.0                  /* Number of standard deviations
                                          used to identify an outlier */
#define LINEAR_DF 2.00002              /* Maximum effective df that is 
                                          treated as a straight line */
#define INTERSECT_RATIO  0.3333333     /* Proportion threshold for
                                          intersections between smoothed 
                                          spline and unsmoothed spline.
                                          Fewer intersections may lead
                                          to conclusion that the smoothed
                                          spline has been smoothed too 
                                          aggressively, and the unsmoothed
                                          curve will be selected instead. */

typedef struct
{
  double x, y, a, b, c, d;
} SplineParams;
 
static double *NewVector(int nl, int nh) 
   {
    double *v = (double *)calloc((unsigned) (nh - nl + 1), sizeof(double));
    if (!v)
       MTFPRINT("Could not allocate memory in NewVector().");
    return (v - nl);
   }
 
static double **NewMatrix(int nrl, int nrh, int ncl, int nch)
   {
    int i;
    double **m;
    
    m = (double **) calloc((unsigned) (nrh - nrl + 1),sizeof(double));
    if (!m)
       MTFPRINT("Could not allocate memory in newMatrix().");
    m -= nrl;
    
    for (i = nrl; i <= nrh; i++)
      {
       m[i] = (double *) calloc((unsigned) (nch - ncl + 1),sizeof(double));
       if (!m[i])
          MTFPRINT("Could not allocate memory for matrix rows.");
       m[i] -= ncl;
      }
    return (m);
   }
 /*__________________________________________*/
 
static void FreeVector(double *v, int nl, int nh) 
   {
    free((char*) (v + nl));
   }
 
 /*__________________________________________*/
 
static void FreeMatrix(double **m, int nrl, int nrh, int ncl, int nch)
   {
    int i;
    for (i = nrh; i >= nrl; i--)
       free((char*) (m[i] + ncl));
    free((char*) (m + nrl));
   }
   
 /*__________________________________________*/
 

static SplineParams *NewSplineVector(int nl, int nh) 
{
  SplineParams *s = (SplineParams *)calloc((unsigned) (nh - nl + 1), sizeof(SplineParams));
  if (!s)
    MTFPRINT( "Could not allocate memory in NewSplineVector().");
  return (s - nl);
}

 /*__________________________________________*/
   
static void FreeSplineParams(SplineParams *S, int nl, int nh) 
   {
    free((char*) (S + nl));
   }
 
/* 
   Solve a linear equation using Choleski decomposition.
*/
static int Quincunx(int n, double *u, double *v, double *w, double *q)
{
    int j; 
    
    u[-1] = 0;
    u[0] = 0;
    v[0] = 0;
    w[-1] = 0;
    w[0] = 0;
    
    for (j = 1; j < n - 1; j++)
      {
       u[j] -=  u[j - 2] *  (w[j - 2] * w[j - 2]);
       u[j] -=  u[j - 1] * (v[j - 1] * v[j - 1]);
       if (u[j] <= 0) return (-j);
       v[j] = (v[j] - u[j - 1] * v[j - 1] * w[j - 1]) / u[j];
       w[j] /=  u[j];
      }
      
    q[0] = 0;
    q[-1] = 0;
    
    for (j = 1; j < n - 1; j++)
       q[j] -= (v[j - 1] * q[j - 1] + w[j - 2] * q[j - 2]);
       
    for (j = 1; j < n - 1; j++)
       q[j] /=  u[j];
       
    q[n] = 0;
    q[n-1] = 0;
    
    for (j = n - 2; j >= 1; j--)
       q[j] -= (v[j] * q[j + 1] + w[j] * q[j + 2]);

    return(0);
} 

/* ------------------------------------------------------------------------
   An O(n) spline smoother using cubic splines.

   This program either forces a cubic spline through the given points,
   or smooths the data point y-values to force a particular degrees of 
   freedom (df), or optimizes the generalized criterion value (gcv), 
   The resultant cubic spline is output and an indication of it's
   monotonicity is output.

   Inputs are:
   
   n        number of data points
   x        vector (strictly increasing) of length n
   w        weights for each of data points
   y        y-values to be smoothed

   method   method for computing the smoothing parameter:
            1  Use lambda=0.0, i.e. no smoothing
            2  force a specific df value
            3  optimize gcv criterion

   Outputs are:
   S        spline coefficients
   df       equivalent degrees of freedom in output spline 
                 (df=2 ==> linear)
   monotone flag indicating if output spline data points are monotone

   Return Error conditions:
      0     no error
      1     n < 5    (at least 5 unique density patches must exist to perform spline)
      5     x is not strictly increasing
      6     w contains nonpositive elements
     -j     nonpositive value of diagonal element
            at index j during factorization

*/

/* 
  Invert a band-structured matrix of order n and bandwidth 2
  that has been replaced by its rational Choleksi decomp.
  On completion the inverse overwrites the input matrix x
*/

static int  bdinvspl (int n, double **x)
{
  int i, k, l;
  int ilim, ipl, ipk;
  double sum;

  /*  check for zero diagonal entries */

  for (i=0; i<n; i++) 
    if (x[0][i] <= 0.0) 
      return (10+i);

  ilim = 2;
  x[0][n-1] = 1.0/x[0][n-1];
  for (i=n-2; i>=0; i--) {
    for (l=1; l<ilim; l++) {
      sum = 0.0;
      ipl = i + l;
      for (k=1; k<ilim; k++) {
        ipk = i + k;
        if     (k == l) 
          sum -=  x[k][ipk]*x[0][ipl];
        else if (k > l) 
          sum -=  x[k][ipk]*x[k-l][ipk];
        else
          sum -=  x[k][ipk]*x[l-k][ipl];
      }
      x[2][l-1] = sum;
    }
    sum = 1.0/x[0][i];
    for (l=1; l<ilim; l++) 
      sum = sum - x[l][i+l]*x[2][l-1];
    x[0][i] = sum;
    for (l=1; l<ilim; l++)
      x[l][i+l] = x[2][l-1];
    if (ilim < 3) ilim = ilim + 1;
  }

  /*  clear upper triangle */
  for (l=0; l<2; l++)
    x[2][l] = 0.0;

  return(0);
}

/* 
  
  An O(n) smooth spline computation with penalty based on lambda.
  Assumes band-structured matrices are pre-computed.

  yhat       array of smoothed y-values
  gcv,df     values of assorted optimization criterion
                df = 2 ==> Straight line, i.e. linear regression.
  lambda     smoothing parameter.  0.0 ==> no smoothing
                                   > 0 ==> some smoothing

  Other arrays are work space.

*/

static int splcal (int n, double *w, double *y, 
            double *gcv, double *df, double lambda, int *monotone,
            double *p, double *h, double *r, double *f,
            double *u, double *v, double *z, double *Qty, 
            double *b, double *doff, double  *u1, double *v1, double *z1,
            double **A )

{
  int i,j,ier;
  int nmnorder;
  double xn, trace,sse,sum,y1,y2;

  nmnorder = n - 2;

/*  Generate bands of band-structured matrix A = M + LAMBDA * Qt Wd Q
    u1 = diagonal, v1 = first off-diagonal, z1 = secondd off-diagonal.
    b, the returned result after processing starts as Qty */

  for (i=1;i<n-1;i++) {
    u1[i] = p[i] + lambda*u[i];
    v1[i] = h[i]/3. + lambda*v[i];
    z1[i] = lambda*z[i];
    b[i] = Qty[i];
  }

  ier = Quincunx(n,u1,v1,z1,b);
  if (ier != 0) return (ier);
 

  for (i=0; i<nmnorder; i++) {
    A[0][i] = u1[i+1];
    A[1][i] = v1[i];
    A[2][i] = z1[i-1];
  }

  /*  Compute band-structured portion of inverse of matrix A
      Computed in-place. */
  ier = bdinvspl (nmnorder, A);
  if (ier != 0) return (ier);

  /* Compute the  doff = (y-d)/(w*lambda) */
  doff[0] = r[0] * b[1];
  *monotone = 1;
  for (j = 1; j < n; j++) {
    doff[j] = r[j - 1] * b[j - 1] + f[j] * b[j] + r[j] * b[j + 1];
    y2 = y[j] - lambda * w[j]*doff[j];
    y1 = y[j-1] - lambda * w[j-1]*doff[j-1];
    if( y2-y1==0.0 || (y2-y1)/(y[n-1]-y[0]) < 0.0 ) {
      *monotone = 0;
    }
  }
  
  /*  compute various evaluation criteria */

  xn    = n;

  trace = (r[0]*r[0]) * A[0][0] * w[0];

  trace += ((r[1]*r[1]) * A[0][1]
         + (f[1]*f[1]) * A[0][0]
     + 2.0*(r[1]*f[1]) * A[1][1] ) * w[1];

  sse   = doff[0]*doff[0] + doff[1]*doff[1];

  for (i=2; i<n-2; i++) {
    sum =   (r[i]  *r[i]  ) * A[0][i]
          + (f[i]  *f[i]  ) * A[0][i-1]
          + (r[i-1]*r[i-1]) * A[0][i-2];

    sum   += 2.0*r[i]*(f[i]*A[1][i]+r[i-1]*A[2][i]);
    sum   += 2.0*f[i]*r[i-1]*A[1][i-1];
    sum   *= w[i];
    trace += sum;
    sse   += doff[i]*doff[i];
  }

  trace += ((f[n-2]*f[n-2]) * A[0][n-3]
        + (r[n-3]*r[n-3]) * A[0][n-4]
    + 2.0*(f[n-2]*r[n-3]) * A[1][n-3] ) * w[n-2];

  trace +=  (r[n-2]*r[n-2]) * A[0][n-3]   * w[n-1];
  
  sse   += doff[n-2]*doff[n-2] + doff[n-1]*doff[n-1];

  *gcv   = (sse/xn)/((trace/xn)*(trace/xn));
  trace *= lambda;
  *df    = xn - trace;
  
  return (ier);
  
}     

/* 
   Golden ratio search routine, which uses two techniques for
   optimization.  Attempts to ensure that final output
   is a monotone function.
*/

static int fmm (int n, double *wvec, double *yvec, 
         double *gcv, double *df, double *lambda, int method, 
         double tol, int *monotone,
         double *pp, double *hp, double *rp, double *fp, 
         double *up, double *vp, double *zp, double *Qty, 
         double *bp, double *doff, 
         double *u1, double *v1, double *z1, double **A )

{
  int i;
  double targdf;
  int nmnorder;
  int mmx;
  double mx, fmx;
  double ratio,t1,t2,eps,tol1,a,b;
  double x,w,v,e,fu=0.0,fx=0.0,fv,fw,xm, tol2;
  double d=0.0,r,q,p,u;

  int ier;
  int XUP = 3;
  double XDN = 1.0E-10;

  targdf   = *df;
  nmnorder = n - 2;
  t1 = 0.0;
  t2 = 0.0;
  for (i=0; i<nmnorder; i++) {
    t1 = t1 + pp[i+1]/2;
    t2 = t2 + up[i+1];
  }
  ratio = t1/t2;

  eps = 1.0;
  tol1 = 2.0;
  while (tol1 > 1.0) {
    eps = eps/2.0;
    tol1 = 1.0 + eps;
  }
  eps = sqrt(eps);

  a = XDN;
  b = XUP;

  if (*lambda <= 0)
    v = .75;
  else
    v = (2.0 + log(*lambda/ratio)/2.772589)/6.0;
  w = v;
  x = v;
  e = 0.0;
  *lambda = ratio*exp(2.772589*(6.0*x - 2.0));


  ier = splcal (n, wvec, yvec, 
                gcv, df, *lambda, monotone,
                pp, hp, rp, fp, up, vp, zp, Qty, bp, doff,
                u1, v1, z1, A);
  if (ier != 0) return (ier);
      
  if (method == 3) fx = *gcv;
  if (method == 2) fx = (targdf-*df)*(targdf-*df);

  /* fprintf(stdout, "Lambda = %g   Df=%g  GCV = %.2e   Mono = %d  \n", 
          *lambda,  *df, *gcv, *monotone);  */

  fv = fx;
  fw = fx;
  fmx = fx;
  mx = x;
  mmx = *monotone;

  /*  --------------------  start main loop ----------------------- */

  xm   = 0.5*(a + b);
  tol1 = eps*fabs(x) + tol/3.0;
  tol2 = 2.0*tol1;

  while (fabs(x - xm) > (tol2 - 0.5*(b - a))) {
      
    /* is golden-section necessary? */

    if (fabs(e) <= tol1) {
      if (x >= xm) e = a - x;
      else e = b - x;
      d = 0.382*e;
    }

    else {      

      /*  fit parabola */

      r = (x - w)*(fx - fv);
      q = (x - v)*(fx - fw);
      p = (x - v)*q - (x - w)*r;
      q = 2.0*(q - r);
      if (q > 0.0) p = -p;
      q =  fabs(q);
      r = e;
      e = d;

      /*  is parabola acceptable? */

      if ( (fabs(p) >= fabs(0.5*q*r)) ||
            (p <= q*(a - x))  ||
            (p >= q*(b - x))  ) {
        if (x >= xm) e = a - x;
        else e = b - x;
        d = 0.382*e;
      }
      else {
        
        /*  a parabolic interpolation step */

        d = p/q;
        u = x + d;

        /*  f must not be evaluated too close to a or b  */

        if ((u - a) < tol2 || (b - u) < tol2) {
          if (xm - x >= 0.0) 
            d =  tol1; 
          else
            d = -tol1;
        }
      }
    }
    
    /*  f must not be evaluated too close to x */

    if (fabs(d) >= tol1) 
      u = x + d;
    else {
      if (d >= 0.0) 
        u = x + tol1; 
      else
        u = x - tol1; 
    }
    *lambda = ratio*exp(2.772589*(6.0*u - 2.0));

    ier = splcal (n, wvec, yvec, 
                  gcv, df, *lambda, monotone,
                  pp, hp, rp, fp, up, vp, zp, Qty, bp, doff,
                  u1, v1, z1, A);

    if (ier != 0) return (ier);

    if (method == 3) fu = *gcv;
    if (method == 2) fu = (targdf-*df)*(targdf-*df);

    /* fprintf(stdout, "Lambda = %g   Df=%g  GCV = %.2e  Mono = %d \n", 
       *lambda,  *df, *gcv, *monotone);  */

    /*  update  a, b, v, w, and x  */

    if (fu > fx) {
      if (u  <  x)  a = u;
      if (u  >= x)  b = u;
      if (fu <= fw  ||  w == x) {
        v  = w;
        fv = fw;
        w  = u;
        fw = fu;
      }
      else if (fu <= fv || v == x || v == w) {
        v  = u;
        fv = fu;
      }
    } else {
      if (u  >= x)  a = x;
      if (u  <  x)  b = x;
      v  = w;
      fv = fw;
      w  = x;
      fw = fx;
      x  = u;
      fx = fu;
    }
    if (*monotone == 1) {
      if (fu <= fmx || mmx == 0 ) {
        mmx = *monotone;
        mx  = u;
        fmx = fu;
      }
    }
    
    xm   = 0.5*(a + b);
    tol1 = eps*fabs(x) + tol/3.0;
    tol2 = 2.0*tol1;
  }

  if (*monotone == 0 && mmx == 1) {
    *lambda = ratio*exp(2.772589*(6.0*mx - 2.0));

    ier = splcal (n, wvec, yvec, 
                  gcv, df, *lambda, monotone,
                  pp, hp, rp, fp, up, vp, zp, Qty, bp, doff,
                  u1, v1, z1, A);

    if (ier != 0) return (ier);

    /* fprintf(stdout, "Reprise  Lambda = %g   Df=%g  GCV = %.2e  Mono = %d \n",
            *lambda,  *df, *gcv, *monotone); */
  }

  return (0);

} 

/*
  Initialize working arrays for spline smoothing computations. 
  Several band diagonal matrices are computed.
*/ 
static void initialize_fixed_spline_data
                            (int n, double *x, double *w, double *y, 
                             double **pp, double **hp, double **rp,
                             double **fp, double **up, double **vp,
                             double **zp, double **Qtyp,
                             double **bp, double **doffp)
{
  double *h, *r, *f, *p, *u, *v, *z, *Qty;
  int i, nm1;

  nm1 = n-1;

  h  =  NewVector(0, nm1);
  r  =  NewVector(0, nm1);
  f  =  NewVector(0, nm1);
  p  =  NewVector(0, nm1);
  u  =  NewVector(0, nm1);
  v  =  NewVector(0, nm1);
  z  =  NewVector(0, nm1);
  Qty=  NewVector(0, n);

  h[0] = (x[1] - x[0]);
  r[0] = 1. / h[0];
  f[0] = 0;
  p[0] = 0;
    
  for (i = 1; i < nm1 ; i++) {
    h[i] = (x[i+1] - x[i]);
    r[i] = 1. / h[i];
    f[i] = -(r[i - 1] + r[i]);
    p[i] = 2. * (x[i+1] - x[i-1]) / 3.;
    Qty[i] =  (y[i+1] - y[i]) / h[i] 
                 - (y[i] - y[i-1]) / h[i - 1];
  }
  r[nm1] = 0;
  f[nm1] = 0;

  for (i = 1; i < nm1; i++) {
    u[i] = r[i - 1] * r[i - 1] * w[i - 1];
    u[i] += f[i] * (f[i] * w[i]);
    u[i] += r[i] * (r[i] * w[i + 1]);
    v[i] = f[i] * r[i] * w[i]; 
    v[i] += r[i] * f[i + 1] * w[i + 1];
    z[i] =  r[i] * r[i + 1] * w[i + 1];
  }
  
  *bp = NewVector(-1, n);
  *doffp = NewVector(0, nm1);

  *hp = h;
  *rp = r;
  *fp = f;
  *pp = p;
  *up = u;
  *vp = v;
  *zp = z;
  *Qtyp = Qty;
}

static int find_smooth_spline (int n, double *x, double *w, double *y,
                        SplineParams *S, double *df, int *monotone, 
                        int method)
{

  int nm1, i;
  double xi;
  double *h, *r, *f, *p, *u, *v, *z, *Qty, *b, *doff;
  double *u1, *v1, *z1;
  double **A;
  double gcv;
  

  double lambda = 0.0;
  double tol = 1.0E-3;
  int ier = 0;

  /*  check arguments  */
  if (n < 5) 
    return (1);

  /*  Check for x strictly increasing */

  for (i=0; i<n; i++) {
    if (w[i] <= 0.0)
      ier = 6;
    xi = x[i];
    if (i >= 1 && xi <= x[i-1]) 
      ier = 5;
  }
  
  if (ier != 0) return(ier);

  initialize_fixed_spline_data(n, x, w, y, &p, &h, &r, &f, &u, &v,
                               &z, &Qty, &b, &doff);
  
  u1 = NewVector(-1, n-1);
  v1 = NewVector(-1, n-1);
  z1 = NewVector(-1, n-1);
  A = NewMatrix(0,2,0,n-3);

  if (method == 1) {
    ier = splcal (n, w, y, &gcv, df, lambda, monotone,
                  p, h, r, f, u, v, z, Qty, b, doff, u1, v1, z1, A);

    /* fprintf(stdout, "Lambda = %g   Df=%g  GCV = %.2e  Mono = %d\n", 
            lambda,  *df, gcv, *monotone);  */
  }
  else {
    method = 3;
    ier = fmm (n, w, y, &gcv, df, &lambda, method, tol, monotone,
               p, h, r, f, u, v, z, Qty, b, doff, u1, v1, z1, A);
    if (ier != 0) return (ier);
    
    if (method > 2 && *df > n) {
      *df = n;
      ier = fmm (n, w, y, &gcv, df, &lambda, 2, tol, monotone,
                 p, h, r, f, u, v, z, Qty, b, doff, u1, v1, z1, A);
    }
  }

  S[0].d = y[0] - lambda * doff[0] * w[0];
  S[1].d = y[1] - lambda * doff[1] * w[1];
  S[0].a = b[1] / (3 * h[0]);
  S[0].b = b[0];
  S[0].c = (S[1].d - S[0].d) / h[0] - b[1] * h[0] / 3;
  S[0].x = x[0];
  S[0].y = y[0];
  
  nm1 = n-1;
  for (i = 1; i < nm1 ; i++)
    {
      S[i].a   = (b[i + 1] - b[i]) / (3 * h[i]);
      S[i].b   =  b[i];
      S[i].c   = (b[i] + b[i - 1]) * h[i - 1] + S[i - 1].c;
      S[i+1].d = y[i+1] - lambda * doff[i+1] * w[i+1];
      S[i].x   = x[i];
      S[i].y   = y[i];
    }
  S[nm1].c   = (b[nm1] + b[n-2]) * h[n-2] + S[n-2].c;
  S[nm1].b   = 0.0;
  S[nm1].a   = 0.0;
  S[nm1].x   = x[nm1];
  S[nm1].y   = y[nm1];

  FreeVector(h, 0, nm1);
  FreeVector(r, 0, nm1);
  FreeVector(f, 0, nm1);
  FreeVector(p, 0, nm1);
  FreeVector(u, 0, nm1);
  FreeVector(v, 0, nm1);
  FreeVector(z, 0, nm1);
  FreeVector(Qty, 0, n);
  FreeVector(b, -1, n);
  FreeVector(doff, 0, nm1);  
  FreeVector(u1,-1, n-1);
  FreeVector(v1,-1, n-1);
  FreeVector(z1,-1, n-1);
  FreeMatrix(A,0,2,0,n-3);


  return (ier);
}     

/* 
   Compare original data points with actual smoothed data points,
   and determine if any one point has been smoothed an abnormally 
   large amount.  Points in the original data set that move an
   undue amount are considered outliers, and removed from the 
   original data set.
*/

static int find_remove_outliers
                    (SplineParams *S, double *x, double *w, double *y, 
                     double *ind,
                     int *n) 
{
  
  int m1,m2,m3, not_finished;
  int i, num_points, index;
  double val1,val2,val3;
  double sumy, sumysq, meany, std, diff;

  num_points = *n;

  /* Find 3 highest offsets of original points from smoothed curve */
  m1 = m2 = m3 = 0;
  val1 = val2 = val3 = fabs(y[0] - S[0].d)/w[0];
  sumy = sumysq = 0;
  for(i=0;i<num_points;i++) { 
    diff =  (y[i] - S[i].d)/w[i];
    sumy += diff;
    sumysq += diff*diff;
    if ( fabs(diff) > val1 ) {
      m3 = m2;
      val3 = val2;
      m2 = m1;
      val2 = val1;
      m1 = i;
      val1 = fabs(diff);
    }
    else if ( fabs(diff) >= val2  || m1 == m2){
      m3 = m2;
      val3 = val2;
      m2 = i;
      val2 = fabs(diff);
    }
    else if ( fabs(diff) >= val3 || m3 == m1 || m2 == m3){
      m3 = i;
      val3 = fabs(diff);
    }
  }
  meany = sumy/(double)num_points;
  std = sumysq/(double)num_points - meany*meany;
  std = sqrt(std);

  /*  Print out for debugging 
      fprintf(stdout, "Max Diff at %d ,  meandiff = %f,  stddiff = %f\n", m1,
      meany, std);
      fprintf(stdout, "Max diffs: %f, %f, %f   Remove Last? %d\n",
      val3, val2, val1, ((val1-val3)>std &&
      (y[m1]-y[m1+1])*(y[m1-1]-y[m1])<0));
  */
  
  not_finished = 0;
  /* Check and see if an actual point is too far off the smoothed curve, 
     if so remove it */
  if ( ((val1-val3)> STD_FRAC*std && (y[m1]-y[m1+1])*(y[m1-1]-y[m1])<0) && (num_points > 4))
  {
    MTFPRINT2("Removed patch %s from smoothed spline curve fit (appears to be a flyer).\n", 
            g_test_pattern[(int)ind[m1]].new_patch_label);

    not_finished = 1;
    index = (int)(ind[m1]);
    for (i=m1; i<num_points-1; i++) {
      x[i] = x[i+1];
      y[i] = y[i+1];
      w[i] = w[i+1]; 
      ind[i] = ind[i+1];
    }
    ind[num_points-1] = index;

    (*n)--;
  }
  return(not_finished);

}

/* This function determines if the smoothed or unsmoothed spline
   is preferable.  The unsmoothed data must pass 2 tests to
   be better.  1) it is monotone,  2) it crosses the smoothed
   spline infrequently.
*/

static int choose_curve(SplineParams *S, double *y, int num_points, 
                        int *monotone) 
{
  
  int i, cnt;
  int diff, lastdiff;
  double ydiff, endtoend, error;
  
  
  *monotone = 1;
  cnt = 0;
  error = y[0] - S[0].d;
  if (error > 0) lastdiff = 1;
  else if (error<0) lastdiff = -1;
  else lastdiff = -2;
  endtoend =y[num_points-1]-y[0];
  
  for(i=1; i<num_points; i++) {
    ydiff = y[i]-y[i-1];
    if( ydiff==0.0 || ydiff/endtoend < 0 )
      *monotone = 0;
    error = y[i] - S[i].d;
    if (error > 0) diff = 1;
    else if (error<0) diff = -1;
    else diff = 0;

    if( lastdiff != -2 || diff != 0 ){
      if( diff != lastdiff && diff != 0 ) {
	cnt++;
	lastdiff = diff;
      }
    }
  }

  /* Debug smooth vs unsmoothed 
  printf( "Monotone? %d  #Intersections %d\n", *monotone, cnt);
  */

  if ( *monotone == 1 && cnt < num_points*INTERSECT_RATIO ) 
    return (0); /* Unsmooth */
  return(1);  /* Smooth */
}

static void remove_squiggly_bits(SplineParams *S, int num_points) 
{
  int i;
  double root, low_turn=0.0, high_turn=0.0;
  
  for(i=0;i<num_points-1;i++) {
    root = S[i].b*S[i].b - 3*S[i].a*S[i].c;
    if ( S[i].a == 0 ) {
      low_turn = -S[i].c/(2*S[i].b);
      high_turn = -S[i].c/(2*S[i].b);
    }
    else {
      if (root >= 0.0) {
        root = sqrt(root);
        if ( S[i].a > 0 ) {
          low_turn = (-S[i].b - root)/(3*S[i].a);
          high_turn = (-S[i].b + root)/(3*S[i].a);
        }
        else if ( S[i].a < 0 ) {
          high_turn = (-S[i].b - root)/(3*S[i].a);
          low_turn = (-S[i].b + root)/(3*S[i].a);
        }
      }
      else {
        low_turn = high_turn = HUGE_VAL;
      }
    }

    if( (S[i+1].x-S[i].x) <= low_turn || 0 >= high_turn ||
        (low_turn <= 0 && S[i+1].x-S[i].x<=high_turn) ) continue;
    else {
      /*      spline fit has doubled back on itself at intermediate values so
	      substitute a local linear fit  */
      S[i].a = 0;
      S[i].b = 0;
      S[i].c = (S[i+1].d-S[i].d)/(S[i+1].x-S[i].x);
      /*
      MTFPRINT3("Spline segment between %7.3f and %7.3f ImgGrays turned into straight line\n",
      		 sort_reflectance[i].grey, sort_reflectance[i+1].grey );
      */
    }
  }
}

static int compare_reflectance_pairs(const void *first, const void *second)
{

  if (((struct spline_reflectance_pair *)first)->reflectance < 
      ((struct spline_reflectance_pair *)second)->reflectance)
    return -1;
  
  if (((struct spline_reflectance_pair *)first)->reflectance > 
      ((struct spline_reflectance_pair *)second)->reflectance)
    return 1;

  return 0;
}

/* Sort reflectance-grey-origindex data for density patches. */
void set_up_sort_reflectance_pairs(void)
{
  unsigned int i;
  short j;

    sort_reflectance
    = (struct spline_reflectance_pair *)calloc(number_unique_density_patches, 
                                        sizeof(struct spline_reflectance_pair));
  
  if (sort_reflectance == NULL)
    {
      fprintf(stderr, "Error: Could not calloc sufficient space for reflectance sorted data\n");
      perror("set_up_sort_reflectance_pairs");
    }

  /* Initialize spline_reflectance_pairs */
  j = 0;
  for (i = 0;  i < g_number_of_patches;  i++)
    if (g_test_pattern[i].patch_type == DENSITY_PATCH)
      {
        sort_reflectance[j].grey = g_test_pattern[i].final_grey;
        sort_reflectance[j].reflectance = g_test_pattern[i].final_reflectance;
        sort_reflectance[j].index = i;
        j++;
      }

  /* Sort spline_reflectance_pairs */
  qsort(sort_reflectance, number_unique_density_patches,
        sizeof(struct spline_reflectance_pair), compare_reflectance_pairs);
}


int smooth_spline(int *monotone, double *chisq)
{
  int i,j, num_points, ier;
  double *x, *y, *w, *ind;
  SplineParams *S;
  int not_finished;
  int use_smooth;
  int smooth_linear;
  double df, temp, val, xdiff;
  double min_std;

  set_up_sort_reflectance_pairs();
  
  num_points = number_unique_density_patches;

  x = NewVector(0,num_points-1);
  y = NewVector(0,num_points-1);
  w = NewVector(0,num_points-1);
  ind = NewVector(0,num_points-1);

  for(i=0;i<num_points;i++) {
    x[i] = sort_reflectance[i].reflectance;
    y[i] = sort_reflectance[i].grey;
    w[i] = num_points-1;  /* Any constant value is the same */
    ind[i] = sort_reflectance[i].index;
  }

  /* Having the very first point be non-monotonic for printer curves,
     tends to confuse the spline fit.  To avoid this problem, either
     drop the point, or de-weight it (as is shown here).  */
  /* Not used for now (V5.0), but might be useful in the future.
     if( (sort_reflectance[0].grey - sort_reflectance[1].grey)/
         (sort_reflectance[0].grey - sort_reflectance[num_points-1].grey) < 0 )
       w[0] = 1./(num_points-1);
  */

  S = NewSplineVector(0,num_points-1);
  
  /* Find smooth spline while checking for outliers in original data.*/
  not_finished = 1;
  while (not_finished){
    
    ier = find_smooth_spline (num_points, x, w, y, S, 
                              &df, monotone, 3);

    if (ier == 5) 
      MTFPRINT("ERROR: Multiple patches at same reflectance\n" )
    else if (ier == 1) 
      MTFPRINT ("Cannot do spline fit, requires at least 5 density patches.\n")
    else if (ier) 
      MTFPRINT2("ERROR (%d): During computation of spline smoothing\n", ier)
    if (ier > 1 || ier < 0) mtfexit(-1);

    not_finished = find_remove_outliers(S, x, w, y, ind, &num_points);

  }  /* End of while loop */

  /* See if smoothed curve is actually linear. */
  smooth_linear = (df < LINEAR_DF && 
                   num_points == (int)number_unique_density_patches)? 1:0;

  if (smooth_linear) {
    for(i=0;i<num_points;i++) {
      S[i].a = 0;
      S[i].b = 0;
      S[i].c = g_slope;
      S[i].d = x[i]*g_slope + g_intercept;
    }
  }

  /* See if smooth or unsmooth is better */
  use_smooth = choose_curve(S, y, num_points, monotone);

  /* Set Sspline values, and determine minumum non-zero std_deviation */
  min_std = HUGE_VAL;
  for(i=0;i<num_points;i++) {
  	if (g_test_pattern[i].std_deviation < min_std && g_test_pattern[i].std_deviation > 0.0) 
  		min_std = g_test_pattern[i].std_deviation;
    g_test_pattern[(int)ind[i]].Sspline = S[i].d;
  }
  min_std /= 10.0; 

  /* Figure out smoothed spline points for any points that were
     removed during spline calculation */
  for( ;i<(int)number_unique_density_patches; i++) {
    temp =  g_test_pattern[(int)ind[i]].final_reflectance;
    j = num_points-1;
    while (j>0 && temp < S[j].x) j--;

    xdiff = temp-S[j].x;
    if(xdiff<0)
      val = S[0].d + S[0].c*(xdiff);
    else
      val = S[j].d + S[j].c*xdiff + S[j].b*xdiff*xdiff +
        S[j].a*xdiff*xdiff*xdiff;
    
    g_test_pattern[(int)ind[i]].Sspline = val;
  }
  
  FreeVector(x,0,num_points-1);
  FreeVector(y,0,num_points-1);
  FreeVector(w,0,num_points-1);
  FreeVector(ind,0,num_points-1);
  FreeSplineParams(S,0,num_points-1);

  /* Set automatic curve type */
  *chisq = 0.0;
  if (smooth_linear && use_smooth) {
    for(i=0;i<(int)number_unique_density_patches; i++) {
      if( g_test_pattern[i].std_deviation != 0.0)
      	xdiff = (g_test_pattern[i].final_reflectance*g_slope + g_intercept 
	     	 - g_test_pattern[i].final_grey)/g_test_pattern[i].std_deviation;
	  else
	  	xdiff = (g_test_pattern[i].final_reflectance*g_slope + g_intercept 
	     	 - g_test_pattern[i].final_grey)/min_std;
      *chisq += (xdiff*xdiff);
    }
    *chisq /= (double)number_unique_density_patches;
    return (6);   /* Linear */
  }
  else if (use_smooth) {
    for(i=0;i<(int)number_unique_density_patches; i++) {
      if( g_test_pattern[i].std_deviation != 0.0)
        xdiff = (g_test_pattern[i].Sspline - g_test_pattern[i].final_grey)/
			g_test_pattern[i].std_deviation;
	  else 
	    xdiff = (g_test_pattern[i].Sspline - g_test_pattern[i].final_grey)/
	    		min_std;
      *chisq += (xdiff*xdiff);
    }
    *chisq /= (double)number_unique_density_patches;
    return (8);   /* Smooth spline */
  }
  else if (*monotone) 
    return (7);  /* Unsmooth spline */
  else 
    return (9);  /* S-curve */

  
}

/* 
   Invert the spline so that reflectance can be computed from grey value.
   This requires that new spline parameters be computed with the chosen set
   of data points (either smooth or unsmooth grey values).  This spline fit
   is done without any additional smoothing.
   
   sort_reflectance[].grey must be correctly set here so that 
   compute_spline_grey_to_reflectance() can be evaluated correctly later.
*/

void invert_spline (void)
{
  int num_points, i;
  double val;
  double *x, *y, *w;
  double df;
  int monotone, ier;
  SplineParams *S;

  num_points = number_unique_density_patches;

  /* Invert the computation */
  x = NewVector(0,num_points-1);
  y = NewVector(0,num_points-1);
  w = NewVector(0,num_points-1);
  
  for(i=0;i<num_points;i++) {
    x[i] = sort_reflectance[i].reflectance;
    if( g_followcurve == 8)
      sort_reflectance[i].grey = 
        g_test_pattern[sort_reflectance[i].index].Sspline;
    y[i] = sort_reflectance[i].grey;
    w[i] = 1;
  }

  /* If yhat is decreasing then reverse data order */
  if(y[0] > y[num_points-1]) {
    for(i=0; i<num_points/2; i++) {
      val = x[i];
      x[i] = x[num_points-1-i];
      x[num_points-1-i] = val;
      val = y[i];
      y[i] = y[num_points-1-i];
      y[num_points-1-i] = val;
      val = sort_reflectance[i].grey;
      sort_reflectance[i].grey = sort_reflectance[num_points-1-i].grey;
      sort_reflectance[num_points-1-i].grey = val;
      val = sort_reflectance[i].reflectance;
      sort_reflectance[i].reflectance = sort_reflectance[num_points-1-i].reflectance;
      sort_reflectance[num_points-1-i].reflectance = val;
      val = sort_reflectance[i].index;
      sort_reflectance[i].index = sort_reflectance[num_points-1-i].index;
      sort_reflectance[num_points-1-i].index = (int)val;
    }
  }

  if( g_followcurve == 7 || g_followcurve == 8 ) {
    S = NewSplineVector(0,num_points-1);
    /* Force fit (y,x) plot to the points found above */
    ier = find_smooth_spline (num_points, y, w, x, S, &df, &monotone, 1);
    if (ier == 5) 
      MTFPRINT("\nERROR2: Spline Graylevels are not strictly increasing\n" )
    else if (ier) 
      MTFPRINT("\nERROR2: Invert_Spline computation error\n")
    if (ier) { 
      g_followcurve = 9;
      MTFPRINT("Switching to S-curve Fit\n")
    }
    else {
      remove_squiggly_bits(S, num_points);
      for(i=0;i<num_points;i++) {
        sort_reflectance[i].slope = S[i].c ;
        sort_reflectance[i].quadratic = S[i].b ;
        sort_reflectance[i].cubic = S[i].a ;
      }
      FreeSplineParams(S,0,num_points-1);
    }
  }
  FreeVector(x,0,num_points-1);
  FreeVector(y,0,num_points-1);
  FreeVector(w,0,num_points-1);

}


/* 
   Given a set of spline coefficients or a linear equation for
   the reflectance curve, compute the reflectance at the given 
   grey_level.  Reflectance values are clipped to be within the 
   range [0,1].
*/

double compute_spline_grey_to_reflectance (double grey_level)
{
  int i;
  double gdiff, val;

  i = number_unique_density_patches-1;
  while(i>0 && grey_level < sort_reflectance[i].grey) i--;

  gdiff = grey_level - sort_reflectance[i].grey;
    
  if(gdiff<0)
    val = sort_reflectance[0].reflectance + 
      sort_reflectance[0].slope * gdiff;
  else
    val = sort_reflectance[i].reflectance + 
      sort_reflectance[i].slope * gdiff +
      sort_reflectance[i].quadratic * gdiff * gdiff +
      sort_reflectance[i].cubic * gdiff * gdiff * gdiff;
  
  if(val > 1.0) val = 1.0;
  else if(val < 0.0) val= 0.0;
  
  return (val);
}


