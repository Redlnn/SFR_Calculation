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
/*********************************************************************/
/*                                                                   */
/* Revisions to scurve.c                                             */
/*                                                                   */
/* Original generated  mal 11/03   mtf V5.0                          */
/*                                                                   */
/* V5.5     Cleaned up compiler warnings.                            */
/*                                             mal, 01/04            */
/* V5.8.1   Fixed bug in scurve expanded MTF printout.               */
/*                              (leading 1.0s replaced by 0.0)       */
/*                                             mal, 03/05            */
/*                                                                   */
/* V6.0     Fixed bug caused when patch std_deviation = 0           */
/*                                             mal, 06/05            */
/*                                                                   */
/*********************************************************************/

/*********************************************************************/
/* This s-curve fitting code is based upon.                          */
/* Meyer, Yung, and Ausubel, "A Primer to Logistic Growth and        */
/* Substitution: The Mathematics of the Loglet Lab Software",        */
/* Technogical Forecasting and Social Change 61(3), pp 247-271, 1999.*/
/* http://phe.rockefeller.edu/LogletLab/logletlab.pdf                */
/* and                                                               */
/* Lampton, "Damping-undamping strategies for the Levenberg-Marquardt*/
/* nonlinear least-squares method", Computers in Physics, Vol 11,    */
/* no 1, pp 110-115, 1997.                                           */
/*********************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TIFF char
#include "mtf.h"
#include "xterndcl.h"

#define NPARAM 4

double g_a, g_k, g_m, g_c;
 
/* Compute the logistic function and partial deriviatives at 
 *  position x. 
 *  
 *  Input:
 *      x = evaluation point
 *      p = vector of logistic parameters
 *         p[0] = kappa;
 *	   p[1] = c;
 *	   p[2] = m;
 *	   p[3] = a;
 *
 * Output:
 *     y = logistic value at x
 *     dydp = partial of logistic for each parameter p.
 *      
 *
 * The logistic formula used is
 * 
 *    y = a + kappa / [1 + exp(-mx-b)]
 */
void logistic(double x, double *p, double *y, double *dydp)
{
  double G, H, H2;

  G = exp(-(p[2]*x+p[1]));
  H = 1./(1.+G);
  *y = p[0]*H;
  H2 = (*y)*G*H;
  dydp[0] = H;
  dydp[1] = H2;
  dydp[2] = x*H2;
  dydp[3] = 1.;
  *y +=  p[3];
  
}

/* Compute the curvature matrix, gradient vector, and goodness of 
 *   fit at a particular model parameter setting.
 *   Inputs
 *      x = fixed tgt reflectances
 *      y = measured grey levels
 *      std = std deviation of measured grey levels
 *      n = number of points
 *      p = model parameter setting
 *   Intermediate space
 *      dydp = Jacobian 
 *   Outputs
 *      alpha = curvature matrix (symmetric, only filled in upper left)
 *      beta = gradient vector
 *      chisq = goodness of fit
 */
void partials(double *x, double *y, double *std, int n,
	      double *p, double *dydp,   /* Intermediate space */
	      double **alpha, double *beta, double *chisq)
{

  int i,j,k;
  double yhat, var, res, wt;

  for(j=0; j<NPARAM; j++) {
    for(k=0;k<NPARAM; k++) 
      alpha[j][k] = 0.;
    beta[j] = 0.;
  }
  
  *chisq = 0.0;
  for(i=0; i<n; i++) {
    var = 1./(std[i]*std[i]);
    logistic(x[i],p,&yhat,dydp);
    res = y[i] - yhat;
    for(j=0;j<NPARAM;j++) {
      wt = dydp[j]*var;
      beta[j] += res*wt;
      for(k=j;k<NPARAM;k++)
	alpha[j][k] += wt*dydp[k];
    }
    *chisq += res*res*var;
  }
}

/* Levenberg-Marquardt Non-Linear Regression to a logistic S-curve
 * model.
 *      f(x;a,kappa,m,c) = a + kappa / [1 + exp (-mx-c)]
 *
 *  where more typical logistic parameters dt and tm are
 *      dt = log(81)/m;
 *      tm = -c/m;
 *
 *
 *  Inputs:
 *     x = sample target reflectances
 *     y = measured image grey levels
 *     std = std deviation at each measured grey level
 *     n = number of sample points
 *  Input/Output logistic curve parameters:
 *     see equation above
 *     The input values of there parameters are used as a starting
 *     point. After iterating to convergence, they are reset to
 *     their output value.
 *     
 */
double levmarq(double *x, double *y, double *std, int n, 
	   double *a, double *b, double *kappa, double *m, double *c)
{
  int j,k,first,*ps, noerr;
  double chisq, newchi, ochisq, lambda;
  double *p, **alpha, *beta, *dydp, *scales;
  double *tryp, **mat, *tryb;


  /* Allocate all the temporary work space needed for iterative 
     computations */
  ps = (int *)malloc(NPARAM*sizeof(int));
  scales = (double *)malloc(NPARAM*sizeof(double));
  p = (double *)malloc(NPARAM*sizeof(double));
  tryp = (double *)malloc(NPARAM*sizeof(double));
  tryb = (double *)malloc(NPARAM*sizeof(double));
  dydp = (double *)malloc(NPARAM*sizeof(double));
  beta = (double *)malloc(NPARAM*sizeof(double));
  mat = (double **)malloc(NPARAM*sizeof(double *));
  mat[0] = (double *)malloc(NPARAM*NPARAM*sizeof(double));
  alpha = (double **)malloc(NPARAM*sizeof(double *));
  alpha[0] = (double *)malloc(NPARAM*NPARAM*sizeof(double));
  for(j=1;j<NPARAM;j++) {
    alpha[j] = alpha[j-1] + NPARAM;
    mat[j] = mat[j-1] + NPARAM;
  }
  
  /* Set the starting parameters */
  p[0] = *kappa;
  p[2] = *m;
  p[1] = *c;
  p[3] = *a;
  
  lambda = 0.01; 

  /* Find goodness of fit, curvature and gradient at initial 
     parameters */
  partials(x, y, std, n, p, dydp, alpha, beta, &chisq);

  /*
  MTFPRINT8("\nChi = %g  A=%g B=%g K=%g Tm=%g dT=%g   lambda=%g\n",
	  sqrt(chisq), p[3], p[3],p[0],-p[1]/p[2],log(81)/p[2], lambda);
  */

  /* Search for model parameters which improve the goodness of fit. */
  /* Continue search until either:
       . lambda becomes extremely large or extremely small 
       . a zero-error fit found
       . last improvement was small, and it is the 2nd small improvement
  */
  ochisq = chisq;
  first = 0;
  while (lambda > 1e-100 && lambda < HUGE_VAL && chisq>0 && (first<2 || ochisq<=chisq || (ochisq-chisq)>0.001) ) {

    if( ochisq > chisq )
      ochisq = chisq;
    
    /* Setup linear equation for in-place solution.  
       mat = alpha*(1+lambda)I.  tryp = beta.
       Fill in entire mat matrix from the symmetric info in alpha. */
    for(j=0;j<NPARAM;j++) {
      for(k=0;k<j;k++) mat[j][k] = alpha[k][j];
      for(k=j;k<NPARAM;k++) mat[j][k] = alpha[j][k];
      mat[j][j] *= (1+lambda);
      tryp[j]=beta[j];
    }
  
    /* Solve the linear equation for a new potential parameter set 'tryp' */
    noerr = lu_decomp(mat,NPARAM,ps,scales);
    if(noerr)
      lu_solve(mat,tryp,NPARAM,ps);
    else { 
      /* Singular matrix found. Quit as soon as sufficient info is printed. */
      chisq = HUGE_VAL;
      break;
    }
    for(j=0;j<NPARAM;j++) tryp[j] += p[j];
 

    /* Compute goodness of fit (and curvature/gradient) at this new location */
    partials(x, y, std, n, tryp, dydp, mat, tryb, &newchi);

    /*
    MTFPRINT8("Chi = %g  A=%g B=%g K=%g Tm=%g dT=%g   lambda=%g\n",
	    sqrt(newchi), tryp[3], tryp[3] ,tryp[0],-tryp[1]/tryp[2],
	    log(81)/tryp[2], lambda);
    */

    /* Is the fit better ? */
    if(newchi < chisq) {
      /* If yes, then decrease lambda (less gradient contribution)
	 and keep the parameter set, gradiant and curvature info */
      if (chisq - newchi < 0.001) first++;  /* Count number of very small
					       increments */
      lambda *= 0.1;
      chisq = newchi;
      for(j=0;j<NPARAM;j++) {
	p[j] = tryp[j];
	beta[j] = tryb[j];
	for(k=j;k<NPARAM;k++) alpha[j][k] = mat[j][k];
      }
    }
    else /* Increase lambda and try again from previous starting point */
      lambda *= 9.;
  }


  /* Set output parameters to their final values */
  *a = p[3];
  *kappa = p[0];
  *m = p[2];
  *c = p[1];

  /* Free storage space */
  free(p);
  free(beta);
  free(alpha[0]);
  free(alpha);
  free(mat[0]);
  free(mat);
  free(tryb);
  free(tryp);
  free(dydp);
  
  return(chisq);
  
}

/*
 * Perform the S-curve fit.
 * Compile data that requires fitting.  Compute a first estimate
 * of parameters, and then perform non-linear regression to approach
 * a reasonable solution.  
 */
double scurve_fit(void)
{
  int i,num_points;
  double *x, *y, *z;
  double a, b, k, dt, min, max, m, c, chisq;
  int *ind;
  double ym, ylo, yhi;
  double tm, tlo, thi;
  int ilo=0, im=0, ihi=0;
  double min_std;
  
  /* Put points in order of increasing reflectance */
  set_up_sort_reflectance_pairs();
  
  num_points = number_unique_density_patches;

  x = (double *)malloc(num_points*sizeof(double));
  y = (double *)malloc(num_points*sizeof(double));
  z = (double *)malloc(num_points*sizeof(double));
  ind = (int *)malloc(num_points*sizeof(int));

  /* Set up x,y,z arrays of ordered data.  Find min/max grey */
  min = sort_reflectance[0].grey;
  max = min;
  min_std = HUGE_VAL;
  for(i=0;i<num_points;i++) {
    x[i] = sort_reflectance[i].reflectance;
    y[i] = sort_reflectance[i].grey;
    ind[i] = (int)sort_reflectance[i].index;
    z[i] = g_test_pattern[ind[i]].std_deviation; 
    if (min_std > g_test_pattern[ind[i]].std_deviation && g_test_pattern[ind[i]].std_deviation > 0.0)
    	min_std = g_test_pattern[ind[i]].std_deviation;
    if (y[i] > max) max = y[i];
    else if (y[i] < min) min = y[i];
  }
  min_std /= 10.0;
  
  for(i=0;i<num_points;i++) 
  	if (z[i] == 0.0) z[i] = min_std;
  	
  /* Set initial values for a, k, tm, and dt */
  /* a = baseline offset of entire logistics curve
     k = max logistics height
     tm = logistics half-way point (inflection)
     dt = distance between 10% and 90% logistics point
  */
  a = min;
  k = max - min;
  ym = 0.5*k + a;
  ylo = 0.1*k + a;
  yhi = 0.9*k + a;
  for(i=0;i<num_points-1;i++) {
    if( y[i] <= ym && y[i+1] > ym ) im = i;
    if( y[i] <= ylo && y[i+1] > ylo ) ilo = i;
    if( y[i] <= yhi && y[i+1] > yhi ) ihi = i;
  }
  
  tm  = x[im]  + (x[im+1]-x[im])/(y[im+1]-y[im])*(ym-y[im]);
  tlo = x[ilo] + (x[ilo+1]-x[ilo])/(y[ilo+1]-y[ilo])*(ylo-y[ilo]);
  thi = x[ihi] + (x[ihi+1]-x[ihi])/(y[ihi+1]-y[ihi])*(yhi-y[ihi]);
  dt = thi-tlo;

  /* Find approx y-intercept of first 'linear' segment of curve */
  if(ilo<3) ilo = 3;
  least_squares(x,y,0,(ilo+1)/2,&b,&a);

  /* Reset a & k to account for logistic values at 0 and 1 */
  a -= k/(1.+exp(log(81.)/dt*tm));
  k = (max - a /* - b */) * (1.+exp(-log(81.)/dt*(1.-tm)));

  m = log(81.)/dt;
  c = -m*tm;

  /* Do non-linear regression */
  chisq = levmarq(x, y, z, num_points, &a, &b, &k, &m, &c );

  if ( chisq < HUGE_VAL ) {
    /* A solution was found, so set working points. */
    for(i=0;i<num_points;i++) {
      g_test_pattern[ind[i]].Scurve = a + k/(1+exp(-(m*x[i]+c)));
    }

    /* Set static parameter values, so the function 
       compute_scurve_greylevel_to_reflectance has access */
    g_a = a;
    g_k = k;
    g_m = m;
    g_c = c;

    /* Uncomment this if you want to see model parameters/goodness of fit in MTFOUT 
    MTFPRINT2("\nS-CURVE FIT with chi-squared goodness of fit = %g\n", chisq)
    MTFPRINT5("ImgGray = %.2f + %.2f/{1 + exp[-log(81)/%.2f(TgtReflect - %.2f)]}\n", 
	     g_a, g_k, log(81.)/g_m, -g_c/g_m)
    MTFPRINT("\nS-CURVE FIT (ImageGray vs TgtReflectance):\n")
    MTFPRINT5("offset = %6.2f  kappa = %6.2f  delta_t = %5.3f  tm = %5.3f\n", 
	     g_a, g_k, log(81.)/g_m, -g_c/g_m)
    MTFPRINT2("chi-squared goodness of fit = %g\n", chisq)
    */
  }
  else {
    MTFPRINT("ERROR:  S-curve Fit failed.  Program unable to continue.");
    MTFPRINT("Please select a different curve fit option");
    /* Program will exit as soon as the density point data is plotted.*/
  }
  
  free(x);
  free(y);
  free(z);
  free(ind);

  return (chisq);
  
}

/* 
 *  Linear regression computation.  This is used for initial computation of the 
 *  logistic constant offset.
 */
void least_squares(double *x, double *y, int first, int n, double *slope, double *intercept)
{
  unsigned int i;
  double sumxy,meanx,meany;
  double meanxsq,sumxsq;
  double sumysq,sumy;
  /* double meanysq, coef_of_det,coef_of_cor; */

  sumy = sumysq = sumxy = meanx = sumxsq = 0.0;
  for (i=first;i<first+n; i++) {
    sumy += y[i]; sumysq += y[i]*y[i];
    sumxy += x[i]*y[i]; meanx += x[i];
    sumxsq += x[i]*x[i];
  }

  meanx /= (double)n;
  meany = sumy/(double)n;
  meanxsq = meanx*meanx;

  *slope = (sumxy - n*meanx*meany) /
    (sumxsq - n*meanxsq);

  *intercept = meany - *slope*meanx;

  /*
  meanysq = meany*meany;
  coef_of_det = (*intercept*sumy + *slope*sumxy - n*meanysq) /
                   (sumysq - n*meanysq);

  coef_of_cor = sqrt(coef_of_det);

  fprintf(stderr, "coef_det=%g  coef_cor=%g  slope=%g inter=%g\n",
	  coef_of_det, coef_of_cor, *slope, *intercept);
  */
}


/* 
 *  Given the scurve model parameters, compute the reflectance at 
 *  the given grey_level.  Reflectance values are clipped to be within 
 *  the range [0,1].
 */

double compute_scurve_grey_to_reflectance (double grey_level)
{
  double val;
  if(grey_level <= g_a) return(0.0);
  val = g_k / (grey_level - g_a);
  
  if (val > 1) 
    val = -(log(val-1.) + g_c)/g_m;
  else
    val = 1;
  
  if (val > 1)  val = 1;
  if (val < 0)  val = 0;
  
  return (val);
}


/* The following code very closely follows the code found at 
   http://www.graphviz.org/pub/graphviz/CURRENT/doxygen/html/lu_8c-source.html
   A few minor changes were made when integrating with the above source. */

/*
 * This code was (mostly) written by Ken Turkowski, who said:
 *
 * Oh, that. I wrote it in college the first time. It's open source - I think I
 * posted it after seeing so many people solve equations by inverting matrices
 * by computing minors naïvely.
 * -Ken
 *
 * The views represented here are mine and are not necessarily shared by
 * my employer.
        Ken Turkowski                   turk@apple.com
        Immersive Media Technologist    http://www.worldserver.com/turk/
        Apple Computer, Inc.
        1 Infinite Loop, MS 302-3VR
        Cupertino, CA 95014
 */

/* Function calls modified slightly to allow temporary working space
 * to be passed in, rather than allocated at each iteration.
 */

/* This module solves linear equations in several variables (Ax = b) using
 * LU decomposition with partial pivoting and row equilibration.  Although
 * slightly more work than Gaussian elimination, it is faster for solving
 * several equations using the same coefficient matrix.  It is
 * particularly useful for matrix inversion, by sequentially solving the
 * equations with the columns of the unit matrix.
 *
 * lu_decompose() decomposes the coefficient matrix into the LU matrix,
 * and lu_solve() solves the series of matrix equations using the
 * previous LU decomposition.
 *
 *      Ken Turkowski (apple!turk)
 *      written 3/2/79, revised and enhanced 8/9/83.
 */

/* lu_decomp() decomposes the coefficient matrix A into upper and lower
 * triangular matrices, the composite being the LU matrix.
 *
 * The arguments are:
 *
 *      a  - the (n x n) coefficient matrix, and output LU matrix
 *      n  - the order of the matrix
 *
 *  1 is returned if the decomposition was successful,
 *  and 0 is returned if the coefficient matrix is singular.
 */

int lu_decomp(double **a, int n, int *ps, double *scales)
{
  register int i, j, k;
  int pivotindex, indx;
  double pivot, biggest, mult, tempf;
 

  pivotindex = 0;
  for (i = 0; i < n; i++) {       /* For each row */
    /* Find the largest element in each row for row equilibration */
    biggest = 0.0;
    for (j = 0; j < n; j++)
      if (biggest < (tempf = fabs(a[i][j])))
	biggest  = tempf;
    if (biggest != 0.0)
      scales[i] = 1.0 / biggest;
    else {
      return(0);      /* Zero row: singular matrix */
    }
    ps[i] = i;          /* Initialize pivot sequence */
  }
  
  for (k = 0; k < n-1; k++) {     /* For each column */
    /* Find the largest element in each column to pivot around */
    biggest = 0.0;
    for (i = k; i < n; i++) {
      indx = ps[i];
      if (biggest < (tempf = fabs(a[indx][k]) * scales[indx])) {
	biggest = tempf;
	pivotindex = i;
      }
    }
    if (biggest == 0.0)
      return(0);      /* Zero column: singular matrix */
    
    indx = ps[k];
    if (pivotindex != k) {      /* Update pivot sequence */
      j = indx;
      ps[k] = ps[pivotindex];
      ps[pivotindex] = j;
    }

    /* Pivot, eliminating an extra variable  each time */
    indx = ps[k];
    pivot = a[indx][k];
    for (i = k+1; i < n; i++) {
      a[ps[i]][k] /= pivot;
      mult = a[ps[i]][k];
      if (mult != 0.0) {
	for (j = k+1; j < n; j++)
	  a[ps[i]][j] -= mult * a[indx][j];
      }
    }
  }

  if (a[ps[n-1]][n-1] == 0.0)
    return(0);          /* Singular matrix */
  return(1);
}

/* lu_solve() solves the linear equation (Ax = b) after the matrix A has
 * been decomposed with lu_decompose() into the lower and upper triangular
 * matrices L and U.   Computations are performed in place so the solution
 * is returned in b.
 *
 * The arguments are:
 *
 *      lu - the lu decomp array created by lu_decomp
 *      b - the constant vector, replaced by solution
 *      n - the order of the equation
 *      ps - pivot permutation vector
*/

void lu_solve(double **lu, double *b, int n, int *ps)
{
  register int i, j;
  double dot;

  /* Vector reduction using U triangular matrix */
  for (i = 0; i < n; i++) {
    dot = b[ps[i]];
    b[ps[i]] = b[i];
    for (j = 0; j < i; j++)
      dot -= lu[ps[i]][j] * b[j];
    b[ps[i]] = b[i];
    b[i] = dot;
  }

  /* Back substitution, in L triangular matrix */
  for (i = n-1; i >= 0; i--) {
    dot = b[i];
    for (j = i+1; j < n; j++)
      dot -= lu[ps[i]][j] * b[j];
    b[i] = dot/lu[ps[i]][i];
  }
}
