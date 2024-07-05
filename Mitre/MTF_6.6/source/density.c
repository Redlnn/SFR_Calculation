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

/*****************************************************************************/
/*                                                                           */
/* Revisions to density.c                                                    */
/*                                                                           */
/* V1.1- replaced some code with macros and did general clean-up             */
/*                                                      brp 10-MAR-94        */
/*                                                                           */
/* V1.3- replaced nint with NINT and M_PI with PI; made change in            */
/*      compute_densities to only include r1 once for negative skew          */
/*      case (see code).  Now cast as double all ints in rc and lc calcs.    */
/*      Also, use floor and ceil functions for rc and lc values respectively.*/
/*      The above corrected a error of up to 3 pixels in the row endpoints.  */
/*                                                      brp 15-JUN-94        */
/*                                                                           */
/* V2.0- A common endpoint routine, next_ep_pair, does all endpoint calcs    */
/*       for both density and sine routines regardless of skew.              */
/*                                                      brp 8-AUG-94         */
/*                                                                           */
/*      g_low_throw_out and g_high_throw_out are set in compute_density_stats*/
/*                                                      brp 15-AUG-94        */
/*                                                                           */
/* V2.3 Allow 0 and 255 pixel values for the black and white patches         */
/*      (F and I)                                                            */
/*      Revised output printout.                                             */
/*      Revised problem report printout.                                     */
/*                                                      DJB 6/7/95           */
/*                                                                           */
/* V3.0 Added Lagrangian interpolation for use when correlation coefficient  */
/*      of linear regression on grey patches is less than CORR_LOWER_LIMIT.  */
/*      Indentation changes.                                                 */
/*                                                      DJB 8/16/95          */
/*      Conditionalize on USE_TIFF to not depend on tiff code.               */
/*                                                      DJB 12/28/95         */
/*      Fix problem with extra carriage returns in printout.                 */
/*                                                      DJB 1/24/96          */
/*                                                                           */
/* V3.1 Removed copyright and added data rights notice.                      */
/*                                                      DJB 6/25/96          */
/*      Added local linear interpolation for use when correlation coefficient*/
/*      of linear regression on grey patches is less than CORR_LOWER_LIMIT.  */
/*      Commented out Lagrangian interpolation stuff.                        */
/*                                                      DJB 8/2/96           */
/*                                                                           */
/* V4.0 Removed low dynamic range check: check not useful in practice.       */
/*      Use density values to determine which patch is white and black,      */
/*               do not rely on labels 'F' and 'I'                           */
/*      General cleanup to adhere to gcc's -ansi, -pedantic, and -Wall       */
/*               switches                                                    */
/*      Removed factory and computed reflectance from verbose printout       */
/*                                                      WE 9/99              */
/*                                                                           */
/*                                                                           */
/*  V4.1 Removed extrapolation beyond gray patch endpoints when nonlinear    */
/*       I/O curve following is invoked via option C or P (sometimes produced*/
/*       erratic extrapolated values); is now clamped to endpoint value.     */
/*                                                      nbn, 9/00            */
/*                                                                           */
/*  V5.0 S-curve fitting option added. Spline fits now revert to s-curve     */
/*                                     if there is an error.                 */
/*                                                      mal, 11/03           */
/*                                                                           */
/*  V5.5    Fixed a bug introduced in v5.0 which prevented using linear      */
/*          regression when fewer than 5 density patches were present.       */
/*          Made alias problems report even when option 'n' selected         */
/*          Added bartarget handling to print_density_data.                  */
/*          Corrected AVGing code, so that 4-patch AVG is retained if other  */
/*          pairs meet averaging threshold too (other pairs not averaged)    */
/*                                                      mal, 02/04           */
/*                                                                           */
/* V5.6     All extended output written to MTFOUT only, not stdout.          */
/*                                                      mal, 06/04           */
/*                                                                           */
/* V5.8     Cleaned up code.  Sections no longer used were removed.          */
/*                                                      mal, 09/04           */
/*                                                                           */
/* V6.2     Added non-linear reflectance data for bartarget mode   	     */
/*                                                      MAL, 08/06           */
/*                                                                           */
/* V6.3     Don't force reflectances to [0,1] for bartarget mode   	     */
/*          Only perform lin-reg bound check for Appendix F spec. 	     */
/*          Remove bartarget dependencies in print_density_data since        */
/*          it's never called on bartargets.                                 */
/*                                                      MAL, 10/06           */
/* V6.5    Added colors to debug image and removed full area fill.           */
/*                                                      AWU, 07/14           */
/*****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#define TIFF char   /* To allow mtf.h processing without hiccup */
#include "mtf.h"
#include "xterndcl.h"

static void set_up_grey_reflectance_pairs(void);
static int compare_grey_reflectance_pairs(const void *, const void *);

typedef struct Grey_Reflectance_Pair
{
  double grey;
  double reflectance;
} grey_reflectance_pair;

static grey_reflectance_pair *grey_reflectance_pairs;
static int number_grey_pairs;

/******************************************************************************/
/*                                                                            */
/******************************************************************************/

/* return successive lc, rc scanline enpoints for a patch. With each successive
   call we increment the global row counter (g_rc) and compute a global g_lc,
   g_rc pair of endpoints for the row. We return a one status until we have
   exhausted the patch, at which time we return zero.  When zero is returned
   the g_lc, g_rc pair is invalid. The caller initializes for a new patch by
   assigning the patch corners to g_R1, g_C1, etc. thru G_ASSIGN_CORNERS and
   sets g_cr to -1 before calling next_ep_pair.
*/

/* some macros to keep things readable (more or less) */

#define J     /* g_R1, g_C1, etc. are assigned by G_ASSIGN_CORNERS in mtf.h */ \
  (((g_test_pattern[i].height_in_mm/MM_PER_INCH)*g_ppi_horz*                   \
    (1.0-2.0*(double)BORDER)) - (double)(g_cr-g_R1)/cos(g_skew_avg))

#define RC_CALC_FOR_POS_SKEW                                                  \
  if (g_cr< g_R2) {                                                           \
    g_rc=(int)floor((double)g_C1 + (double)(g_cr-g_R1)/tan(g_skew_avg));      \
    if (g_rc > g_C2) g_rc = g_C2;                                             \
  }                                                                           \
  else {                                                                      \
    g_rc = (int)floor((double)g_C1 + g_l*cos(g_skew_avg) -                    \
                 (double)(g_cr-g_R2)*tan(g_skew_avg));                        \
    if (g_rc > g_C2) g_rc = g_C2;                                             \
  }

#define RC_CALC_FOR_NEG_SKEW                                                  \
  if (g_cr< g_R4) {                                                           \
    g_rc = (int)floor((double)g_C1 + g_l*cos(g_skew_avg) -                    \
                 (double)(g_cr-g_R2)*tan(g_skew_avg));                        \
    if (g_rc > g_C4) g_rc = g_C4;                                             \
  }                                                                           \
  else {                                                                      \
    g_rc = (int)floor((double)g_C1-(double)(g_cr-g_R1)*tan(g_skew_avg) -      \
                 (double)(J)/sin(g_skew_avg));                                \
    if (g_rc > g_C4) g_rc = g_C4;                                             \
  }

short next_ep_pair(register int i)
{
  /* the best way to see what's going on here is to draw a square, label its
     corners with g_R1,g_C1 at the upper left, g_R2,g_C2 at upper right, g_R3,g_C3
     at lower left and then skew it- you'll see why this gets complicated - 
     each scan line starts and ends at a unique position.  We needed to do
     it this way because rotating and translating the image to a non-skewed,
     normalized orientation would have left some holes and caused some pixels
     to map to the same location. */

  if (g_cr == -1) {
    if (g_skew_avg >= 0)
      g_cr = g_R1;
    else
      g_cr = g_R2;
  }
  else g_cr++;


  if (g_skew_avg >= 0) {
    if (g_cr > g_R4) return(0);
    if (g_cr<g_R3) {
      g_lc = (int)ceil((double)g_C1-
                  (double)(g_cr-g_R1)*tan(g_skew_avg));
      if (g_lc < g_C3) g_lc = g_C3;
      RC_CALC_FOR_POS_SKEW
    }
    else {
      g_lc = (int)ceil((double)g_C3+
                  (double)(g_cr-g_R3)/tan(g_skew_avg));
      if (g_lc < g_C3) g_lc = g_C3;
      RC_CALC_FOR_POS_SKEW
    }
  }
  else {                        /* negative skew cases */
    if (g_cr >= g_R3) return(0);
    if (g_cr<g_R1) {
      g_lc = (int)ceil((double)g_C1-
                  (double)(g_R1-g_cr)/tan(g_skew_avg));
      if (g_lc < g_C1) g_lc = g_C1;
      RC_CALC_FOR_NEG_SKEW
    }
    else {
      g_lc = (int)ceil((double)g_C1-
                  (double)(g_cr-g_R1)*tan(g_skew_avg));
      if (g_lc < g_C1) g_lc = g_C1;
      RC_CALC_FOR_NEG_SKEW
    }
  }
  return(1);
}


/******************************************************************************/
/*                                                                            */
/******************************************************************************/

/* look at all the density patches and compute linear regression line */

void compute_densities(unsigned int white_patch, unsigned int black_patch)
{
  register unsigned int i, j;
  register int k;
  int monotone, counter, number_density_patches;
  int *val_array;           /* histogram of a density patch stored here */
  double chisq, spline_chisq;

  val_array = (int *)calloc((g_max_pixel_value+1),sizeof(int));

  number_density_patches = 0;
  for (i=0;i<g_number_of_patches;i++) {
    if (g_test_pattern[i].patch_type == DENSITY_PATCH) {

	number_density_patches = number_density_patches + 1;
      /* the vertices of the box 10% within the density patch */

      G_ASSIGN_CORNERS          /* g_R1 = g_test_pattern[i].r[0],etc.*/

      g_cr = -1;                /* tell next_ep_pair to initialize */
      for (j=0;j<=g_max_pixel_value;j++) val_array[j]=0; /* histogram */
      /*g_l is the unskewed 80% patch width in pixels. 80% because of the .1 border*/
      /* note that if the value of BORDER changes from .1 this should still work */
      g_l = (g_test_pattern[i].width_in_mm/MM_PER_INCH)*g_ppi_horz*
        (1.0-2.0*(double)BORDER);


      while(next_ep_pair(i)) {
	/* make a histogram of the patch on which we will do statistics
	   this is easy, accurate and flexible - see compute_density_stats()
	   to see how it works  */

	switch (g_bytes_per_pixel) {
	  unsigned short *sbuf;
	  /*the following commented out code writes the filled in boxes on the density patches. since we now draw a border on the calculated region, the code is no longer necessary*/
	case 1:
	  for(k=g_lc;k<=g_rc;k++) {
	    val_array[g_image_array[g_cr*g_total_image_width+k]]++;
	    /*if (g_debug)
	      g_image_array_modified[g_cr*g_total_image_width+k] = g_max_pixel_value;*/
	  }
	  break;
	case 2:
	  sbuf = (unsigned short *)g_image_array;
	  for(k=g_lc;k<=g_rc;k++) {
	    val_array[sbuf[g_cr*g_total_image_width+k]]++;
	    /*if (g_debug)
	      g_image_array_modified[g_cr*g_total_image_width+k] = 255;*/
	  }
	  break;
	}
      }
      /* If we have an image with a break in the middle (from a target scanned in
         sections on a small scanner and pieced together) then we must dump any
         data from the middle break.  This is done for ALL patches.
         The break may be 255 or 0 - dump both values.  If a break is 
         detected, then the data will be dumped from the break  */

      /* Allow 0 and 255 pixel values for the black and white patches.
         DJB 6/7/95 */

      if ((i != white_patch) && (i != black_patch))
        {
          val_array[0]=0;
          val_array[g_max_pixel_value]=0;
        }

      compute_density_stats(val_array,i,white_patch,black_patch);
    }
  }

  free(val_array);

/* 
   Combine all density patches that have nearly the same densities,
   i.e., average the reflectances and average the grey levels; this new set of
   reflect-grey values is then used in curve fit.
   In any given tgt, either 2, 3, or 4 tgt patches may have nearly the same 
   density (by design). BUT if following algorithm finds 2 or more SEPARATE 
   groups of equidensity patches then NO patches are averaged, because it is 
   assumed a tgt doesnt contain such separate groups
   (i.e., the patches were meant to be different, even if nearly the same). 
*/

  {
  double xmax, xmin, x1, x2, r_ratio; 
  counter = number_unique_density_patches = number_of_averages = 0;
  avg_index = (int **)calloc(10,sizeof(int *));

  if(number_density_patches < 17) 
  	r_ratio = R_RATIO1;	  /* R_RATIO defined in mtf.h */
  else  
  	r_ratio = R_RATIO2;
   
  for (i=0;i<g_number_of_patches;i++) {
    if(g_test_pattern[i].patch_type == DENSITY_PATCH) {
      counter = 1;              
      number_unique_density_patches++;  
      g_test_pattern[i].final_reflectance = g_test_pattern[i].reflectance; 
      g_test_pattern[i].final_grey = g_test_pattern[i].grey_level;
      for (j=i+1;j<g_number_of_patches;j++) {
        if(g_test_pattern[j].patch_type == DENSITY_PATCH) {
          x1 = g_test_pattern[i].reflectance;           
          x2 = g_test_pattern[j].reflectance;           
          if(x1 < x2) {
            xmax = x2;
            xmin = x1;
          }                                                             
          else {
            xmax = x1;
            xmin = x2;
          }
          if ((xmax / xmin) < r_ratio) {               
            /* MTFPRINT4("%s %s  %f ", g_test_pattern[i].patch_label, 
               g_test_pattern[j].patch_label, xmax/xmin); */
            if ( counter == 1 ) {
              number_of_averages++;
              strncpy( g_test_pattern[i].new_patch_label, "AVG\0", 4);
              avg_index[number_of_averages] = (int *)calloc(10,sizeof(int));
              avg_index[number_of_averages][1] = i;
            }
            g_test_pattern[i].final_reflectance =
              (g_test_pattern[i].final_reflectance * counter + 
               g_test_pattern[j].reflectance )
              / (double)(counter+1);
            g_test_pattern[i].final_grey =
              (g_test_pattern[i].final_grey * counter +
                                      g_test_pattern[j].grey_level) 
              / (double)(counter+1);
            g_test_pattern[j].patch_type = DUPLICATE_DENSITY;
            counter += 1;
            avg_index[number_of_averages][0] = counter;
            avg_index[number_of_averages][counter] = j;

          }                                             
        }
      }
    }
  }
   

  /* If > 1 AVG set was found, reset the averaged array sets back to original, 
     unaveraged values.  Assumption is that tgt never contains more than one true
     set of equidensity patches, so if multiple sets are found, their constituents
     are assumed to be Different densities by design, not equidensity, 
     so avoid averaging entirely if multiple sets are found */
  /* But avoid clearing any solo 4-pt AVG v5.5 */
     if(number_of_averages > 1) {
       int keep = 0;
       
       counter = 0;
       for (i=1; i<=number_of_averages; i++) {
	 if(avg_index[i][0]==4) { 
	   counter++;
	   keep = i;
	 }
       } 
       /* Anything other than exactly one 4-pt average was found, then
	  we won't keep any of the averages */
       if (counter != 1) keep = 0;
       
       for (i=1; i<=number_of_averages; i++) {
	 if (i!=keep) { 
	   /* Revert to original patch content */
	   for (j=1;j<=avg_index[i][0];j++) { 
	     k = avg_index[i][j];
	     g_test_pattern[k].final_reflectance = 
	       g_test_pattern[k].reflectance; 
	     g_test_pattern[k].final_grey = g_test_pattern[k].grey_level;      
	     strncpy(g_test_pattern[k].new_patch_label, 
		     g_test_pattern[k].patch_label, 4);
	     g_test_pattern[k].patch_type = DENSITY_PATCH; 
	   }
	   number_unique_density_patches += avg_index[i][0]-1;
	 }
       }
       /* Set the number_of_averages back to correct value */
       if (keep > 0) number_of_averages = 1;
       else number_of_averages = 0;
     }
	
/* Identify Min & Max TgtReflectances; used in local_linear_interpolation */   
     minreflect = 1.1;
     maxreflect = -0.1;
     for (i=0; i < number_unique_density_patches-1;i++) {
       if(g_test_pattern[i].final_reflectance < minreflect)
	 minreflect = g_test_pattern[i].final_reflectance;
       if(g_test_pattern[i].final_reflectance > maxreflect)
	 maxreflect = g_test_pattern[i].final_reflectance;  
     }
  }
  
  /* do the linear regression */
  regress();


  /* 
  If doing printer MTF and no specific fitting is forced, use S-curve fit. Else,
  do the smoothed & unsmoothed spline curve fits, then determine which is best,
  either linreg, smoothspline, or unsmoothspline. Unsmoothed spline goes thru 
  all input data points; smoothed spline smooths the input data before applying
  spline.  Compare the results to an S-curve fit via chi-square goodness of
  fit.  The best of these 4 (linreg, smooth/unsmooth spline, s-curve) becomes
  the chosen fit.

  Initial value of g_followcurve is defined in mtf.c, based upon g_user_requestcurve.
  if initial g_followcurve = g_user_requestcurve = 0, then program determines best curve to use; 
  if initial g_followcurve = 5,6,7,8 or 9 then user or previous options set the 
                                          appropriate curve type to use
  
  if initial g_followcurve = 0 then it is reset as follows: 
     g_followcurve = 6   program using linreg 
     g_followcurve = 7   program using unsmooth spline
     g_followcurve = 8   program using smooth spline
     g_followcurve = 9   program using s-curve fit

  If initial g_followcurve = g_user_requestcurve = 7 or 8, there is a chance
  that the spline will not be able to fit the data monotonically.  When this
  occurs, g_followcurve is reset to 9 and the s-curve fit is used instead.

  The piecewise curve fit is only available as a last resort via manual user
  setting.    *** New in V5.0.  Previously piecewise linear was the fallback
  when the spline fit didn't work due to non-monotonicity.  ***
     g_followcurve = 5     user selected piecewise linear fit 
  */
  
  chisq = 0.0;
  if (g_followcurve == 9) {
    /* Only do the s-curve fit */
    chisq = scurve_fit();
  }
  else if (g_followcurve != 6) {
    /* Do the spline fit and check against the s-curve fit */
    curve_type = smooth_spline(&monotone, &spline_chisq);
    chisq = scurve_fit();
    if (curve_type == 9) 
      g_followcurve = 9;
    else if(g_followcurve == 0 && chisq < spline_chisq)
      if (chisq < spline_chisq) curve_type = 9;
  }

  if (number_unique_density_patches < 5 && g_followcurve > 6) {	
    /* not enough density patches to perform spline or s-curve fit, force linear */
    g_followcurve = 6;
    for(i=0; i<number_unique_density_patches; i++) {
      g_test_pattern[i].Sspline = 0.0;
      g_test_pattern[i].Scurve = 0.0;
    } 
  }
  else if(number_unique_density_patches > 4) {
    /* program determines best curve type to use and g_followcurve is reset */
    if (g_followcurve == 0) 
      g_followcurve = curve_type;
    if(g_followcurve == 7 || g_followcurve == 8)
      invert_spline();
  }

  print_density_data();

  if( g_followcurve < 5 || g_followcurve > 9) {
    MTFPRINT2("\nWARNING: Unknown curve fit (%d) applied. Should be in range [5,9]\n", g_followcurve)
    MTFPRINT("           LINEAR REGRESSION will be used instead!\n")
    g_followcurve = 6;
  }
  
  if( g_followcurve == 5 && !monotone ) {
    /* Piecewise Linear Interpolation used on non-monotone data */
    MTFPRINT("\nWARNING: MeasImgGray vs. TgtReflectance data is not monotonic and this\n")
    MTFPRINT("erratic data is not smoothed by piecewise linear interp.\n")
  }
  
  if (chisq >= HUGE_VAL) /* Non-recoverable error, (caused by a singular
			matrix during s-curve fit. Exit now that
		        density data has been printed. */
      mtfexit(-1);
  
} 



/******************************************************************************/
/*                                                                            */
/******************************************************************************/

void print_density_data(void)
{
  char problem_string[82];
  unsigned int i, j;
  double lrig,ler, patchdensity;
  int lin_reg_problem;


  MTFPRINT("\n")
  if(g_followcurve==9 || (g_printer && g_followcurve <8))
     MTFPRINT3("%s%s\n","Density   Tgt     Tgt       Meas.     ",
            "S-CurveFit   LinReg    LinReg")
  else
    MTFPRINT3("%s%s\n","Density   Tgt     Tgt       Meas.     ",
              "SmthSpline   LinReg    LinReg")
  MTFPRINT3("%s%s\n","Patch   Density   Reflect   ImgGray   ",
            "ImgGray      ImgGray   Error")

  lin_reg_problem = 0;
  for (i=0;i<g_number_of_patches;i++) {
    if (g_test_pattern[i].patch_type == DENSITY_PATCH) {
      lrig = g_test_pattern[i].final_reflectance*g_slope +
        g_intercept;
      ler = g_test_pattern[i].final_grey - lrig;
      if (g_compare==APPF && g_bytes_per_pixel== 1 && fabs(ler) - (double)DENSITY_MAX_LERR > 0) {
        sprintf(problem_string,
                "The lin reg error for Density Patch %s exceeds %.2f\n",
                g_test_pattern[i].new_patch_label,
                (double)DENSITY_MAX_LERR);
        put_problem(problem_string, IQS);
        lin_reg_problem = 1;
      }
     
     if(i == 0) 
     	patchdensity = log10(1. / g_test_pattern[0].final_reflectance);
     else
     	patchdensity = g_test_pattern[i].patch_density;
     
     if(g_followcurve==9 || (g_printer && g_followcurve <8))
	MTFPRINT8("%-9.7s%5.3f%9.3f%12.3f%11.3f%11.3f%10.3f\n"
                ,g_test_pattern[i].new_patch_label,
                patchdensity,
                g_test_pattern[i].final_reflectance,
                g_test_pattern[i].final_grey,
                g_test_pattern[i].Scurve, 
                lrig,
                ler)
     else
	MTFPRINT8("%-9.7s%5.3f%9.3f%12.3f%11.3f%11.3f%10.3f\n"
                ,g_test_pattern[i].new_patch_label,
                patchdensity,
                g_test_pattern[i].final_reflectance,
                g_test_pattern[i].final_grey,
                g_test_pattern[i].Sspline, 
                lrig,
                ler)
    }
  }
  
  MTFPRINT2("%s","ImgGray-to-TgtReflect curve used: ");
  if(g_user_requestcurve == 0) {
    if(g_followcurve == 6) 
      MTFPRINT2("%s\n","Linear Regression (program selected)")
    if(g_followcurve == 7) 
      MTFPRINT2("%s\n","Unsmoothed Spline (program selected)")
    if(g_followcurve == 8) 
      MTFPRINT2("%s\n","Smoothed Spline (program selected)")
    if(g_followcurve == 9 && g_printer) 
      MTFPRINT2("%s\n","S-curve Fit (printer default)")
    else if(g_followcurve == 9 && !g_printer) 
      MTFPRINT2("%s\n","S-curve Fit (program selected)")
  } else if (g_followcurve != g_user_requestcurve) { 
    /* Something caused the program to override user selection */
    /* Typically non-monotonicity */
    if(g_followcurve == 5) 
      MTFPRINT2("%s\n","Piecewise Linear Interp. (program forced)")
    if(g_followcurve == 6) 
      MTFPRINT2("%s\n","Linear Regression (program forced)")
    if(g_followcurve == 7) 
      MTFPRINT2("%s\n","Unsmoothed Spline (program forced)")
    if(g_followcurve == 8) 
      MTFPRINT2("%s\n","Smoothed Spline (program forced)")
    if(g_followcurve == 9) 
      MTFPRINT2("%s\n","S-curve Fit (program forced)")
  } else { /* User gets what he requested */
    if(g_followcurve == 5) 
      MTFPRINT2("%s\n","Piecewise Linear Interp. (user selected)")
    if(g_followcurve == 6) 
      MTFPRINT2("%s\n","Linear Regression (user selected)")
    if(g_followcurve == 7) 
      MTFPRINT2("%s\n","Unsmoothed Spline (user selected)")
    if(g_followcurve == 8) 
      MTFPRINT2("%s\n","Smoothed Spline (user selected)")
    if(g_followcurve == 9) 
      MTFPRINT2("%s\n","S-curve Fit (user selected)")
  }
  
  for( i=1; i<=number_of_averages; i++) {
    MTFPRINT2("%s","AVG=avg of EquiDensityPatches: ");
    for (j=1; j<=(unsigned int)avg_index[i][0]; j++) {
      MTFPRINT2("%s ", g_test_pattern[avg_index[i][j]].patch_label);
    }
    MTFPRINT2("%s\n","     MeasImgGray,TgtDensity pairs:");
    for (j=1; j<=(unsigned int)avg_index[i][0]; j++) {
      MTFPRINT3("%11.3f,%6.3f", g_test_pattern[avg_index[i][j]].grey_level,
      		g_test_pattern[avg_index[i][j]].patch_density); 
    }
    MTFPRINT("\n");
  }
  
  if (lin_reg_problem != 0)
    put_problem("\n", IQS);
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/


/* compute the reflectance linear regression line */
void regress(void)
{
  unsigned int i;
  double sumxy,meanx,meany;
  double meanxsq,meanysq,sumxsq;
  double sumysq,sumy,coef_of_det,coef_of_cor;
  double r,gl;

  sumy = sumysq = sumxy = meanx = sumxsq = 0.0;
  for (i=0;i< g_number_of_patches; i++) {
    if (g_test_pattern[i].patch_type == DENSITY_PATCH) {
      gl = g_test_pattern[i].final_grey;
      r = g_test_pattern[i].final_reflectance;
      sumy += gl; sumysq += gl*gl;
      sumxy += (r*gl); meanx += r;
      sumxsq += r*r;
    }
  }

  meanx /= (double)number_unique_density_patches;
  meany = sumy/(double)number_unique_density_patches;
  meanxsq = meanx*meanx;
  meanysq = meany*meany;

  g_slope = (sumxy - number_unique_density_patches*meanx*meany) /
    (sumxsq - number_unique_density_patches*meanxsq);

  g_intercept = meany - g_slope*meanx;

  coef_of_det = (g_intercept*sumy + g_slope*sumxy - 
                 number_unique_density_patches*meanysq) /
                   (sumysq - number_unique_density_patches*meanysq);

  coef_of_cor = sqrt(coef_of_det);

  
  set_up_grey_reflectance_pairs();


  MTFPRINT("LINEAR REGRESSION (ImageGray vs TgtReflectance):")
  MTFPRINT4("\nslope = %f  intercept = %f  corrcoeff = %f\n",g_slope,g_intercept,coef_of_cor)

}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/

/* this routine is used only if the user wants extended data output */
void print_histo(unsigned int i)
{

  unsigned int low_index,high_index;

  if (!g_extended) return;

  if (g_test_pattern[i].patch_type == DENSITY_PATCH
	|| g_test_pattern[i].patch_type == DUPLICATE_DENSITY) {
    /*  fprintf(g_mtfout,"\nHistogram of density patch after flier removal\n");
	fprintf(g_mtfout, "%s%s\n","(High values are above mean, lows are ",
	                  "equal or below mean)"); */
  }
  else {
    fprintf(g_mtfout,"\nGray Level Histograms at Image Sine Peaks & Valleys:\n");
    
    fprintf(g_mtfout,"\t    PEAK\t\t    VALLEY\n");
    fprintf(g_mtfout,"\tvalue count\t\tvalue count\n");
  }
  
  high_index=g_max_pixel_value;
  low_index=0;

  while(1 /*CONSTANTCONDITION*/) {
    while(g_test_pattern[i].high_histo[high_index]==0
          && high_index) high_index--;
    if(high_index) {
      if(high_index==(unsigned int)(g_test_pattern[i].sine_peak)){
	fprintf(g_mtfout,"\t%5d*%5d",high_index,
		g_test_pattern[i].high_histo[high_index]);
	  }
      else {
	fprintf(g_mtfout,"\t%5d %5d",high_index,
		g_test_pattern[i].high_histo[high_index]);
	  }
      high_index--;
    }
    else fprintf(g_mtfout, "\t\t");
    while(g_test_pattern[i].low_histo[low_index]==0 &&
          low_index <= g_max_pixel_value) low_index++;
    if(low_index <= g_max_pixel_value) {
      if(low_index==(unsigned int)(g_test_pattern[i].sine_valley)) {
	fprintf(g_mtfout, "\t\t%5d*%5d\n",low_index,
		g_test_pattern[i].low_histo[low_index]);
	  }
      else {
	fprintf(g_mtfout,"\t\t%5d %5d\n",low_index,
		g_test_pattern[i].low_histo[low_index]);
	  }
      low_index++;
    }
    else fprintf(g_mtfout,"\n");
    if (high_index <= 0 && low_index > g_max_pixel_value) {
      fprintf(g_mtfout,"\n");
      break;
    }
  }
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/

/* this is called only if the user wants additional data in the printout */
void print_density_stats(void)

{
  unsigned int i;

  if (g_bartarget) return;

  fprintf(g_mtfout, "\n\n\tDENSITY PATCHES: IMAGE GRAY LEVEL STATISTICS\n\n");
  
  for (i=0;i<g_number_of_patches;i++) {
    if (g_test_pattern[i].patch_type == DENSITY_PATCH
	|| g_test_pattern[i].patch_type == DUPLICATE_DENSITY) {
      fprintf(g_mtfout, "Patch = %s\tlow = %d high = %d\n\
\ttotal+fliers = %d total = %d fliers = %d\n\
\tfirst mean = %f final mean = %f\n\
\tstandard deviation after flier removal = %f\n\n",
                 g_test_pattern[i].patch_label,g_test_pattern[i].low,
                 g_test_pattern[i].high,
                 g_test_pattern[i].total+g_test_pattern[i].fliers,
                 g_test_pattern[i].total,g_test_pattern[i].fliers,
                 g_test_pattern[i].first_mean,
                 g_test_pattern[i].grey_level,
	      g_test_pattern[i].std_deviation);
      
    }
  }
}




/******************************************************************************/
/*                                                                            */
/******************************************************************************/

void compute_density_stats(int *val_array,register unsigned int i, unsigned int white_patch, unsigned int black_patch)
{
  /* i is an index into g_test_pattern (the array of patch data */

  int low = 0,high,k;
  int rid_count_high,rid_count_low;
  double m_sum,m_total,sum_x_sqrd,mean;

  high = g_max_pixel_value;
  
  /* val_array has a histogram of the density patch under examination */
  /* first find the lowest and highest values we found in the patch */

  for (k=0;k<=(int)g_max_pixel_value;k++) {
    if (val_array[k] != 0) {
      g_test_pattern[i].low = low = k;
      break;
    }
  }
  for (k=g_max_pixel_value;k>=0;k--) {
    if (val_array[k] != 0) {
      g_test_pattern[i].high = high = k;
      break;
    }
  }


  /* get the total number of pixels we looked at and the sum of their values */
  m_sum=m_total = 0.0;
  for (k=low;k<=high;k++) {
    m_sum += (double)val_array[k]*(double)k;
    m_total += (double)val_array[k];
  }

  g_test_pattern[i].first_mean = m_sum/m_total;
  g_test_pattern[i].fliers = (int)(2.0*(double)THROW_OUT*m_total);

  /* get rid of highest and lowest 16 % of samples */
  /* we already know the sum of everything, just subtract away contributions
     made by the highest and lowest 16%  - remember that val_array is
     a histogram and so values are the index into val_array  (k*val_array[k] =
     the sum of contributions made by pixels of intensity = k) */

  rid_count_high = rid_count_low = NINT((double)THROW_OUT*m_total);
  m_total -= (double)(rid_count_high+rid_count_low);
  g_test_pattern[i].total = NINT(m_total);

  for (k=low;k<=high;k++) {
    if (val_array[k] != 0) {
      rid_count_low -= val_array[k];
      m_sum -= (double)val_array[k]*(double)k;
      /* we went too far, add some back in */
      if (rid_count_low <= 0) {
        val_array[k] = abs(rid_count_low);
        m_sum += (double)val_array[k]*(double)k;
        low=k;
        break;
      }
      else
        val_array[k]=0;
    }
  }
  for (k=high;k>=low;k--) {
    if (val_array[k] != 0) {
      rid_count_high -= val_array[k];
      m_sum -= (double)val_array[k]*(double)k;
      if (rid_count_high <= 0) {
        val_array[k] = abs(rid_count_high);
        m_sum += (double)val_array[k]*(double)k;
        high=k;
        break;
      }
      else
        val_array[k]=0;
    }
  }

  g_test_pattern[i].grey_level = mean = m_sum/m_total;

  /* These values are the cutoffs for sine.c and the peak finders to determine
     any "gaps" in the data introduced by piecing together quadrant sub-images
     of a target that was scanned on a small area scanner */

  /* First set the high value to be above the average value of the "white"
     patch.  We set it halfway between the avg. and 255 */
  if (i == white_patch) {
    g_high_throw_out = g_test_pattern[i].high;
  }
  /* The low value is halfway between the "blackest" grey and zero */
  if (i == black_patch) {
    g_low_throw_out = g_test_pattern[i].grey_level/2.0;
  }

  sum_x_sqrd = 0.0;

  if(g_extended) {
    g_test_pattern[i].high_histo = (int *)calloc(g_max_pixel_value+1,sizeof(int));
    g_test_pattern[i].low_histo = (int *)calloc(g_max_pixel_value+1,sizeof(int));
  }
  
  for (k=low;k<=high;k++) {
    if (val_array[k] !=0) {
      sum_x_sqrd +=
        (double)val_array[k]*((double)k-mean)*((double)k-mean);
      if (g_extended) {
	if (k > mean)
	  g_test_pattern[i].high_histo[k]+=
	    val_array[k];
	else
	  g_test_pattern[i].low_histo[k]+=
	    val_array[k];
      }
    }
  }
  g_test_pattern[i].std_deviation = sqrt(sum_x_sqrd/(m_total-1));
}




/******************************************************************************/
/*   Local Linear Interpolation, aka Piecewise Linear Interpolation:           */
/******************************************************************************/

double local_linear_interpolate_grey_to_reflectance(double value)
     /* Find reflectance value corresponding to given grey value by fitting
        a straight line to the grey-level patches on either side of the
        given grey value; if given grey value is less than min grey patch
        clamp reflectance to lowest reflectance patch, likewise if given grey value
	is greater than max grey patch, clamp relectance to highest reflectance patch.
        The output value is clamped to be between 0 and 1 (inclusive). */
{
  int i, lo = 0, hi = -1;
  double m, b, reflectance = 0.0;

  if (value < grey_reflectance_pairs[0].grey)
    /* Must be '<' comparison to remove bug.  (mal 10/06) */
    /* Supplied grey value is less than avg grey of lowest grey patch, clamp Reflectance
           to lowest reflectance gray patch */
    {
          return minreflect;
    }
  else if (value > grey_reflectance_pairs[number_grey_pairs - 1].grey)
    /* Must be '>' comparison to remove bug.  (mal 08/06) */
    /*Supplied grey value is greater than avg grey of highest grey patch, 
      clamp Reflectance to highest reflectance gray patch */
    {
     return maxreflect;     
    }
  else
    {
      /* Find grey level patches on either side of supplied grey value */
      for (i = 1;  i < number_grey_pairs;  i++)
        {
          if (grey_reflectance_pairs[i].grey >= value)
            {
              lo = i - 1;
              hi = i;
              break;
            }
        }
    }

  /* Sanity check */
  if (hi == -1)
    {
      hi = number_grey_pairs - 1;
      MTFPRINT("Trouble in local_linear_interpolate_grey_to_reflectance\n")
      MTFPRINT7("value = %f  low patch 0 %f %f  high patch %d %f %f\n",
                value,
                grey_reflectance_pairs[0].grey, grey_reflectance_pairs[0].reflectance,
                hi, grey_reflectance_pairs[hi].grey, grey_reflectance_pairs[hi].reflectance)
      mtfexit(-1);
    }

  /* Calculate reflectance based on straight-line interpolation */
  if ( grey_reflectance_pairs[hi].grey != grey_reflectance_pairs[lo].grey ) {
    m = ((grey_reflectance_pairs[hi].reflectance - grey_reflectance_pairs[lo].reflectance)
       / (grey_reflectance_pairs[hi].grey - grey_reflectance_pairs[lo].grey));
    b = (((grey_reflectance_pairs[hi].grey * grey_reflectance_pairs[lo].reflectance)
        - (grey_reflectance_pairs[lo].grey * grey_reflectance_pairs[hi].reflectance))
       / (grey_reflectance_pairs[hi].grey - grey_reflectance_pairs[lo].grey));
    reflectance = m * value + b;
  }
  else {
      MTFPRINT("Zero division in local_linear_interpolate_grey_to_reflectance\n")
      MTFPRINT8("value = %f  patch %d %f %f   patch %d %f %f\n",
                value, lo,
                grey_reflectance_pairs[lo].grey, grey_reflectance_pairs[lo].reflectance,
                hi, grey_reflectance_pairs[hi].grey, grey_reflectance_pairs[hi].reflectance)
      mtfexit(-1);
  }

  if (!g_bartarget) { /* Don't do this for bar targets */
    /* Make sure returned reflectance is between 0.0 and 1.0 */
    if (reflectance > 1.0)
      return 1.0;
    else if (reflectance < 0.0)
      return 0.0;
  }
  return reflectance;

}


void get_grey_reflectance_pairs(int number) {
    grey_reflectance_pairs
    = (grey_reflectance_pair *)calloc(number, sizeof(grey_reflectance_pair));
  
  if (grey_reflectance_pairs == NULL)
    {
      fprintf(stderr, "Error: Could not calloc sufficient space for grey sorted data\n");
      perror("get_grey_reflectance_pairs");
    }

  number_grey_pairs = number;
}

void set_grey_reflectance_pair(int j, double grey, double reflectance) {
  if( j < number_grey_pairs) {
    grey_reflectance_pairs[j].grey = grey;
    grey_reflectance_pairs[j].reflectance = reflectance;
  }
  else {
      fprintf(stderr, "Error: Attempt to access grey sorted data beyond valid data\n");
      perror("set_grey_reflectance_pairs");
  }
}

/* Sort pairs so that min grey level is first, max last */
void sort_grey_reflectance_pairs() {
  qsort(grey_reflectance_pairs, number_grey_pairs,
        sizeof(grey_reflectance_pair), compare_grey_reflectance_pairs);
  /* Set these values so that obvious errors will appear if logic is 
     wrong for bar targets */
  minreflect = 1.1;
  maxreflect = -0.1;
}

int bad_grey_check(double grey, int exitflag) {

  if (grey < grey_reflectance_pairs[0].grey || 
      grey > grey_reflectance_pairs[number_grey_pairs-1].grey) {
    if(exitflag) {
      fprintf(stderr, "Error: accessing grey value %f, beyond nonlinear range [%f,%f] \n", 
	     grey, grey_reflectance_pairs[0].grey,
	      grey_reflectance_pairs[number_grey_pairs-1].grey);
      mtfexit(-1);
    }
    return 1;
  }

  return 0;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/

static void set_up_grey_reflectance_pairs(void)
     /* Set up grey-reflectance pairs for finding reflectance values corresponding
        to grey values via interpolation. */
{
  unsigned int i;
  short j;

  get_grey_reflectance_pairs(g_number_of_patches); 

  /* Initialize grey_reflectance_pairs */
  j = 0;
  for (i = 0;  i < g_number_of_patches;  i++)
    if (g_test_pattern[i].patch_type == DENSITY_PATCH)
      {
        set_grey_reflectance_pair(j, g_test_pattern[i].final_grey,
				  g_test_pattern[i].final_reflectance);
        j++;
      }
  number_grey_pairs = j;

  /* Sort grey_reflectance_pairs, i.e., sort on grey patch values, low to high */
  sort_grey_reflectance_pairs();
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
static int compare_grey_reflectance_pairs(const void *first, const void *second)
{
  /* 
     MOD 1 DAVE G. This code is guaranteed to work when the two
     values differ by less than 1.0 
  */

  if (((grey_reflectance_pair *)first)->grey < ((grey_reflectance_pair *)second)->grey)
    return -1;

  if (((grey_reflectance_pair *)first)->grey > ((grey_reflectance_pair *)second)->grey)
    return 1;

  return 0;

}


