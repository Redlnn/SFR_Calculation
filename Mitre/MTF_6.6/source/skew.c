/*********************************************************************
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
/*									     */
/* Revisions to skew.c							     */
/*									     */
/* V1.1 used SKIP_BLANK and GET_LINE macros in get_target_data() 	     */
/*      used EUCLID_DIS for sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))	     */
/*	cleaned up get_orientation()					     */
/*	added read_scan_line()						     */
/*	added do_white_patch()						     */
/*	made all source lines fit in 80 columns				     */
/*	now use read-in target dimensions instead of defines		     */
/*					brp - 7-MAR-94			     */
/*									     */
/* V1.3 replaced exp with pow and nint with NINT and M_PI with PI  7-JUN-94  */
/*	if skew angle is exactly 0.0; make it .00001 to avoid divide by zero */
/*									     */
/* V2.0 if the tag for the tif field PHOTOMETRICINTERPRETATION has been set  */
/*      to "white is zero", then we invert the greyscale values of the pixels*/
/*      						gtk 22-JUL-94        */
/*									     */
/* V2.1 the user now inputs the inner coordinates of the four patches at the */
/*      corners of the target rather than the outer coordinates.  From these */
/*      inputs, the outer coordinates are calculated based on the skew of the*/
/*      target.                                                              */
/*                                                      gtk 23-AUG-94        */
/*									     */
/* V2.3	Added additional message on error opening data file.		     */
/*	Moved free(buf) in get_orientation().				     */
/*	Added a close for the data file.				     */
/*							DJB 6/12/95	     */
/*									     */
/* V3.0 Indentation changes.						     */
/*							DJB 7/27/95	     */
/*	Conditionalize on USE_TIFF to not depend on tiff code.		     */
/*							DJB 12/28/95	     */
/*	Added out-of-bounds check for lower right outside corner.	     */
/*	Added tolerance to testing that outside corners are within image,    */
/*	 because the patch widths and heights aren't always exactly right.   */
/*	User-input inner corner coordinates are now float.  Use floating     */
/*	 point values for ppi and skew calculation, round to ints for	     */
/*	 everything else.						     */
/*							DJB 1/22/96	     */
/*	For DOS, open debug image file in BINARY mode.			     */
/*	Delete all orient.raw (g_orient_debug) code.			     */
/*							DJB 2/22/96	     */
/*									     */
/* V3.1	Removed copyright and added data rights notice.			     */
/*							DJB 6/25/96	     */
/*	In get_target_data, don't let read of patch_label read beyond end of */
/*	string in data file.						     */
/*							DJB 8/5/96	     */
/*	Added messages for bad file closes.				     */
/*							DJB 8/14/96	     */
/*									     */
/* V4.0 Added support for 1000 ppi.					     */
/*	Replaced average ppi with horizontal- or vertical-specfic ppi in     */
/*      Added test for image polarity                                        */
/*      get_target_data follows read-in density values and sets which patch  */
/*              is white and which patch is black                            */
/*      General cleanup to adhere to gcc's -ansi, -pedantic, and -Wall       */
/*		switches                                                     */
/*							WE 9/99		     */
/*									     */
/* V4.1 When printer option P is invoked, printer tgt data file mm and cy/mm */
/*      values are scaled by PRINTER_TARGET_BASE, allowing use of the same   */
/*      digital printer tgt and data file for MTF at any printer resolution. */
/*                                                      nbn, 9/00            */
/*									     */
/* V5.0 Several key changes/additions:                                       */
/*      1. Check for a common labelling error within the datafile.           */
/*         If present, info is automatically corrected at runtime            */
/*      2. When 'b' option is selected, the input corner points are used     */
/*         as starting points for an auto corner selection process.          */
/*      3. Orientation is discovered earlier and used in skew/ppi calculation*/
/*      4. External corners determined using PPI*Rot matrix mult, rather than*/
/*                                           Rot*PPI                         */
/*                                                      mal, 10/03           */
/*									     */
/* V5.5     Added capability to handle bar target images without density     */
/*          patches. The code automatically switches to bar target mode      */
/*	    when patch type 'b' is encountered in the target data file and   */
/*          automatically adds 4 virtual corner density patches to allow     */
/*	    processing to continue. get_target_data is called from mtf.c now.*/
/*          print_upright_info added to print description of 'upright'       */
/*          position.                                                        */
/*	    Comparison of target data file inputs (stated total width/height */
/*	    vs area covered by patches, and white_density vs black_density). */
/*          Problem printed if inconsistency found.                          */
/*                                                      mal, 02/04           */
/* V5.6     Added CTF spec lower bounds for bar targets                      */
/*	    Modified the target data file correction, so that it does an     */
/*	    extra test to determine relative patch positions on the page.    */
/*	    This leaves certain manufacturor target data files alone since   */
/*	    the will work OK.                                                */
/*                                                      mal, 06/04           */
/* V5.7.5   Changed MTF spec bounds to equation computation rather than      */
/*	    table lookup.  All spec bounds now based on effective freq.      */
/*	    Effective frequency set here rather than sine.c                  */
/*                                                      mal, 09/04           */
/* V5.8     Cleaned up code.  Sections no longer used were removed.          */
/*                                                      mal, 09/04           */
/* V5.8.1   Changed include info to allow CodeWarrior v9 compilation         */
/*          Moved NINT definition here to accommodate compilation.           */
/*          Moved TIFF manipulation definitions elsewhere for compilation.   */
/*          Allow MTF bounds test to occur when cos_skew_ang is slightly     */
/*          less than 1.                                                     */
/*                                                      mal, 03/05           */
/* V6.0     Added 'g' option handling for digital images.  This scales       */
/*          the effective frequency by the obj_to_img_scale factor.          */
/*                                                      mal, 06/05           */
/* V6.2     Added non-linear reflectance data for bartarget mode   	     */
/*                                                      MAL, 08/06           */
/* V6.3     Added PIV spec checks                                 	     */
/*          Print error message if target modulation is outside (0,1)        */
/*          Fixed print_upright_info bug.                                    */
/*                                                      MAL, 10/06           */
/* V6.3.01  Fixed minor range checking bug introduced in v6.3       	     */
/*                                                      MAL, 12/06           */
/* V6.4     Fixed 6.3 bug in AppF CTF 1000 eqn coefficients                  */
/*                                                      MAL, 11/07           */
/* V6.5     Added colors and additional outlines to debug images.	     */
/*                                                      AWU, 07/14           */
/* V6.6     Added capability to pre-set corner density patches (via 'c')     */
/*          when they are significantly inset from minimal rectangle         */
/* 		  Fixed bug in PPI computation.							  */
/*           Changed corner parameters to floating point.                    */
/*                                                      MAL, 04/16           */
/*****************************************************************************/


#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <fcntl.h>
#include <errno.h>

#if defined(MSDOS)
#include <io.h>
#else
#include <unistd.h>
#endif

typedef char TIFF;
#include "mtf.h"
#include "xterndcl.h"
extern long TIFFScanlineSize(TIFF *);

/* Local globals */
int density_order_correction;

/* 
   the next subroutine was included to adapt the existing code, 
   which used the following non-ansi C routine -  5/11/94  GTK  
   This always rounds 0.5 values up towards +infinity.
*/

int NINT(double x)
{
  if ((x-floor(x)) < 0.5) 
    return ((int) floor(x));
  else return ((int) ceil(x));
}


/* 
   Round() routine removed since NINT does same job, and rounds consistently
   positive, which is better for this code.  MAL-03/05
*/

/* V6.3  Coefficients for lower bound formulas for AppF(0) and PIV(1)
   specs for both 500 and 1000 ppi, MTF and CTF */
double lower_MTF_500[2][5] = {
  {1.01306, -1.062835e-1, -2.818546e-3,  1.092026e-3, -5.603957e-5},
  {1.02829, -1.67473e-1 ,  1.06255e-2 , -2.80874e-4,   0.0 }};
double lower_CTF_500[2][5] = {
  {1.02774, -7.99095e-2, 3.04105e-4,  0.0,         0.0},
  {1.00838, -8.05399e-2, -8.94631e-3, 1.43781e-3, -5.71711e-5}};
/* V6.4 Coeff lower_CTF_1000[0][2] and [0][3] corrected (bad typo during v6.3 conversion) */
double lower_CTF_1000[2][5] = {
  {1.01341, -5.73701e-2, 1.41666e-3, -1.85487e-5, 0.0},
  {0.0, 0.0, 0.0, 0.0, 0.0} /* Never used !! */};
/* V6.3  Coefficients for upper bounds AppF(0) and PIV(1) specs */
double upper[2] = { 1.05, 1.12 };

short initialize(char *data_filename, unsigned int white_patch, TIFF *tif)
{
  unsigned int i,j;
  int TgtType;
  double test_pattern_width, test_pattern_height;
  short rotation = -1;
  double cos_skew_avg;
  double *lower_coeff;
  int spec;
  double sum, power;

  /* Initializations. */
  g_image_array = NULL;		/* holds the input image */
  g_image_array_modified = NULL; /* copy of image for debug output */

  /* Compute to typical patch size in mm and inches for later use */
  test_pattern_width  = g_test_pattern[1].width_in_mm;
  test_pattern_height = g_test_pattern[1].height_in_mm;
  g_patch_height = test_pattern_height/MM_PER_INCH;
  g_patch_width  = test_pattern_width /MM_PER_INCH;

  if (g_refine_corners && !g_bartarget) {
    /* auto refine the corner positions */
    approx_ppi(g_scan_orient);
    refine_corner_coords(tif);
  }

  /* compute the skew angles and the pixels per inch */
  compute_skew_angles(g_scan_orient);

  /* determine the actual orientation and white patch position
     on the input image.  This allows linkage of input corners to 
     test_pattern patch info. V5.0 
     use a different method of determining orientation when 
     processing bar targets V5.5*/

  if (!g_bartarget)
    rotation = get_orientation(tif, white_patch); 
  else
    switch(g_scan_orient) {
    case LAND:  
      if (g_upright == 0 || g_upright == 90) rotation = RIGHT;
      else  rotation = LEFT;
      break;
    case PORT:  
      if (g_upright == 0 || g_upright == -90) rotation = TOP;
      else rotation = BOTTOM;
      break;
    }

  /* Compute skew and ppi using precise patch dimension info. V5.0   */
  /* Eliminates minor perturbations caused by patch size differences */
  compute_correct_skews(rotation); 

  /* Since v2.1, the user is inputting the coordinates of the inner corners of 
     the four corner patches of the target.  Since the code was originally 
     written to use the outer corners of these four patches, we must calculate 
     the locations of these outer corners based on inner coordinates, and the 
     physical characteristics of the target - gtk 8/19/94 */

  internal_to_external_corners(test_pattern_width, test_pattern_height,
			       rotation); 

  /* These are kept as floats to avoid using precision for later computations */
  g_test_pattern_image_height = g_test_pattern_yll-g_test_pattern_yul;
  g_test_pattern_image_width  = g_test_pattern_xur-g_test_pattern_xul;

  g_orig_test_pattern_xul = g_test_pattern_xul;
  g_orig_test_pattern_xur = g_test_pattern_xur;
  g_orig_test_pattern_xll = g_test_pattern_xll;
  g_orig_test_pattern_xlr = g_test_pattern_xlr;

  g_orig_test_pattern_yul = g_test_pattern_yul;
  g_orig_test_pattern_yur = g_test_pattern_yur;
  g_orig_test_pattern_yll = g_test_pattern_yll;
  g_orig_test_pattern_ylr = g_test_pattern_ylr;


/* If sine tgt model is M6, M13, or M15 (TgtType=1) then leftedge of leftmost 
   graypatch is aligned with leftedge of leftmost sine pattern (in program 
   standard orientation) so can run test to see if Entire sine tgt is in image.
   For sine tgt model M5, M14, or M19 (TgtType=-1) the leftend sine pattern 
   extends beyond the leftend graypatch (so xul_in_mm has negative value for 
   this sine pattern in TgtDataFile) so test must be skipped to allow program 
   to complete.  
*/  
   
  TgtType = 1;							 
  for (i=0;i<g_number_of_patches;i++) {
    if (g_test_pattern[i].patch_type == SINE_PATCH &&
	g_test_pattern[i].xul_in_mm < 0.0)  
      TgtType = -1;      
  }
    
  if (TgtType == 1) {

  /* Test to see if calculated outer corners are within image.  We add in a 
     little tolerance because the patch widths and heights aren't always
     exactly right. */
    if ( (g_orig_test_pattern_xul < -WITHIN_IMAGE_TOLERANCE)
	 || (g_orig_test_pattern_yul < -WITHIN_IMAGE_TOLERANCE)
	 || (g_orig_test_pattern_xur >= (g_scan_cols + WITHIN_IMAGE_TOLERANCE))
	 || (g_orig_test_pattern_yur < -WITHIN_IMAGE_TOLERANCE)
	 || (g_orig_test_pattern_xll < -WITHIN_IMAGE_TOLERANCE)
	 || (g_orig_test_pattern_yll >= (g_scan_rows + WITHIN_IMAGE_TOLERANCE))
	 || (g_orig_test_pattern_xlr >= (g_scan_cols + WITHIN_IMAGE_TOLERANCE))
	 || (g_orig_test_pattern_ylr >= (g_scan_rows + WITHIN_IMAGE_TOLERANCE)))
      {
	MTFPRINT2("\nImage input file name = \t%s\n",image_filename)           
        MTFPRINT2("Data input file name = \t\t%s\n",data_filename)           
        MTFPRINT("\nInput coordinates of INNER CORNERS of sine image:\n")      
        MTFPRINT3("UpperLeftCol =\t\t%g\t   UpperLeftRow =\t%g\n",       
		  g_input_test_pattern_xul,g_input_test_pattern_yul)        
        MTFPRINT3("UpperRightCol =\t\t%g\t  UpperRightRow =\t%g\n", 
		  g_input_test_pattern_xur,g_input_test_pattern_yur)     
        MTFPRINT3("LowerLeftCol =\t\t%g\t   LowerLeftRow =\t%g\n",   
		  g_input_test_pattern_xll,g_input_test_pattern_yll)    
      
        MTFPRINT("\nCalculated coordinates of OUTER CORNERS of sine image:\n");
        MTFPRINT3("\nUpperLeftCol =\t\t%.1f\t   UpperLeftRow =\t%.1f\n",
		  g_orig_test_pattern_xul,g_orig_test_pattern_yul);
        MTFPRINT3("UpperRightCol =\t\t%.1f\t   UpperRightRow =\t%.1f\n",
		  g_orig_test_pattern_xur,g_orig_test_pattern_yur);
        MTFPRINT3("LowerLeftCol =\t\t%.1f\t   LowerLeftRow =\t%.1f\n",
		  g_orig_test_pattern_xll,g_orig_test_pattern_yll);
        MTFPRINT3("LowerRightCol =\t\t%.1f\t   LowerRightRow =\t%.1f\n",
		  g_orig_test_pattern_xlr,g_orig_test_pattern_ylr);

        MTFPRINT("\nQUITTING! Complete Sine Target not imaged or Wrong Orientation.\n");                                                                          
	MTFPRINT2("         End of Report for %s\n",image_filename);      
	MTFPRINT3("%s%s\n","****************************************",    
		  "****************************************\n\n");          
	mtfexit(-1);
      }
  }


#ifdef DEBUG_COORDS
  fprintf(stderr,"\nCalculated coordinates of outer corners of target:\n");
  fprintf(stderr,"\nUpper Left Col =\t\t%.1f\tUpper Left Row =\t%.1f\n",
	  g_orig_test_pattern_xul,g_orig_test_pattern_yul);
  fprintf(stderr,"Upper Right Col =\t\t%.1f\tUpper Right Row =\t%.1f\n",
	  g_orig_test_pattern_xur,g_orig_test_pattern_yur);
  fprintf(stderr,"Lower Left Col =\t\t%.1f\tLower Left Row =\t%.1f\n",
	  g_orig_test_pattern_xll,g_orig_test_pattern_yll);
  fprintf(stderr,"Computed Lower Right Col =\t%.1f\tLower Right Row =\t%.1f\n",
	  g_orig_test_pattern_xlr,g_orig_test_pattern_ylr);

#endif

  /* The values below are used by the read_in_image routine.  They define
     what area of the image to read in and store in a normalized orientation.
     Right now these values define the full image, which is usually larger
     than the test pattern.  This is safe to do since the image will probably
     be slightly skewed and this way we won't crop away any of the test
     pattern.  These can be refined to use a tighter bounding box, based on
     how skewed the image is.  Although these are used in the read_in_image
     routine, they are defined here to keep all the skew related stuff
     grouped together. */

  g_input_column_start = 0;	/* start of valid data in input scan line */
  g_input_last_column = g_scan_cols-1; /* end of valid data */
  g_input_row_start = 0;	/* row start of valid data */
  g_input_last_row = g_scan_rows-1; /* end of valid data */
  /* The next four items depend on the previous four and are just a shorter
     way to use the right hand side of the assignment. mo stands for 
     minus one */
  g_input_width_mo = g_input_last_column - g_input_column_start;
  g_input_width = g_input_width_mo + 1;
  g_input_height_mo = g_input_last_row - g_input_row_start;
  g_input_height = g_input_height_mo + 1;

  /* Upper and lower allowable MTF's for each patch - vary with eff. cy/mm */
  cos_skew_avg = cos(g_skew_avg);

  /* V6.3 Default spec computation is AppF.  Reset if PIV was requested. */
  spec = 0;
  if (g_compare_userset == PIV) spec = 1;
  if (g_bartarget) {
    if (g_target_res == PPI_BASE)
      lower_coeff = lower_CTF_500[spec];
    else
      lower_coeff = lower_CTF_1000[spec];
  }
  else {
    if (g_target_res == PPI_BASE)
      lower_coeff = lower_MTF_500[spec];
    else  /* Lower_coeff is not used. Direct computation used instead */
      lower_coeff = NULL;
  }

  /* Diagnostic messages for bounds formulas v6.3
  if (g_compare) {
    MTFPRINT2("\nUpper Bound Spec: %f", upper[spec])
    if (g_bartarget || (g_target_res == PPI_BASE)) {
      MTFPRINT("\nLower Bound Spec Equation: C0 + C1*f + C2*f^2 + C3*f^3 + C4*f^4\n")
      MTFPRINT6("Coefficients: %g %g %g %g %g\n\n", lower_coeff[0],lower_coeff[1],lower_coeff[2],lower_coeff[3],lower_coeff[4])
    } else {
      MTFPRINT("\nLower Bound Spec Equation: C0 * exp(C1*f)\n")
      MTFPRINT3("Coefficients: %g %g\n\n", MTF_A_1000, MTF_B_1000)
    }
  }
  */

  for (i=0;i<g_number_of_patches;i++) {
    if (g_test_pattern[i].patch_type == SINE_PATCH) {

      double f;

      /* V5.7.5 use effective frequency instead of target frequency */
      f = g_obj_img_scale * g_test_pattern[i].cycles_per_mm*cos_skew_avg;
      g_test_pattern[i].effective_frequency = f;
      
      /* V6.3 compute the formula */
      if (!g_bartarget && !(g_target_res == PPI_BASE))
	  g_test_pattern[i].mtfll = MTF_A_1000*exp(MTF_B_1000*f);
      else { /* A + B*f + C*f*f + D*f*f*f + E*f*f*f*f */
	power = 1.0;
	sum = 0.0;
	for(j=0; j<5; j++) {
	  if (lower_coeff[j] == 0.0) break;
	  sum += lower_coeff[j]*power;
	  power *= f;
	}
	g_test_pattern[i].mtfll = sum;
      }	  
      g_test_pattern[i].mtful = (double)upper[spec];
      /* V5.8.1 set lower bound for eff freq f, to accommodate 
	 effects of cos_skew_ang < 1 */
      if (f < .9998 || f > 20 || ((g_target_res == PPI_BASE) && f>10) ) {
	g_test_pattern[i].mtfll = (double)MTFLL_SMALL;
	g_test_pattern[i].mtful = (double)MTFUL_BIG;
      }
    }
  }

  return(rotation);

}

/******************************************************************************/
/*									      */
/******************************************************************************/

/* This is done on standardized image (horizontal, white patch at top) - the
   internal storage format of the image. We compute the corner coords of 
   a box 10% within each patch.  This internal area is used for measurements
   instead of the entire patch; a safety margin. */

void compute_patch_boundaries(void)
{
  short j;
  unsigned int i;
  double row,col;
  double sinsa,cossa;
  double tolerance;


  sinsa = sin(g_skew_avg);
  cossa = cos(g_skew_avg);

  for (i=0;i<g_number_of_patches;i++) {
    /* xul_in_mm, width_in_mm, etc. are read in from the patch data file for 
       each patch.  We locate a box 10% within the patch.  his is an imaginary 
       box on the test pattern itself; not on the scanned-in image.  The 
       imaginary box is then translated into row,col coords on the actual, 
       scanned image. */

    tolerance = (double)BORDER*g_test_pattern[i].width_in_mm;
    if (tolerance*g_ppi_horz < 3) tolerance = 3.0/g_ppi_horz;
    g_test_pattern[i].x[0] =
      g_test_pattern[i].x[2] =
	g_test_pattern[i].xul_in_mm + tolerance;
    g_test_pattern[i].x[1] =
      g_test_pattern[i].x[3] =
	g_test_pattern[i].xul_in_mm +
	  g_test_pattern[i].width_in_mm - tolerance;


    tolerance = (double)BORDER*g_test_pattern[i].height_in_mm;
    if (tolerance*g_ppi_vert < 3) tolerance = 3.0/g_ppi_vert;
    g_test_pattern[i].y[0] =
      g_test_pattern[i].y[1] =
	g_test_pattern[i].yul_in_mm + tolerance;
    g_test_pattern[i].y[2] =
      g_test_pattern[i].y[3] =
	g_test_pattern[i].yul_in_mm +
	  g_test_pattern[i].height_in_mm - tolerance;

    /* taking ppi and skew into account locate the row,col coords of the 10% 
       box on the actual image */
    for(j=0;j<4;j++) {
      row = (g_test_pattern[i].y[j]/MM_PER_INCH)*g_ppi_vert;
      col = (g_test_pattern[i].x[j]/MM_PER_INCH)*g_ppi_horz;

      g_test_pattern[i].c[j] =
	NINT(col*cossa - row*sinsa + g_test_pattern_xul);
      g_test_pattern[i].r[j] =
	NINT(col*sinsa + row*cossa + g_test_pattern_yul);

	}
	if(g_debug){
	/*	draws boundaries on patches (all of them)*/
		draw_bound(g_test_pattern[i].c[0], g_test_pattern[i].r[0], g_test_pattern[i].c[1], g_test_pattern[i].r[1],COLORB,0);
		draw_bound(g_test_pattern[i].c[1], g_test_pattern[i].r[1], g_test_pattern[i].c[3], g_test_pattern[i].r[3],COLORB,0);
		draw_bound(g_test_pattern[i].c[3], g_test_pattern[i].r[3], g_test_pattern[i].c[2], g_test_pattern[i].r[2],COLORB,0);
		draw_bound(g_test_pattern[i].c[2], g_test_pattern[i].r[2], g_test_pattern[i].c[0], g_test_pattern[i].r[0],COLORB,0);
    }
  }
}

/******************************************************************************/
/*									      */
/******************************************************************************/

/* Approximate skew angles and ppi scaling without knowing white_patch
   location */
void compute_skew_angles(short orient)
{

  /* We can have positive or negative skew.  Positive skew is when the upper
     right row of the image is greater than the upper left row (system has
     the origin at the upper left corner of the "screen"). */

  g_skew_horz = asin((double)(g_input_test_pattern_yur - g_input_test_pattern_yul)/
		     EUCLID_DIS(g_input_test_pattern_xul,g_input_test_pattern_yul,
				g_input_test_pattern_xur,g_input_test_pattern_yur));

  g_skew_vert = asin((double)(g_input_test_pattern_xul - g_input_test_pattern_xll)/
		     EUCLID_DIS(g_input_test_pattern_xul,g_input_test_pattern_yul,
				g_input_test_pattern_xll,g_input_test_pattern_yll));


  g_skew_avg = (g_skew_horz+g_skew_vert)/2.0;

  /* avoid divide by zero */
  if (g_skew_avg == (double)0.0) g_skew_avg = (double)0.000001; 

  approx_ppi(orient);
  
}


/******************************************************************************/
/*									      */
/******************************************************************************/

/* Approximate the ppi from the input corner points.  Must be called
   before automatic corner refinement. */

void approx_ppi(short orient)
{

  if (orient == LAND) {
    g_ppi_horz = EUCLID_DIS(g_input_test_pattern_xul,g_input_test_pattern_yul,
			    g_input_test_pattern_xur,g_input_test_pattern_yur)/(g_target_width-2.0*g_patch_width);
    g_ppi_vert = EUCLID_DIS(g_input_test_pattern_xul,g_input_test_pattern_yul,
			    g_input_test_pattern_xll,g_input_test_pattern_yll)/(g_target_height-2.0*g_patch_height);
  }
  else {
    g_ppi_vert = EUCLID_DIS(g_input_test_pattern_xul,g_input_test_pattern_yul,
			    g_input_test_pattern_xur,g_input_test_pattern_yur)/(g_target_height-2.0*g_patch_height);
    g_ppi_horz = EUCLID_DIS(g_input_test_pattern_xul,g_input_test_pattern_yul,
			    g_input_test_pattern_xll,g_input_test_pattern_yll)/(g_target_width-2.0*g_patch_width);
  }

  if (((int)g_ppi_horz > (int)PPI_THRESHOLD) &&
      ((int)g_ppi_vert > (int)PPI_THRESHOLD))
      g_target_res = PPI_HR_BASE;
  else
      g_target_res = PPI_BASE;
}


/******************************************************************************/
/*									      */
/******************************************************************************/

/* Sum up the pixels in a "trial" white patch - this is called from
	get_orientation() */

void do_white_patch(unsigned char *buf,TIFF *tif,double *sum,
		    unsigned int *count,short rstart,short rend,short cstart,short cend)
{
  short j,k;
  unsigned short *sbuf;

  for(j=rstart;j<rend;j++) {
    read_scan_line(tif,buf,j);
    
    switch (g_bytes_per_pixel) {
    case 1:
      for(k=cstart;k<=cend;k++) {
	(*sum)+=(double)buf[k];
	(*count)++;
      }
      break;
    case 2:
      sbuf = (unsigned short *)buf;
      for(k=cstart;k<=cend;k++) {
	(*sum)+=(double)sbuf[k];
	(*count)++;
      }
      break; 
    default:
      fprintf(stderr,"\nCannot handle more than 2 bytes per pixel!!\n\n");
      mtfexit(-1);
    }
  }
}


/******************************************************************************/
/*									      */
/******************************************************************************/
/* test for image polarity --WE 9/99 V4.0 */
/* get_orientation() is called *before* we have computed the mean grey levels
   so test for image polarity must be done here instead */
void test_for_image_polarity(unsigned int white_patch)
{
  unsigned int i;
  unsigned int whitest_patch;
  double max_gray;

  max_gray = 0.0;
  whitest_patch = 0;
  for (i=0;i<g_number_of_patches;i++) {
    if (g_test_pattern[i].grey_level > max_gray) {
      max_gray = g_test_pattern[i].grey_level;
      whitest_patch = i;
    }
  }

  if (density_order_correction == 1) {
    MTFPRINT("\nATTN: Density ordering in datafile was automatically corrected\n");
  }
  
   
   if (white_patch != whitest_patch) {                                                                                                       
     MTFPRINT("\nPOSSIBLE PROBLEM WITH IMAGE PIXEL VALUES:");
     MTFPRINT("\n  if TgtReflect vs MeasImgGray sequence is correct, probably no problem."); 
     MTFPRINT("\n  otherwise,if 8bpp, try reversing image polarity (Menu option f).");
     MTFPRINT("\n            if 16bpp, image may be signed integer, which prg cant handle.\n");
   }  
  
}

/******************************************************************************/
/*									      */
/******************************************************************************/
/* test for bar target image polarity -- MAL  01/04 V5.5 */
void test_for_bar_target_polarity(unsigned int lowf_patch)
{
  int i, j;
  double effective_freq, ppp, ppqp;
  double sinsa, cossa;
  unsigned char *img_ptr;
  unsigned short *simg_ptr;
  double test_row, test_col;
  unsigned int row, col;
  int qs, es, ci;
  int sum[2];

  /* Choose a row that is within the lowest frequency bar target. */  
  test_row = ((g_test_pattern[lowf_patch].yul_in_mm
              +(double)BORDER*g_test_pattern[lowf_patch].height_in_mm)
	      /MM_PER_INCH)*g_ppi_vert;
  /* Find the column at the leftmost edge of the lowest freq bar patch */
  test_col = (g_test_pattern[lowf_patch].xul_in_mm/MM_PER_INCH)*g_ppi_horz;

  sinsa = sin(g_skew_avg);
  cossa = cos(g_skew_avg);

  /* Find actual positions on the normalized orientation skew image. */
  col =	NINT(test_col*cossa - test_row*sinsa + g_test_pattern_xul);
  row =	NINT(test_col*sinsa + test_row*cossa + g_test_pattern_yul);

  /* pixels per period */
  effective_freq = g_test_pattern[lowf_patch].effective_frequency;
  ppp = g_ppi_horz / (effective_freq * MM_PER_INCH);
  /* pixels/quarter-period */
  ppqp = ppp / 4.0;

  qs = NINT(ppqp);
  if (qs == 0) qs = 1;
  es = qs/2;


  img_ptr = g_image_array + row*g_total_image_width;
  simg_ptr = (unsigned short *)g_image_array + row*g_total_image_width;

  /* Look area around quarter period to either side of leftmost
     edge (col) of the low frequency pattern.  This should include
     the background in sum[0] and the first bar in sum[1].
     If the polarity is correct the background (space) will be brighter
     than the bar. */
  ci = col - qs;
  for (i=0; i<2; i++,ci+=2*qs)
  {
     sum[i] = 0;
     for (j=ci-es;j<=ci+es;j++) 
     {
       if (g_bytes_per_pixel == 1) 
           sum[i] += (int)(img_ptr[j]);
      else
           sum[i] += (int)(simg_ptr[j]);
     }
  }

  if (sum[0] <= sum[1]) {                                                                                                       
     MTFPRINT("\nPOSSIBLE IMAGE POLARITY PROBLEM:");
     MTFPRINT("\n  Expected polarity is black bars and white spaces, but program is "); 
     MTFPRINT("\n  reading opposite on snippet of lowest frequency bar pattern, i.e."); 
     MTFPRINT3("\n  bar gray = %d and space gray = %d.", sum[1]/(2*es+1), sum[0]/(2*es+1))  
     MTFPRINT("\n  Not a problem if this bar pattern's location in TgtDataFile starts") 
     MTFPRINT("\n  with white space adjacent to leftmost bar, instead of leftmost bar.") 
     MTFPRINT("\n  Otherwise,if 8bpp, reverse image polarity (Menu option f).");
     MTFPRINT("\n            if 16bpp, image may be signed integer, which program cant handle.\n");
   }  
  
}

void read_in_nonlinear(char *line, FILE * fd) {
  int i;
  int tot;
  long pos;
  double grey, reflect;

  sscanf(line, "%lf %lf", &reflect, &grey);

  pos = ftell(fd); /* Remember where we are in the file */

  tot = 1;
  /* Count the number of non-blank/non-comment lines */
    while(1 /*CONSTANTCONDITION*/) {
      if (fgets(line, 81, fd) == NULL) {
	break;
      }
      SKIP_BLANK
      /* Any extra lines are for nonlinear system I/O */
      tot++;
    }

  get_grey_reflectance_pairs(tot);

  set_grey_reflectance_pair(0, grey, reflect);

  fseek(fd, pos, SEEK_SET); /* Go back to marked place in the file */

  i = 1;
  while ( fscanf(fd, "%lf %lf", &reflect, &grey) == 2 && i < tot ) {
    set_grey_reflectance_pair(i, grey, reflect);
    i++;
  }

  sort_grey_reflectance_pairs();
}


/******************************************************************************/
/*									      */
/******************************************************************************/

/* Read in the patch data file - store patch dimensions into g_test_pattern 
   Return index of the white_patch, black_patch, and the lowest frequency
   sine patch (v5.5) for use in later orientation and polarity tests. 
   Return flag indicating whether file could be opened. */

int get_target_data(char *data_filename, unsigned int *white_patchp, unsigned int *black_patchp, unsigned int *lowf_patchp)
{
  FILE *fd;
  char line[82];
  unsigned short i,entries;
  unsigned int j;
  float val[7],tmp;
  unsigned char patch_type;
  float d;
  unsigned char patch_label[82];
  double min_x, max_x, min_y, max_y;
  int extent_error;
  int bounded_modulation;
  int num_density_patches = 0;
  int rel_white_patchp = 0, rel_black_patchp = 0;
  float white_density = 99.99f, black_density = 0.0f, lowest_freq = FLT_MAX;
  *white_patchp = 0, *black_patchp = 0;

  /* Initialize check for patch extent v5.5 */
  max_x = max_y = 0.0;

  /* Make g_bartarget undefined, until patch data is read */
  g_bartarget = 2;

  /* Start out assuming the modulation is bounded  (0,1) */
  bounded_modulation = 1;

  memset((char *)patch_label, 0, 82);

  fd = fopen(data_filename,"r");
  if (fd == NULL) {
    fprintf(stderr,"\nCould not open %s\n",data_filename);
    fprintf(stderr,"%s\n\n", strerror(errno));
    return(0);
  }

  /* we don't have any keywords to look for; the patch data file is based
     on the position of the entries - we look for them in order. It may be
     better to have a system based on keywords.  If the maker of the patch
     data file gets anything out of order ... */
	
  while (1 /*CONSTANTCONDITION*/) {	/* loop until we get the first entry (the serial #) */
    GET_LINE
    /* skip over blank lines and comment lines */
    SKIP_BLANK
    if (sscanf(line,"%s", g_serial_number) == 1) break;
    else {
      TgtFileErrMsg;
      mtfexit(-1);
    }
  }

  /* number of patches in file (on target) */
  while (1 /*CONSTANTCONDITION*/) {
    GET_LINE
    SKIP_BLANK
    if (sscanf(line,"%u", &g_number_of_patches) == 1) break;
    else {
      TgtFileErrMsg;
      mtfexit(-1);
    }
  }

  /* target width */
  while (1 /*CONSTANTCONDITION*/) {
    GET_LINE
    SKIP_BLANK
    if (sscanf(line,"%f", &tmp) == 1) {
      g_target_width_in_mm = tmp;
      if(g_printer) {
      	g_target_width_in_mm = g_target_width_in_mm * PRINTER_TARGET_BASE/g_printer_ppi;
      }
      g_target_width = g_target_width_in_mm / MM_PER_INCH;
      break;
    }
    else {
      TgtFileErrMsg;
      mtfexit(-1);
    }
  }

  /* target height */
  while (1 /*CONSTANTCONDITION*/) {
    GET_LINE
    SKIP_BLANK
    if (sscanf(line,"%f", &tmp) == 1) {
      g_target_height_in_mm = tmp;
      if(g_printer) {
      	g_target_height_in_mm = g_target_height_in_mm * PRINTER_TARGET_BASE/g_printer_ppi;
      }
      g_target_height = g_target_height_in_mm / MM_PER_INCH;
      break;
    }
    else {
      TgtFileErrMsg;
      mtfexit(-1);
    }
  }

  min_x = g_target_width_in_mm;
  min_y = g_target_height_in_mm;

  /* get enough memory to cover the data for all the patches */
  g_test_pattern
    = (struct test_pattern *) calloc((int)g_number_of_patches+4,(int)(sizeof(struct test_pattern)));
  if ( g_test_pattern == NULL ) {
    fprintf(stderr,"\nCan't allocate memory for patch data!!\n\n" );
    mtfexit(-1);
  }

  j=0;
  g_corners_read = 0;
  while(1 /*CONSTANTCONDITION*/) {
    while(1 /*CONSTANTCONDITION*/) {
      GET_LINE
      SKIP_BLANK
      /* the patch label */
      if (sscanf(line,"%s",patch_label) != 1) {
	TgtFileErrMsg;
	mtfexit(-1);
      }
      else {
	for (i=0;i<8;i++){
	  g_test_pattern[j].patch_label[i] = patch_label[i];
	  g_test_pattern[j].new_patch_label[i] = patch_label[i];
	}
	g_test_pattern[j].patch_label[9] = '\0';
	g_test_pattern[j].new_patch_label[9] = '\0';
	break;
      }
    }
    while(1 /*CONSTANTCONDITION*/) {
      GET_LINE
      SKIP_BLANK
      /* patch type density or sine ? */
      if (sscanf(line,"%c",&patch_type) != 1) {
	TgtFileErrMsg;
	mtfexit(-1);
      }
      if (patch_type == 'S' || patch_type == 's') {
        if (g_bartarget == 1) {
	   TgtFileBarMsg;
	   mtfexit(-1);
	}
        g_bartarget = 0;
	g_test_pattern[j].patch_type = SINE_PATCH;
	entries = 6;
	break;
      }
      else if (patch_type == 'd' || patch_type == 'D'  ||
      	       patch_type == 'c' || patch_type == 'C') {
        if (g_bartarget == 1) {
	   TgtFileBarMsg;
	   mtfexit(-1);
	}
        g_bartarget = 0;
	g_test_pattern[j].patch_type = DENSITY_PATCH;
	entries = 5;
	g_test_pattern[j].corner_flag = 0;
	if (patch_type == 'c' || patch_type == 'C') {
		g_corners_read += 1;      
		if (g_corners_read > 4) {
			TgtFileCnrMsg;
			mtfexit(-1);
		}	
		g_test_pattern[j].corner_flag = 1;
	}
	break;
      }
       else if (patch_type == 'B' || patch_type == 'b') {
        if (g_bartarget == 0) {
	   TgtFileBarMsg;
	   mtfexit(-1);
	}
        g_bartarget = 1;
	g_test_pattern[j].patch_type = SINE_PATCH;
	entries = 6;
	break;
      }
      else {
	TgtFileErrMsg;
	mtfexit(-1);
      }
    }

    i=0;
    while(i<entries) {
      GET_LINE
      SKIP_BLANK
      if (sscanf(line,"%f",&val[i]) != 1) {
	TgtFileErrMsg;
	mtfexit(-1);
      }
      else i++;
    }
      
    if(g_printer) {
    	val[0] = val[0] * PRINTER_TARGET_BASE/g_printer_ppi; 
	val[1] = val[1] * PRINTER_TARGET_BASE/g_printer_ppi;
	val[2] = val[2] * PRINTER_TARGET_BASE/g_printer_ppi;
	val[3] = val[3] * PRINTER_TARGET_BASE/g_printer_ppi;
    }

    g_test_pattern[j].xul_in_mm = val[0];
    g_test_pattern[j].yul_in_mm = val[1];
    g_test_pattern[j].width_in_mm = val[2];
    g_test_pattern[j].height_in_mm = val[3];

    /* Determine maximum patch coverage in width and height. 
       UL corner from minimum patch ul corners.  LR corner
       from maximum of patch ul corner plus patch width/height. */

    if (entries != 6 || val[0] >= 0.0) {
      /* Ignore sine patches with negative x location, but
	 otherwise find minimum x location. */
      if (val[0] < min_x) min_x = val[0];
    }
    if (val[1] < min_y) min_y = val[1];
    if (val[0] + val[2] > max_x) max_x = val[0] + val[2];
    if (val[1] + val[3] > max_y) max_y = val[1] + val[3];

    /* For debugging data file patch info 
    MTFPRINT8("Patch %d at (%g,%g) -> (%g,%g) = %g x %g\n",
	      j, val[0],val[1], val[0]+val[2], val[1]+val[3], val[2],val[3])
    */

    
    if (entries == 5) {
      num_density_patches++;
      g_test_pattern[j].patch_density = (double)val[4];
      d = (float)(g_test_pattern[j].patch_density);
      if (d < white_density) {
        white_density = d;
        *white_patchp = j;
        rel_white_patchp = num_density_patches;
      }
      if (d > black_density) {
        black_density = d;
        *black_patchp = j;
        rel_black_patchp = num_density_patches;
      }
      /* reflectance = 10**-density */
      g_test_pattern[j].reflectance = 1.0/pow((double)10.0,d);
      g_test_pattern[j].modulation = 0.0;
      g_test_pattern[j].cycles_per_mm = 0.0;
    }
    else {
      g_test_pattern[j].patch_density = 0.0;
      g_test_pattern[j].modulation = val[4];
      g_test_pattern[j].cycles_per_mm = val[5];
      /* Find lowest frequency sine patch for bar target processing v5.5 */
      if (val[5] < lowest_freq) {
         lowest_freq = val[5];
	 *lowf_patchp = j;
      }
      if (val[4] > 1 || val[4] < 0)
	bounded_modulation = 0;  /* V5.6 Modulation goes out of bounds? */
       
      /* If calculating printer MTF, adjust patch frequencies for resolution
	 of printer, which is tied to printer target resolution */
      if (g_printer)
	g_test_pattern[j].cycles_per_mm
	  = g_test_pattern[j].cycles_per_mm * g_printer_ppi / PRINTER_TARGET_BASE;
    }
    if (++j >= g_number_of_patches) break;
  }
  
  /* No patches read at all.  Revert to sine/density mode */
  if (g_bartarget == 2) g_bartarget = 0;

  if (g_bartarget)
  {
     g_number_of_patches += 4;
     g_corners_read = 4;

     /* Create four virtual density patches of zero area, one
        at each corner. The only info of importance for later
	processing is the location, size, and type of these patches. V5.5 */
     g_test_pattern[j].xul_in_mm = 0.0;
     g_test_pattern[j].yul_in_mm = 0.0;
     g_test_pattern[j+1].xul_in_mm = g_target_width_in_mm;
     g_test_pattern[j+1].yul_in_mm = 0.0;
     g_test_pattern[j+2].xul_in_mm = g_target_width_in_mm;
     g_test_pattern[j+2].yul_in_mm = g_target_height_in_mm;
     g_test_pattern[j+3].xul_in_mm = 0.0;
     g_test_pattern[j+3].yul_in_mm = g_target_height_in_mm;

     for (i=0; i<4; i++) {
          g_test_pattern[j].patch_type = DENSITY_PATCH;
          sprintf(g_test_pattern[j].patch_label,"FAKE%1d",i);
	  g_test_pattern[j].patch_density = 0.0;
	  g_test_pattern[j].reflectance = 0.0;
	  g_test_pattern[j].modulation = 0.0;
	  g_test_pattern[j].cycles_per_mm = 0.0;
	  g_test_pattern[j].width_in_mm = 0.0;
	  g_test_pattern[j].height_in_mm = 0.0;
	  g_test_pattern[j].corner_flag = 1;
          num_density_patches++;
	  j++;
     }

    while(1 /*CONSTANTCONDITION*/) {
      if (fgets(line, 81, fd) == NULL) {
	g_nonlinear = 0;
	break;
      }
      SKIP_BLANK
      /* Any extra lines are for nonlinear system I/O */
      g_nonlinear = 1;
      break;
    }

    if(g_nonlinear) {
      read_in_nonlinear(line, fd);
    }
  }

  /* Check that if corners are marked, then there are exactly 4 of them */
  if (g_corners_read > 0 && g_corners_read < 4) {
	TgtFileCnrMsg;
	mtfexit(-1);
  }	

  /* Check that patches do not extend beyond target area v5.5
     Print problem report if this is noticed.
   */
  extent_error = 0;
  if (max_x-min_x-g_target_width_in_mm > ZERO_BOUND) extent_error = 1;
  if (max_y-min_y-g_target_height_in_mm > ZERO_BOUND) extent_error = 1;

  if (extent_error) {
    char problem_string[82];

    sprintf(problem_string,
	    "Target data file:\n");
    put_problem(problem_string,0);
    sprintf(problem_string, "  Measurements from patches (%g = width, %g = height) exceeds\n",
	    max_x-min_x, max_y-min_y);
    put_problem(problem_string,0);
    sprintf(problem_string, "  stated target width x target height = %g x %g\n\n",
	  g_target_width_in_mm, g_target_height_in_mm);
    put_problem(problem_string,0);
  
  }
  
  /* Check to see if all density patches have 'same' density.
     If so print error message and stop. v5.5 */
  if (!g_bartarget && (black_density-white_density<ZERO_BOUND)) {
    MTFPRINT("\n\nSTOP - all density patches in TgtDataFile have same density value\n");
    MTFPRINT2("\tTgtDataFile = %s\n\n", data_filename);
    mtfexit(-1);
  }

  /* Check to see if a common data file error is present.  If so, correct it. 
     V5.0  mal 10/03 */
  density_order_correction = 0;
  
  if (!g_bartarget && 
  num_density_patches==14 && rel_white_patchp==2 && rel_black_patchp==13) {
    float ax, ay, bx, by;
    /* V5.6 Add extra test to figure out relative positions of first patch
       vs white & black patch on the page.  If the cross product has positive 
       sign, then the file is improperly constructed.  Otherwise it's probably
       OK. */
    ax = g_test_pattern[rel_white_patchp].xul_in_mm-g_test_pattern[0].xul_in_mm;
    ay = g_test_pattern[rel_white_patchp].yul_in_mm-g_test_pattern[0].yul_in_mm;
    bx = g_test_pattern[rel_black_patchp].xul_in_mm-g_test_pattern[0].xul_in_mm;
    by = g_test_pattern[rel_black_patchp].yul_in_mm-g_test_pattern[0].yul_in_mm;

    if( ax*by > ay*bx) {
    /* Improperly constructed data file.  Attempted to fit sinepatterns density
       info into MITRE labels.  Rearrange densities to correct MITRE order */
    double *density,*reflect;
    int dpatch_num;

    MTFPRINT("\nNOTE - Density patches are not in expected order: Automatically reordered.\n");
    density=(double *) calloc(num_density_patches,(int)(sizeof(double)));
    reflect=(double *) calloc(num_density_patches,(int)(sizeof(double)));
    
    if ( density == NULL || reflect == NULL) {
      fprintf(stderr,"\nCan't allocate memory for temporary density data!!\n\n" );
      mtfexit(-1);
    }
  
    dpatch_num=0;
    for (i=0;i<g_number_of_patches;i++) {
      if (g_test_pattern[i].patch_type == DENSITY_PATCH) {
	density[dpatch_num]=g_test_pattern[i].patch_density;
	reflect[dpatch_num]=g_test_pattern[i].reflectance;
	dpatch_num++;
      }
    }
    
    dpatch_num=0;
    for (i=0;i<g_number_of_patches;i++) {
      if (g_test_pattern[i].patch_type == DENSITY_PATCH) {
	if (dpatch_num<7) {
	  g_test_pattern[i].patch_density=density[6-dpatch_num];
	  g_test_pattern[i].reflectance=reflect[6-dpatch_num];
	}
	else {
	  g_test_pattern[i].patch_density=density[20-dpatch_num];
	  g_test_pattern[i].reflectance=reflect[20-dpatch_num];
	}
	
	dpatch_num++;

	if (g_test_pattern[i].patch_density == white_density)
	  *white_patchp = i;
	else if (g_test_pattern[i].patch_density == black_density)
	  *black_patchp = i;
      }
    }
    
    free(density);
    free(reflect);

    density_order_correction = 1;
    
    }
  }
  
  g_white_patch = *white_patchp;

  define_corner_patches();

  if (fclose(fd) != 0)
    fprintf(stderr,"skew.c: Could not fclose fd: %s\n\n", strerror(errno));


  /* If modulation goes out of bounds, then stop with error */
  if (!bounded_modulation) { /* V6.3 */
    MTFPRINT("ERROR: target modulation in datafile is outside allowed range [0,1]\n")
    mtfexit(-1);
  }
  return (1);
}


/******************************************************************************/
/*									      */
/******************************************************************************/

/* Read in a scan line from a raw image file */

short raw_readscanline(unsigned char *buf, short i)
{
  int pos;

  if ((unsigned)i > g_scan_rows-1) return(-1);
  if (g_bytes_per_pixel != 1) return(-1);
  pos = g_scan_header + i*g_scan_cols;
  if (lseek(g_scan_image_file_id,pos,0) < 0) return(-1);
  if (read(g_scan_image_file_id,buf,(int)g_scan_cols) != (int)g_scan_cols)
    return(-1);
  else return(0);
}


/******************************************************************************/
/*									      */
/******************************************************************************/
/* Find which test_pattern indexes correspond to the four corner gray patches*/
/* This info allows more reliable angle and ppi computations.  New V5.0 mal  */
void define_corner_patches(void)
{
  int i;
  float ulex=0.0, urex=0.0, llex=0.0, lrex=0.0;
  float ul, ur, ll, lr;
  int uli=0, uri=0, lli=0, lri=0;
  int found = 0;
  int flag = 0;

  /* If Corners already marked, then restrict view to ordering just those 4 */
  if (g_corners_read == 4)  flag = 1;
  	  
  /* Find first density patch and set starting extrema. */
  for(i=0; !found && i<g_number_of_patches; i++) {
    if (g_test_pattern[i].patch_type == DENSITY_PATCH && 
    	    g_test_pattern[i].corner_flag == flag) {
      ulex = g_test_pattern[i].xul_in_mm + g_test_pattern[i].yul_in_mm;
      urex = g_test_pattern[i].xul_in_mm - g_test_pattern[i].yul_in_mm
 	   + g_test_pattern[i].width_in_mm;
      llex = - g_test_pattern[i].xul_in_mm + g_test_pattern[i].yul_in_mm
	                                   + g_test_pattern[i].height_in_mm;
      lrex = g_test_pattern[i].xul_in_mm + g_test_pattern[i].yul_in_mm
	 + g_test_pattern[i].width_in_mm + g_test_pattern[i].height_in_mm;
      uli = uri = lli = lri = i;
      found = 1;
    }
  }
  
  /* Scan rest of density patches and find ones that are furthest out. */
  for( ; i<g_number_of_patches; i++) {
    if (g_test_pattern[i].patch_type == DENSITY_PATCH && 
    	    g_test_pattern[i].corner_flag == flag) {
      ul = g_test_pattern[i].xul_in_mm + g_test_pattern[i].yul_in_mm;
      ur = g_test_pattern[i].xul_in_mm - g_test_pattern[i].yul_in_mm
	+ g_test_pattern[i].width_in_mm;
      ll = - g_test_pattern[i].xul_in_mm + g_test_pattern[i].yul_in_mm
	                                 + g_test_pattern[i].height_in_mm;
      lr = g_test_pattern[i].xul_in_mm + g_test_pattern[i].yul_in_mm
	+ g_test_pattern[i].width_in_mm + g_test_pattern[i].height_in_mm;

      if ( ul < ulex ) {
	ulex = ul;
	uli = i;
      }
      if ( ur > urex ) {
	urex = ur;
	uri = i;
      }
      if ( ll > llex ) {
	llex = ll;
	lli = i;
      }
      if ( lr > lrex ) {
	lrex = lr;
	lri = i;
      }
    }
    
  }
  g_test_pattern_corner[0] = uli;
  g_test_pattern_corner[1] = uri;
  g_test_pattern_corner[2] = lri;
  g_test_pattern_corner[3] = lli;

}

/******************************************************************************/
/*									      */
/******************************************************************************/
/* Locate the external pattern corners, based upon the 3 internal corners and
   the white patch location. V5.0 mal */
void internal_to_external_corners(double test_pattern_width, 
				  double test_pattern_height, short rotation) 
{
  double Delta_old[2],Delta_new[2];
  double ppi_horz, ppi_vert;
  float *cx[3], *cy[3];
  float *dx[3], *dy[3];
  int sx[3], sy[3];
  int i, cnr, corner_index;
  double width, height;
  double cos_ang, sin_ang;



  /* These are the two components of a vector denoting the change from the
     input points to the calculated corner points, if there were no skew.
     If there is skew, these components will be transformed to the proper
     values using a rotation matrix.  */

  cos_ang = cos(g_skew_avg);
  sin_ang = sin(g_skew_avg);

  /* Initialize pointers and signs to allow for loop processing of the 
     3 corners ll, ul, ur in order. */
  cx[0] = &g_test_pattern_xll;
  cy[0] = &g_test_pattern_yll;
  cx[1] = &g_test_pattern_xul;
  cy[1] = &g_test_pattern_yul;
  cx[2] = &g_test_pattern_xur;
  cy[2] = &g_test_pattern_yur;
 
  dx[0] = &g_input_test_pattern_xll;
  dy[0] = &g_input_test_pattern_yll;
  dx[1] = &g_input_test_pattern_xul;
  dy[1] = &g_input_test_pattern_yul;
  dx[2] = &g_input_test_pattern_xur;
  dy[2] = &g_input_test_pattern_yur;

  sx[0]=-1;
  sy[0]=1;
  sx[1]=-1;
  sy[1]=-1;
  sx[2]=1;
  sy[2]=-1;

 /* Check that corner patch dimensions are consistent with the 
     default patch. */

  /* Orient ppi values correctly for non-normalized image */
  if (g_scan_orient == LAND){
    ppi_horz = g_ppi_horz;
    ppi_vert = g_ppi_vert;
  }
  else{
    ppi_horz = g_ppi_vert;
    ppi_vert = g_ppi_horz;
  }

  /* Process each corner in order, 0=ll, 1=ul, 2=ur */
  cnr = rotation;  /* Start at the lower left corner */
  for(i=0; i<3; i++, cnr++) {
    cnr %= 4;
    corner_index = g_test_pattern_corner[cnr];
    width = g_test_pattern[corner_index].xul_in_mm;
    height = g_test_pattern[corner_index].yul_in_mm;
    if (cnr==0 || cnr==1) height += g_test_pattern[corner_index].height_in_mm;
    if (cnr==0 || cnr==3) width += g_test_pattern[corner_index].width_in_mm;
    if (cnr==1 || cnr==2) width = g_target_width_in_mm - width;
    if (cnr==3 || cnr==2) height = g_target_height_in_mm - height;

      width = width / MM_PER_INCH;
      height = height / MM_PER_INCH;
      /* fprintf(stderr, "Patch non-def size %f %f vs %f %f\n", 
         width, height, g_patch_width, g_patch_height); <--for debug only */
    
    /* Dimensions re-oriented to input image pose */
    if (g_scan_orient == LAND){
      Delta_old[0] = width;
      Delta_old[1] = height;
    }
    else{ /* these are reverse if the image is oriented in portrait position */
      Delta_old[1] = width;
      Delta_old[0] = height;
    }

    /* Set vector direction according to corner position */
    Delta_old[0] *= sx[i];
    Delta_old[1] *= sy[i];
    
    /* Rotate vector based on skew angle  and then scale */
    Delta_new[0] = (Delta_old[0]*cos_ang-Delta_old[1]*sin_ang)* ppi_horz;
    Delta_new[1] = (Delta_old[0]*sin_ang+Delta_old[1]*cos_ang)* ppi_vert;

    /* Move input corner position out to external corner */
    *(cx[i]) = *(dx[i]) + Delta_new[0];
    *(cy[i]) = *(dy[i]) + Delta_new[1];

    if (i==0) { /* ll corner just found, so compute lr as well */
      /* Compute offset from the ll input corner position */
      /* This computation assumes that the lr external corner is 
	 at same vertical position as the ll external corner. */

      if (g_scan_orient == LAND) {
	width = (g_target_width_in_mm / MM_PER_INCH) - width;
	Delta_old[0] = width;
	Delta_old[1] = height;
      }
      else {
	height = (g_target_height_in_mm / MM_PER_INCH) - height;
	Delta_old[1] = width;
	Delta_old[0] = height;
      }

      Delta_new[0] = (Delta_old[0]*cos_ang-Delta_old[1]*sin_ang)* ppi_horz;
      Delta_new[1] = (Delta_old[0]*sin_ang+Delta_old[1]*cos_ang)* ppi_vert;

      g_test_pattern_xlr = g_input_test_pattern_xll + Delta_new[0];
      g_test_pattern_ylr = g_input_test_pattern_yll + Delta_new[1];
    }

  }
}

/******************************************************************************/
/*									      */
/******************************************************************************/

/* Determine the orientation of the input image.  The orientation is determined
	by the location of the white patch */

short get_orientation(TIFF *tif, unsigned int white_patch)
{
  double sum[2];
  unsigned int count[2];
  int j;
  unsigned char *buf;
  double cossa,sinsa,col,row;
  short r[4],c[4];
  double width_in_mm,height_in_mm,box;
  double x_offset, y_offset;
  short rstart,rend,cstart,cend;
  int corner_patch;

  /* The "box" is the area within the patch from which we read pixel values tp
     determine which of two possible patches is the white patch. Make the
     "box" small enough so that it will be within the patch no matter
     how skewed the image (but big enough to give a good average reading) */

  width_in_mm = g_test_pattern[white_patch].width_in_mm;
  height_in_mm = g_test_pattern[white_patch].height_in_mm;
  box = (width_in_mm < height_in_mm) ? width_in_mm:height_in_mm;
  box  = box*4.0/9.0;

  /* These will be overwritten later in the routine "compute_patch_boundaries";
     for now we store the dimensions of the "box" here. This is an imaginary
     box centered within the white patch on the actual target. This box will
     be rotated and translated to fall within the patch on the actual image. */

  g_test_pattern[white_patch].x[0] =
    g_test_pattern[white_patch].x[2] =
	(width_in_mm-box)/2.0;
  g_test_pattern[white_patch].x[1] =
    g_test_pattern[white_patch].x[3] =
	(width_in_mm+box)/2.0;

  g_test_pattern[white_patch].y[0] =
    g_test_pattern[white_patch].y[1] =
	(height_in_mm-box)/2.0;
  g_test_pattern[white_patch].y[2] =
    g_test_pattern[white_patch].y[3] =
	(height_in_mm+box)/2.0;

  count[0]=count[1]=0;
  sum[0]=sum[1]=0.0;
  
  buf = NULL;
  if (g_scan_image_file_id)
    buf = (unsigned char*) calloc((unsigned int)g_scan_cols,sizeof(unsigned int));
#ifdef USE_TIFF
  else 
    buf = (unsigned char*) calloc((unsigned int)TIFFScanlineSize(tif),sizeof(unsigned int));
#endif
  if ( buf == NULL ) {
    fprintf(stderr,
	    "\nCan't allocate memory for scanline buffer!\n\n");
    mtfexit(-1);
  }

  /* Use corner patch[1] (= UR target patch) as the reference position.
     This patch is typically closest to the white patch, errors due to
     skew or ppi approximations are minimized. */
  corner_patch = g_test_pattern_corner[1];
    
  x_offset = g_test_pattern[white_patch].xul_in_mm -
    g_test_pattern[corner_patch].xul_in_mm;

  y_offset = g_test_pattern[white_patch].yul_in_mm -
    (g_test_pattern[corner_patch].yul_in_mm + 
     g_test_pattern[corner_patch].height_in_mm);
    
  /* For an image in portrait orientation, the white patch can be on the upper
     left or the lower right.  We look at both cases and call the largest
     average pixel value the "white patch". */

  if (g_scan_orient == PORT) {
    /* upper left case  - rotate by PI/2 - g_skew_avg and translate
       to lower left corner */

    cossa = cos(-PI_2 + g_skew_horz);
    sinsa = sin(-PI_2 + g_skew_horz);

    for(j=0;j<4;j++) {
      row=(g_test_pattern[white_patch].y[j]+y_offset)/MM_PER_INCH;
      col=(g_test_pattern[white_patch].x[j]+x_offset)/MM_PER_INCH;
      
      /* Rotate and then scale and then translate */
      c[j] = NINT(g_ppi_vert*(col*cossa-row*sinsa)+g_input_test_pattern_xul);
      r[j] = NINT(g_ppi_horz*(col*sinsa+row*cossa)+g_input_test_pattern_yul);
    }

    /* The box is small enough so that we don't have to worry about skew - 
       i.e. if r3 is left or right of r1, we will still be well within the 
       overall patch. We just take a sample from within the patch and find 
       an average pixel value. */
    if (g_skew_horz >= 0) {
      rstart = r[1]; rend = r[2]; cstart = c[0]; cend = c[3];
    }
    else {			/* negative skew angle */
      rstart = r[3]; rend = r[0]; cstart = c[1]; cend = c[2];
    }
    do_white_patch(buf,tif,sum,count,rstart,rend,cstart,cend);

    /* Since lr corner is not computed shift basis of computation to
       the ll corner (corner 2 in cyclic testpattern corner order) */
    corner_patch = g_test_pattern_corner[2];
    
    x_offset = g_test_pattern[white_patch].xul_in_mm -
      g_test_pattern[corner_patch].xul_in_mm;

    y_offset = g_test_pattern[white_patch].yul_in_mm -
      g_test_pattern[corner_patch].yul_in_mm;
  
    cossa = cos(PI_2 + g_skew_avg);
    sinsa = sin(PI_2 + g_skew_avg);

    for(j=0;j<4;j++) {
      row=(g_test_pattern[white_patch].y[j]+y_offset)/MM_PER_INCH;
      col=(g_test_pattern[white_patch].x[j]+x_offset)/MM_PER_INCH;

      /* Normally would be offset from xlr and ylr, but offset basis
	 was changed to ll since those two values are not yet computed */
      c[j] = NINT(g_ppi_vert*(col*cossa-row*sinsa)+g_input_test_pattern_xll);
      r[j] = NINT(g_ppi_horz*(col*sinsa+row*cossa)+g_input_test_pattern_yll);
    }

    if (g_skew_avg >= 0) {
      rstart = r[2]; rend = r[1]; cstart = c[3]; cend = c[0];
    }
    else {
      rstart = r[0]; rend = r[3]; cstart = c[2]; cend = c[1];
    }
    do_white_patch(buf,tif,&sum[1],&count[1],rstart,rend,cstart,cend);
  }
  else {			/* LAND */
    cossa = cos(g_skew_vert);
    sinsa = sin(g_skew_vert);

    for(j=0;j<4;j++) {
      row=(g_test_pattern[white_patch].y[j]+y_offset)/MM_PER_INCH;
      col=(g_test_pattern[white_patch].x[j]+x_offset)/MM_PER_INCH;

      c[j] = NINT(g_ppi_horz*(col*cossa-row*sinsa)+g_input_test_pattern_xur);
      r[j] = NINT(g_ppi_vert*(col*sinsa+row*cossa)+g_input_test_pattern_yur);
    }

    if (g_skew_vert >= 0) {
      rstart = r[0]; rend = r[3]; cstart = c[2]; cend = c[1];
    }
    else {
      rstart = r[1]; rend = r[2]; cstart = c[0]; cend = c[3];
    }
    do_white_patch(buf,tif,sum,count, rstart,rend,cstart,cend);

    cossa = cos(PI - g_skew_avg);
    sinsa = sin(PI + g_skew_avg);

    for(j=0;j<4;j++) {
      row=(g_test_pattern[white_patch].y[j]+y_offset)/MM_PER_INCH;
      col=(g_test_pattern[white_patch].x[j]+x_offset)/MM_PER_INCH;

      c[j] = NINT(g_ppi_horz*(col*cossa-row*sinsa)+g_input_test_pattern_xll);
      r[j] = NINT(g_ppi_vert*(col*sinsa+row*cossa)+g_input_test_pattern_yll);
    }

    if (g_skew_avg >= 0) {
      rstart = r[3]; rend = r[0]; cstart = c[1]; cend = c[2];
    }
    else {
      rstart = r[2]; rend = r[1]; cstart = c[3]; cend = c[0];
    }
    do_white_patch(buf,tif,&sum[1],&count[1],rstart,rend,cstart,cend);
  }

  /* Free the buffer.  This was after the switch statement below, which was
     meaningless.  DJB 6/7/95 */
  free(buf);

  /* find the largest average value - this should be the white patch */
  sum[0] /= (double)count[0];
  sum[1] /= (double)count[1];

  switch (g_scan_orient) {
  case PORT: 
    if (sum[0] > sum[1]) return(TOP);
    return(BOTTOM);
    /*NOTREACHED*/
    break;
  case LAND:
    if (sum[0] > sum[1]) return(RIGHT);
    return(LEFT);
    /*NOTREACHED*/
    break;
  default:
    fprintf(stderr, "Impossible situation in get_orientation!\n");
    fprintf(stderr, "Orientation neither portrait nor landscape!\n");
    mtfexit(-1);
  }
  /*NOTREACHED*/
  return(TOP);
}

/******************************************************************************/
/*									      */
/******************************************************************************/
/* Computes the angle from target vector to the image vector  */
double angle_between_vectors(double tarvec_x, double tarvec_y,
			     double imgvec_x, double imgvec_y) 
{
  double dot_product;
  double mod_img, mod_tar;
  int sgn;

  dot_product = imgvec_x*tarvec_x + imgvec_y*tarvec_y;
  mod_img = sqrt(imgvec_x*imgvec_x + imgvec_y*imgvec_y);
  mod_tar = sqrt(tarvec_x*tarvec_x + tarvec_y*tarvec_y);

  if ( (imgvec_x * tarvec_y) > (imgvec_y * tarvec_x) ) sgn = -1;
  else sgn = 1;
  
  return ( sgn * acos(dot_product/(mod_img*mod_tar)) );
  
}

/******************************************************************************/
/*									      */
/******************************************************************************/
/*  Use actual corner patch dimensions to compute the skew.  This is important
    if one of the corner patches has very different dimensions. 
    Uses a least-squares fit to compute the ppi values for the new skew angle
*/
void compute_correct_skews(short rotation) 
{
  double cnrx[4], cnry[4];
  int cnr, i;
  int corner_index;
  double *ppi_horz, *ppi_vert;
  double u[3],v[3],w[3],z[3], s,c;
  double sumu, sumv, ssu, ssv, suw, svz, sumw, sumz;
  double woffset, zoffset;

  /* Find target reference positions of each of the 3 capture corners.
     In order: ll corner (0), ul corner (1), ur corner (2). */
  cnr=(int)rotation;
  for(i=0; i<3; i++, cnr++) {
    cnr %= 4;
    corner_index = g_test_pattern_corner[cnr];
    cnrx[i] = g_test_pattern[corner_index].xul_in_mm;
    cnry[i] = g_test_pattern[corner_index].yul_in_mm;
    if (cnr==0 || cnr==1) cnry[i] += g_test_pattern[corner_index].height_in_mm;
    if (cnr==0 || cnr==3) cnrx[i] += g_test_pattern[corner_index].width_in_mm;
  }
  
  /* Selected Image pixels, offset so UL corner is 0,0
      W,Z = Input image space */
  w[0] = (g_input_test_pattern_xll - g_input_test_pattern_xul);
  z[0] = (g_input_test_pattern_yll - g_input_test_pattern_yul);
  w[1] = 0.0;
  z[1] = 0.0;
  w[2] = (g_input_test_pattern_xur - g_input_test_pattern_xul);
  z[2] = (g_input_test_pattern_yur - g_input_test_pattern_yul);
	 
  if( rotation == LEFT ) { 
    /* Flip vectors from input image data by 180 degrees  to avoid 
       +PI vs -PI angle calculations. So the LEFT computed angles are     
       computed relative to the normalized image orientation. */
 	for (i=0; i<3; i++) {
   	 	w[i] *= -1;
    		z[i] *= -1;
     }
  }

  g_skew_horz = angle_between_vectors(cnrx[2]-cnrx[1], cnry[2]-cnry[1], w[2], z[2]);
  g_skew_vert = angle_between_vectors(cnrx[0]-cnrx[1], cnry[0]-cnry[1], w[0], z[0]);

  g_skew_avg = (g_skew_horz+g_skew_vert)/2.0;

  /* Set variables that are named according to input image space */
  if (g_scan_orient == LAND) {
    ppi_horz = &g_ppi_horz;
    ppi_vert = &g_ppi_vert;
  }
  else {
    ppi_vert = &g_ppi_horz;
    ppi_horz = &g_ppi_vert;
  }

  /* Compute the best fit ppi's based upon the average skew */
  /* Formulas come from assuming some error in input corner
     points, and attempting to minimize it in a scale/rotation 
     only system. Ideally one would minimize over the skew
     angle as well, but the formula is nasty. */
    
  c = cos(g_skew_avg);
  s = sin(g_skew_avg);
  
  /* Apply the rotation to the 3 corner positions in target space. */
  for(i=0; i<3; i++) {
	u[i] = (c*cnrx[i]-s*cnry[i])/MM_PER_INCH;
	v[i] = (c*cnry[i]+s*cnrx[i])/MM_PER_INCH;
  }
  /* Linear regression for best fit, [w z] = ppi * [u v] + offset */
  sumu = sumv = sumw = sumz = 0.0;
  ssu = ssv = 0.0;
  suw = svz = 0.0;
  for(i=0; i<3; i++) {
	sumu += u[i];
	sumv += v[i];
	sumw += w[i];
	sumz += z[i];
	ssu += u[i]*u[i];
	ssv += v[i]*v[i];
	suw += u[i]*w[i];
	svz += v[i]*z[i];
  }
  *ppi_horz = (3.0*suw - sumu*sumw)/(3.0*ssu-sumu*sumu);
  *ppi_vert = (3.0*svz - sumv*sumz)/(3.0*ssv-sumv*sumv);
 
   woffset = (sumw*ssu - sumu*suw)/(3.0*ssu-sumu*sumu);
   zoffset = (sumz*ssv - sumv*svz)/(3.0*ssv-sumv*sumv);

 /* printf("Tx_offset = %.1f pixels   Ty_offset = %.1f pixels\n",
	woffset, zoffset); */

 /* Rotate angles to eliminate basic 90 degree component */
  switch(rotation){
  case TOP:
    g_skew_horz += PI_2;
    g_skew_vert += PI_2;
    g_skew_avg  += PI_2;
    break;
  case BOTTOM:
    g_skew_horz -= PI_2;
    g_skew_vert -= PI_2;
    g_skew_avg  -= PI_2;
    break;
  case LEFT:  /* computed with normalized info */
  case RIGHT: /* orientation is already normal */
    break;
  }
  
  /* avoid divide by zero */
  if (g_skew_avg == (double)0.0) g_skew_avg = (double)0.000001; 

  if (((int)g_ppi_horz > (int)PPI_THRESHOLD) &&
      ((int)g_ppi_vert > (int)PPI_THRESHOLD))
      g_target_res = PPI_HR_BASE;
  else
      g_target_res = PPI_BASE;
}



/******************************************************************************/
/*									      */
/******************************************************************************/

/* Print out double comment info in the patch data file. Describing upright   */
/* position. */
/* Return flag indicating whether any info found. */

void print_upright_info(char *data_filename)
{
  FILE *fd;
  char line[82];
  char *err;

  fd = fopen(data_filename,"r");
  if (fd == NULL) {
    fprintf(stderr,"\nCould not open %s\n",data_filename);
    fprintf(stderr,"%s\n\n", strerror(errno));
    mtfexit(-1);
  }

  do {
       err = fgets(line,81,fd);
  }  while (err != NULL && (line[0] != '#' || line[1] != '#')); 

  if (err == NULL) {
   /* No double comment lines found. */
   fclose(fd); 
   return;  /* v6.3 added to avoid printing any lines */
  }
  
  do {
    fprintf(stderr, "%s", line);    
    err = fgets(line,81,fd);
  }  while (err != NULL && line[0] == '#' && line[1] == '#'); 

  fclose (fd);
}

