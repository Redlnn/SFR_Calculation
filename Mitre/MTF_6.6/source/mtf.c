/**********************************************************************
                                Notice
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
/* Revisions to mtf.c                                                        */
/*                                                                           */
/* V1.1 General Clean-up  brp 11-MAR-94                                      */
/*                                                                           */
/* V1.2 Corrected args to "open" by adding a 2nd, 0 arg  brp 15-MAR-94       */
/*                                                                           */
/* V1.3 NINT routine added for portability  gtk 7-JUN-94                     */
/*      changed usage message.  brp 9-JUN-94                                 */
/*                                                                           */
/* V2.0 Added code for tiff polarity = g_photometric (black=255 or black=0)  */
/*                                                              gtk 22-JUL-94*/
/*                                                                           */
/* V2.1 Added subroutines round(), magnitude(), dabs(), sign().              */
/*                                       gtk 23-AUG-94                       */
/*                                                                           */
/* V2.3 Close image files after they are read in so we don't get an error    */
/*       too many open files if processing multiple images.                  */
/*      Target data filename is now required for each image (no default or   */
/*       global setting).                                                    */
/*      Revised printout that prompts user for information.                  */
/*      Revised printout for -h switch.                                      */
/*      Revised output printout.                                             */
/*      Moved switch parsing from get_args() to get_switches(), which is new.*/
/*      Revised problem report printout.                                     */
/*      Added a check for the existence of the data file.                    */
/*                                                      DJB 6/12/95          */
/*                                                                           */
/* V3.0 Added -p (printer MTF) switch and code.                              */
/*      Indentation changes.                                                 */
/*                                                      DJB 8/16/95          */
/*      If -p, don't compare to 500 ppi scanner spec.                        */
/*      If -p, adjust computed ppi based on user-input printer ppi.          */
/*      Initialize g_printer_ppi and g_scan_image_file_id.                   */
/*      Conditionalize on USE_TIFF to not depend on tiff code.               */
/*      Require user input header size on non-tiff images.                   */
/*                                                      DJB 12/28/95         */
/*      For -d switch, rotate output image back to orientation of original   */
/*       input image.                                                        */
/*      No longer print out calculated coordinates of outer corners of target*/
/*       image.                                                              */
/*      If number of problems to report exceeds maximum number of problems   */
/*       program can record, say so at bottom of problem report.             */
/*      Only call compute_skew_angles once per pass.                         */
/*      User-input inner corner coordinates are now float.                   */
/*      For DOS, open image file in BINARY mode.                             */
/*                                                      DJB 1/22/96          */
/*      Give user another chance if an invalid image or data file name was   */
/*      input.                                                               */
/*                                                      DJB 2/16/96          */
/*      Fixed bug--resolution was only being checked in horizontal direction.*/
/*                                                      DJB 2/19/96          */
/*      Removed some blank printout lines and "Problem report follows." line.*/
/*                                                      DJB 2/20/96          */
/*      Delete g_orient_debug code.                                          */
/*                                                      DJB 2/22/96          */
/*                                                                           */
/* V3.1 Removed copyright and added data rights notice.                      */
/*      Changed print_copyright to print_development_notice.                 */
/*      Added -r switch for data rights notice.                              */
/*                                                      DJB 6/25/96          */
/*      Added messages for bad file closes.             DJB 8/14/96          */
/*                                                                           */
/*  V3.2 Power Macintosh version created.               SB & NBN 2/98        */
/*                                                                           */
/*  V3.2 updated Win95/NT & Unix versions to v3.2                            */
/*       & updated Mac memory allocation method         NBN 10/98            */
/*       (chgs in mtf.c, density.c, sine.c, skew.c, mtf.h)                   */
/*      Other chgs in v3.2: 
                menu display replaces command line entry.
                new -c option allows point-to-point st.line conversion from 
                image space gray to target space reflectance.
                cleaned-up some unneeded, confusing I/O statements             
         Mac-specific:
                Mac v3.2 compiled without TIFF format reader modules         */
/*                                                                           */
/*    V4.0  Added new contract numbers.                                      */
/*          Added support for 1000 ppi.                                      */
/*          Added "n" switch: do not compare output to spec.                 */
/*          Added (hidden) "streamline" mode for debuging                    */
/*              assumes all arguments are given in command line.             */
/*          Added test for image polarity                                    */
/*          get_target_data follows read-in density values and sets which    */
/*              patch is white and which patch is black                      */
/*          General cleanup to adhere to gcc's -ansi, -pedantic, and -Wall   */
/*               switches                                                    */
/*          Change orientation error to V/H scale warning, if V/H differ by  */
/*               more than 2% then warn user.                                */
/*                                                      WE 9/99              */
/*                                                                           */
/*                                                                           */
/*  V4.1  For printer option P ("g_printer" is active), utilize a new        */
/*        parameter: PRINTER_TARGET_BASE (defined in mtf.h) which replaces   */
/*        all occurrences of PPI_BASE & PPI_HR_BASE that were associated with*/
/*        calculations involving g_printer_ppi.  PPI_BASE, PPI_HR_BASE are   */
/*        validly used elsewhere in program, but are not relevant to         */
/*        calculations involving g_printer_ppi.                              */
/*                                                     nbn 9/00              */
/*                                                                           */
/* V5.0 Several key changes/additions:                                       */
/*      1. Check for a common labelling error within the datafile.           */
/*         If present, info is automatically corrected at runtime            */
/*      2. When 'b' option is selected, the input corner points are used     */
/*         as starting points for an auto corner selection process.          */
/*      3. Orientation is discovered during initialization                   */
/*      4. More extensive help info added.                                   */
/*      5. S-Curve fitting added.                                            */
/*      6. Piecewise fit added as a user choice, and spline fit now          */
/*	   reverts to s-curve fit when there are errors.                     */
/*                                                      mal, 11/03           */
/*                                                                           */
/* V5.5     Fixed four bugs and added two features.                          */
/*          1. Fixed rare sine calculation bug (peak mtf=0) from v4.52       */
/*          2. Fixed bug introduced in v5.0 preventing linear regression     */
/*             usage when there are fewer than 5 density patches             */
/*          3. Made alias problems report even when option 'n' selected      */
/*          4. Allow single 4-point density average to be retained even if   */
/*             another density average pair is found and removed.            */
/*          Added the ability to process bar target images that do not       */
/*	     contain density patches.  An extra user response                */
/*           is required for each bar target image, giving target orientation*/
/*           (bar target indicated by 'b' patch type in target data file)    */
/*          Target data file contents checked for consistency. Patch         */
/*          measurements are compared to stated target width/height. A       */
/*          problem is reported if patch measurements exceed target dims.    */
/*	    Target data file reading moved up in logic chain. It now occurs  */
/*	     as soon as the filename is entered.			     */
/*          Problem report initialization moved earlier in logic.            */
/*          Also cleaned up uninitialized/unused warnings from several files.*/
/*                                                      mal, 02/04           */
/*                                                                           */
/* V5.6     Corrected interactive target orientation query text.             */
/*          All extended output written to MTFOUT only, not stdout.          */
/*                                                      mal, 06/04           */
/*                                                                           */
/* V5.8     Cleaned up code.  Sections no longer used were removed.          */
/*                                                      mal, 09/04           */
/*                                                                           */
/* V5.8.1   Fixed Scurve MTF printout (1.0's that should show 0.0).          */
/*          Changed include info to allow CodeWarrior v9 compilation         */
/*          Removed unused math functions, and moved used ones to            */
/*          different locale, to accommodate compilation.                    */
/*                                                      mal, 03/05           */
/* V5.8.2   Debug image written as a TIFF (or PGM if TIFF is not available)  */
/*                                                      mal, 05/05           */
/*                                                                           */
/* V6.0     Option 'g' added for digital camera assessment                   */
/*                                                      mal, 06/05           */
/* V6.2     Added non-linear reflectance data for bartarget mode   	     */
/*          Added PGM read capability                                        */
/*                                                      MAL, 08/06           */
/* V6.3     New options: 'j' compare output to PIV spec         	     */
/*                       'k' compare output to AppF spec                     */
/*          Removed option 'n' (don't compare to AppF spec) since that is    */
/*              now the default behavior when neither 'j' nor 'k' is chosen. */
/*          Added PIV spec checking to code.                                 */
/*                                                      MAL, 10/06           */
/* V6.4.01  Compilation without CodeWarrior, and using external TIFF library */
/*                                                      MAL, 09/12           */
/* V6.4.02  Enabled reading PGM with H/W on separate lines.		     */
/*                                                      MAL, 11/13           */
/* V6.5    Added colors and additional outlines to debug images.	     */
/*         Debug image outputs in BMP format if TIFF is not available.       */
/*                                                      AWU, 07/14           */
/*****************************************************************************/


/******************************************************************************/
/* This is the "MAIN" routine of the MTF program                              */
/******************************************************************************/
#include <string.h>
#include <math.h>


#ifdef USE_TIFF
#include <tiffio.h>
TIFF *dtif;  /* Debug TIFF file pointer */
#else
typedef char TIFF;
#endif

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <ctype.h>
#include <fcntl.h>
#include <errno.h>

#if !defined(MSDOS)
#include <unistd.h>    /* chdir() */
#include <libgen.h>	  /* dirname()*/
#endif

#if !defined(O_BINARY)
#define O_BINARY 0
#endif

#include "mtf.h" /* a lot of macros and defines */
#include "xterndef.h" /* a lot of global variables */

/* set this to 1 to produce "orient.raw", a debug image showing how the
        patches were read to determine the input image orientation. */
#define DO_ORIENT_DEBUG 0


static int too_many_problems;
static double g_FL, g_target_distance;

static unsigned char g_pgm;
static int read_pgm_header(char *);
static void write_debug_image(short rotation);
static void write_debug_bounds(short rotation);

FILE *g_ofp; /* Debug non-tiff image file pointer */
int main(int argc, char **argv)
{
  unsigned char original_scan_orient;
  char ans[10];
  char problem_string[82];
  unsigned int white_patch, black_patch, lowf_patch;

  short rotation;
  short i;

  TIFF* tif = NULL;
  double vh_scale_diff;

  int insert_page_break;

  /* Once-only initializations */
  g_scan_image_file_id = 0;

  /* Development notice to the terminal for a few seconds */ 
  print_development_notice();

  get_switches(argc, argv, image_filename, data_filename);              

  while(1 /*CONSTANTCONDITION*/) {      /* loop for 1 or more files to process */

    /* Once-per-pass initializations */
    need_spec_footnote = 0;
    g_printer_ppi = PRINTER_TARGET_BASE;
    too_many_problems = 0;
    g_LL_auto_crnr_err=g_UR_auto_crnr_err=g_UL_auto_crnr_err=0;
    g_compare = g_compare_userset;  

    /* Set curvefitting for this iteration to the global option request */
    if(g_user_requestcurve == 0 && g_printer ) g_followcurve = 9;
    else g_followcurve = g_user_requestcurve;

    /* If MTFOUT already exists, then we are processing another image,
       so we want a page break in MTFOUT so the report for this image
       starts on a new page in the printout. */
    if ((g_mtfout = fopen(MTFOUTNAME, "r")) == NULL)
      insert_page_break = 0;
    else {
      insert_page_break = 1;
      if (fclose(g_mtfout) != 0)
        fprintf(stderr, "mtf.c: Could not fclose g_mtfout: %s\n\n", strerror(errno));
    }

    /* always append output to this file */
    if ((g_mtfout = fopen(MTFOUTNAME,"a")) == NULL) {
      fprintf(stderr,"Can't open %s\n",MTFOUTNAME);
      mtfexit(1);
    }

    if (insert_page_break)
      /* Write a page break to MTFOUT */
      fprintf(g_mtfout, "\n");

    /* command line args */
#ifdef USE_TIFF
    if (tif != NULL)
      {
        TIFFClose(tif);
        tif = NULL;             /* Next image might not be a TIFF */
      }
    else
#endif
      if (g_scan_image_file_id != 0)
        {
          if (close(g_scan_image_file_id) != 0)
            fprintf(stderr, "mtf.c: Could not close g_scan_image_file_id: %s\n\n", strerror(errno));
          g_scan_image_file_id = 0; /* Next image might not be a raw */
        }

    /* Initialize early, so problems can be logged v5.5 */
    g_problem_count = 0;		/* # of problems in our problem report */
    g_IQS_problem_count = 0;	


    /* set which patch is the black patch; set which patch is the white patch
       white patch: used in image polarity test once we determine the mean grey
       levels of each patch  */
    /* Set which sine patch has lowest frequency. Used to check image 
       polarization for bar targets  V5.5 */
    get_args(image_filename,data_filename,&tif,&white_patch,&black_patch,&lowf_patch);

    /* see 'skew.c' module for this routine */
    rotation = initialize(data_filename,white_patch,tif);

    /* Check to make sure that PIV check isn't used on 1000ppi data */
    if (g_target_res == PPI_HR_BASE && g_compare == PIV) {
      sprintf(problem_string,
	      "PIV spec does not apply to 1000ppi data.\nSpec checking not applied to this image\n\n");
      put_problem(problem_string, 0);
      g_compare = 0;
    }

    /* read the image in - we store it internally in a standard way */
    read_in_image(tif,rotation,image_filename);

    original_scan_orient = g_scan_orient;

    g_scan_orient = LAND;       /* normalized image is landscape */

    /* If the resolution is out of range say so in the problem report.
       If we're way off base the user probably indicated a vertical (portrait)
       orientation instead of a horizontal (landscape) orientation (or conversely).
       Problem reports are stored in a string array and output later. */

    if (!g_printer) {
      int ppi_error = PPI_ERROR;
      if (    ((g_compare) && (g_target_res == PPI_HR_BASE))
	   || (g_compare == PIV)) {
          ppi_error = PPI_HR_ERROR;
      }
      if (   (abs((int)g_ppi_horz - g_target_res) > (int)ppi_error)
          || (abs((int)g_ppi_vert - g_target_res) > (int)ppi_error)) {
        sprintf(problem_string,
                "Resolution outside the %d-%d ppi spec** range\n",
                g_target_res-ppi_error,g_target_res+ppi_error);
        put_problem(problem_string, IQS);
        put_problem("\n", IQS);
        need_spec_footnote = g_target_res;
      }
    }
    vh_scale_diff = g_ppi_vert/g_ppi_horz;
    if (vh_scale_diff > 1.02 || vh_scale_diff < 0.98) {
      sprintf(problem_string,
              "V & H Scales Differ by > 2 percent \n");
      put_problem(problem_string, IQS);
      put_problem("\n", IQS);
    }


    /* print the user entered data (corner coords, etc.) in the output */
    print_header(image_filename,data_filename,original_scan_orient);

    /* the x,y coords of the corners of each patch */
    compute_patch_boundaries();

    if (g_bartarget) {
      g_high_throw_out = g_max_pixel_value;
      g_low_throw_out = 0;
      if( g_nonlinear == 0) {
       /* Set up defaults for bar target processing. 
	  1. Flag for straight line fit.
	  2. Linear fit parameters.
       */
       g_followcurve = 6;
       g_slope = (double)g_max_pixel_value;
       g_intercept = 0.0;
      }  else {
       /* Set up nonlinear bar target processing. */
       g_followcurve = 5;
      }
        /* test for image polarity on a bar target V5.5 */
       test_for_bar_target_polarity(lowf_patch);
    }
    else {
        /* linear regression on the densities from the density patches; 
	   also prints out data into output */
        compute_densities(white_patch, black_patch);

        /* test for image polarity --WE 9/99 V3.3 */
        test_for_image_polarity(white_patch);
    }


    compute_sines();            /* also prints data into output */

    /* log any problems we encountered into output */
    if (g_reversepolarity)
      MTFPRINT("\nNOTE: Original image polarity was reversed by user request.\n")

    print_problems();

    /* if the user wanted histograms too */
    if (g_extended) {
      fprintf(stdout, "\n\nEXTENDED output printed to MTFOUT file includes:\n");
      fprintf(stdout, "\tSine Patch data\n");
      fprintf(stdout, "\tAliasing data\n");
	  fprintf(stdout, "\tAliasing summary\n");
      print_sine_stats();
      fprintf(stdout, "\tDensity Patch data\n");
      print_density_stats();
      fprintf(stdout, "\tGray-to-Reflectance curve\n\n");
      if(g_bytes_per_pixel==1) {
        if (g_followcurve==6)
          fprintf(g_mtfout,"GRAY-TO-REFLECTANCE VIA LINEAR,LEAST SQUARES REGRESSION LINE\n");
        else if (g_followcurve==7)
	  fprintf(g_mtfout, "GRAY-TO-REFLECTANCE VIA UNSMOOTHED SPLINE\n");
        else if (g_followcurve==8)
          fprintf(g_mtfout, "GRAY-TO-REFLECTANCE VIA SMOOTHED SPLINE\n");
        else if (g_followcurve==9)
	  fprintf(g_mtfout, "GRAY-TO-REFLECTANCE VIA S-CURVE FIT\n");
        else if (g_followcurve==4 || g_followcurve==5)
          fprintf(g_mtfout, "GRAY-TO-REFLECTANCE VIA PIECEWISE LINEAR INTERPOLATION\n");

	fprintf(g_mtfout, "\nGray     \tReflectance\n\n");
	
        for (i = 0;  i <= (int)g_max_pixel_value;  i++) {
	  if( g_bartarget && g_nonlinear && bad_grey_check((double)i, 0) ) 
	    fprintf(g_mtfout, "%d\t\tUNDEFINED\n",i);
	  else
	    fprintf(g_mtfout, "%d\t\t%f\n",i, GreyToReflectance((double)i));
	}
      }
    }

    /* if the user wants a couple of images showing what we looked at */
    if (g_debug)
      write_debug_image(rotation);




    if (g_test_pattern != NULL) free(g_test_pattern);
    if (g_image_array != NULL) free(g_image_array);
    if (g_image_array_modified != NULL)free(g_image_array_modified);




    MTFPRINT("\n")
    MTFPRINT2("         End of Report for %s\n",image_filename)
    MTFPRINT3("%s%s\n","****************************************",
              "****************************************")
    MTFPRINT("\n")
    MTFPRINT("\n")

    fflush(g_mtfout);
    if (fclose(g_mtfout) != 0)
      fprintf(stderr, "mtf.c: Could not fclose g_mtfout: %s\n\n", strerror(errno));
    g_streamline_args=0;
    fprintf(stderr,"\nProcess another image [y/n]?\t");
    scanf("%s",ans);
    fprintf(stderr,"\n");
        
    if(strcmp(ans,"y") && strcmp(ans,"Y")) break;
  }
  fprintf(stderr,"\n");

    if (argc == 1)
        fprintf(stderr,"Program ended.");

  mtfexit(1);
  return(1);
}                               /* end of main */

/******************************************************************************/
/*                                                                            */
/******************************************************************************/



void get_switches(int argc, char **argv, char *image_filename, char *data_filename)
{
  char switch_entry = ' ';
  int count = 0;
  int read_value = 0;

  g_streamline_args=g_print_average_mtf=g_debug=g_extended=0;
  g_printer=g_user_requestcurve=g_reversepolarity=g_bartarget=0;
  g_refine_corners=g_digital_camera=0;
  g_compare_userset=0;
 
    if (argc == 1) {
#if !defined(MSDOS)
char *fnamep;
fnamep = dirname(argv[0]);
        if (argc == 1 && fnamep != NULL && fnamep[0] == '/') {
     /*  printf ("Changed to directory of executable:  %s\n\n", fnamep); */
            chdir((const char *)fnamep);
        }
#endif
    }
       /* This is a COMMAND LINE mode for tif images, i.e.,
     read in all args from command line rather than waiting for questions;
     usage: mtf image.tif image.dat upper_left_col upper_left_row
            upper_right_col upper_right_row lower_left_col lower_left_row
            scan_orientation
     But jump around this section if mac or windows compile is via Metrowerks CodeWarrior,
     because it is assumed the CodeWarrior is using the sioux gui, which is not command line (incompatibility) */
  if (argc > 9) {
      argv++;
      strcpy(image_filename,argv[0]); argv++;
      strcpy(data_filename,argv[0]);  argv++;
      g_input_test_pattern_xul = (float)atof(argv[0]); argv++;
      g_input_test_pattern_yul = (float)atof(argv[0]); argv++;
      g_input_test_pattern_xur = (float)atof(argv[0]); argv++;
      g_input_test_pattern_yur = (float)atof(argv[0]); argv++;
      g_input_test_pattern_xll = (float)atof(argv[0]); argv++;
      g_input_test_pattern_yll = (float)atof(argv[0]); argv++;
      if ((*argv[0] == 'h') || (*argv[0] == 'H'))
         g_scan_orient = LAND;
      else
         g_scan_orient = PORT;
      argv++;
 
      if (argc > 10){
        /* Determine target rotation info in degrees */
        read_value = atoi(argv[0]);
 
        read_value *= 10;
        if (read_value == 10) read_value = 180;
        g_upright = read_value;
 
       read_value = 0;
      }
      g_streamline_args = 1;
  }
 
 
  fprintf(stderr,"  a      Output average MTF (in addition to peak MTF) \n");
  fprintf(stderr,"  b      Auto refine input corner coordinates \n");
  fprintf(stderr,"  c      Select I/O curve fit type for sine tgt (instead of autofit) \n");
#ifdef USE_TIFF
  fprintf(stderr,"  d      Create diagnostic image (_box.tif) \n");
#else
  fprintf(stderr,"  d      Create diagnostic image (_box.pgm) \n");
#endif
  fprintf(stderr,"  e      Verbose output \n");
  fprintf(stderr,"  f      Reverse image polarity \n");
  fprintf(stderr,"  g      Digital camera assessment \n");
  fprintf(stderr,"  j      Compare output to PIV spec \n");
  fprintf(stderr,"  k      Compare output to App. F spec \n");
  /*  fprintf(stderr,"  n      Don't compare output to App. F spec\n"); */
  fprintf(stderr,"  p      Printer MTF assessment \n");
  fprintf(stderr,"  r      Help & Program Notice \n\n");
  fprintf(stderr,"Enter run options...or...Press RETURN if none:  ");
 
  do
    {
      count++;
      read_value = scanf("%c", &switch_entry);
      switch (switch_entry)
        {
          case 'a':
          case 'A':
            g_print_average_mtf=1;
            break;
          case 'b':
          case 'B':
            g_refine_corners = 1;
            break;
          case 'c':
          case 'C':
            g_user_requestcurve=1;
          break;
          case 'd':
          case 'D':
            g_debug = 1;
            break;
          case 'e':
          case 'E':
            g_extended = 1;
            break;
          case 'g':
          case 'G':
           g_digital_camera = 1;
           g_compare_userset = 0;
           break;
          case 'j':
          case 'J':
           if ( g_compare_userset == APPF )
             {
              fprintf(stderr, "Options 'j' and 'k' cannot be used together\n\n");
              print_help();
              mtfexit(-1);
             }
           g_compare_userset = PIV;
           strcpy(g_spec,"PIV");
           break;
          case 'k':
         case 'K':
           if ( g_compare_userset == PIV )
             {
              fprintf(stderr, "Options 'j' and 'k' cannot be used together\n\n");
              print_help();
              mtfexit(-1);
             }
           g_compare_userset = APPF;
           strcpy(g_spec,"Appendix F");
           break;
          case 'n':
          case 'N':
           /*  g_nocompare_userset = 1; */
           fprintf (stderr,"'n' option is no longer needed.\n");
            break;
          case 'p':
          case 'P':
            g_printer = 1;
            break;
          case 'r':
          case 'R':
            print_data_rights_notice();
           print_help();
           break;
          case 'f':
          case 'F':
            /* Force image polarity to be reversed */
            g_reversepolarity = 1;
            break;
          case ' ':
          case '\n':
          case '\0':
            break;
          default:
            fprintf(stderr," INVALID OPTION INPUT\n");
            mtfexit(1);
            break;
        }
    } while((count < 80) && (switch_entry != '\n') && (switch_entry != '\0')
            && (read_value == 1));
 
  if( g_printer && g_digital_camera )
  {
       fprintf(stderr, "\nERROR: Options 'p' and 'g' are incompatible. Exiting\n\n");
       mtfexit(1);
  }
  if( g_user_requestcurve == 1 ) {
    int temp;
    do {
      fprintf(stderr," \nSelect number for desired input/output curve fit type: \n");
      fprintf(stderr,"  5    Piecewise Linear fit - usually last resort \n");
      fprintf(stderr,"  6    LinReg - best fit, single straight line over all input points \n");
      fprintf(stderr,"  7    Unsmoothed spline - local cubics go through all input points \n");
      fprintf(stderr,"  8    Smoothed spline - input points smoothed before applying local cubics \n");
      fprintf(stderr,"  9    S-curve fit - input points fit to a logistic growth curve \n\n");
      fprintf(stderr,"       Note: if spline fit encounters problem, it reverts to S-curve fit \n");
      scanf("%d",&temp);
      g_user_requestcurve = temp;
    } while(g_user_requestcurve < 5 && g_user_requestcurve > 9);
  }
}
 
 
 
 
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
 
void get_args(char *image_filename,char *data_filename, TIFF** tif,
             unsigned int *white_patchp, unsigned int *black_patchp,
            unsigned int *lowf_patchp)
{
 
  char orientation[10];
  int read_success;
  int angle;
 
  while ((*tif == NULL) && (g_scan_image_file_id <= 0))
    {
      if (!g_streamline_args) {
        fprintf(stderr,"\n\n Enter image filename:   \t");
        scanf("%s",image_filename);
      }
 
#ifdef USE_TIFF
        {
            char *tmp;
            if(  ((tmp = strstr (image_filename, ".pgm")) == NULL ||
                       tmp != strrchr (image_filename, '.')) &&
                     ((tmp = strstr (image_filename, ".PGM")) == NULL ||
                       tmp != strrchr (image_filename, '.') ) )
                *tif = TIFFOpen((char *)image_filename, "r");
        }
#endif
      if (*tif == NULL ) {
 
       g_pgm = 0;
       read_pgm_header (image_filename);
 
       g_scan_image_file_id = open(image_filename, (O_RDONLY | O_BINARY));
        if (g_scan_image_file_id < 0) {
          fprintf(stderr, "\nCould not open %s\n", image_filename);
          fprintf(stderr, "%s\n\n", strerror(errno));
        }
      }
    }
 
  if (*tif == NULL ) {
    if ( !g_pgm )
     {
      fprintf(stderr,"Image file is not in TIFF or PGM format\n");
      fprintf(stderr,"Enter image header bytes, pixel width, pixel height:\t");
      GET_COORD_PAIR(&g_scan_header,&g_scan_cols,&g_scan_rows);
      g_bps=8;
      g_photometric=1;
     }
  }
#ifdef USE_TIFF
  else {
    unsigned short sampleformat;
      unsigned short compressed;
    g_bps = 8;
    sampleformat = 1;
    TIFFGetField( *tif, TIFFTAG_IMAGEWIDTH, &g_scan_cols );
    (void) TIFFGetField( *tif, TIFFTAG_IMAGELENGTH, &g_scan_rows );
    (void) TIFFGetField( *tif, TIFFTAG_BITSPERSAMPLE, &g_bps );
 
    /* reading in the tag on compression. Don't process lossy or b/w formats. */
    (void) TIFFGetField( *tif, TIFFTAG_COMPRESSION, &compressed);
    if ((compressed != COMPRESSION_NONE) &&
        (compressed != COMPRESSION_LZW) &&
        (compressed != COMPRESSION_ADOBE_DEFLATE) &&
        (compressed != COMPRESSION_PACKBITS)) {
        MTFPRINT2("ERROR: B/W or Lossy Compressed Data. Compression type = %d\n", compressed);
        mtfexit(-1);
    }
 
    /* reading in the tag on photometric interpretation - 7/22/94 - gtk */
    (void) TIFFGetField( *tif, TIFFTAG_PHOTOMETRIC, &g_photometric);
    /*    MTFPRINT2("G_photometric = %d\n", g_photometric); */
 
    (void) TIFFGetField( *tif, TIFFTAG_SAMPLEFORMAT, &sampleformat);
    if (sampleformat != 1) {
      MTFPRINT2("ERROR: Signed or float imagery. SampleFormat = %d\n", sampleformat);
      mtfexit(-1);
    }
    if( g_bps != 8 && g_bps != 16){
      MTFPRINT2("ERROR: Cannot handle %d BitsPerPixel TIFF format. Try 8 or 16.\n",g_bps);
      mtfexit(-1);
    }
  }
#endif
  g_bytes_per_pixel = (g_bps+7)/8;
  g_max_pixel_value = (1<<g_bps)-1;
 
  /*
      MTFPRINT2("G_bps = %hd\n", g_bps);
      MTFPRINT2("G_bytes_per_pixel = %d\n", g_bytes_per_pixel);
      MTFPRINT2("G_max_pixel_value = %d\n", g_max_pixel_value);
  */
 
  read_success = 0;
  while (read_success == 0)
    {
      if (!g_streamline_args) {
          fprintf(stderr, "Enter target data filename:   \t");
          scanf("%s", data_filename);
      }
 
      /* Read the target data file immediately, so that g_bartarget is set v5.5 */
      read_success = get_target_data(data_filename, white_patchp, black_patchp, lowf_patchp);
    }
 
   if (!g_streamline_args) {
    fprintf(stderr,"\n     (ref: col,row = 0,0 at upper left corner of ENTIRE image,\n");
    fprintf(stderr,"     cols increase left to right, rows increase top to bottom.)\n");
 
    if (g_bartarget) {
      fprintf(stderr,"\n Enter Col, Row for each of 3 target corners:\n");
 
      fprintf(stderr,
             "  upper left corner:\t");
 
      GET_FLOAT_PAIR(&g_input_test_pattern_xul,&g_input_test_pattern_yul)
 
      fprintf(stderr, "  upper right corner:\t");
      GET_FLOAT_PAIR(&g_input_test_pattern_xur,&g_input_test_pattern_yur)
 
      fprintf(stderr, "  lower left corner:\t");
 
      GET_FLOAT_PAIR(&g_input_test_pattern_xll,&g_input_test_pattern_yll)
 
      fprintf(stderr," \n");
      fprintf(stderr," \n");
      fprintf(stderr,"        | | | | | \n");
      fprintf(stderr,"        | | | | | \n");
      fprintf(stderr,"        | | | | | \n");
      fprintf(stderr,"        | | | | |                  H\n");
      fprintf(stderr,"        | | | | | \n");
      fprintf(stderr,"        | | | | | \n");
      fprintf(stderr,"        | | | | | \n");
      fprintf(stderr," \n");
      fprintf(stderr," \n");
      fprintf(stderr,"     ________________ \n");
      fprintf(stderr,"     ________________ \n");
      fprintf(stderr,"     ________________              V\n");
      fprintf(stderr,"     ________________ \n");
      fprintf(stderr,"     ________________ \n");
      fprintf(stderr," \n");
      fprintf(stderr,"Enter bar measurement direction,  H  or  V : \t");
 
    }
    else {
      fprintf(stderr,"\n Enter Col, Row for each of 3 gray patch corners:\n");
 
      fprintf(stderr,
             "  lower right corner of the UPPER LEFT PATCH:\t");
 
      GET_FLOAT_PAIR(&g_input_test_pattern_xul,&g_input_test_pattern_yul)
 
      fprintf(stderr, "  lower left corner of UPPER RIGHT PATCH:\t");
      GET_FLOAT_PAIR(&g_input_test_pattern_xur,&g_input_test_pattern_yur)
 
      fprintf(stderr, "  upper right corner of LOWER LEFT PATCH:\t");
 
      GET_FLOAT_PAIR(&g_input_test_pattern_xll,&g_input_test_pattern_yll)
 
      fprintf(stderr," \n");
      fprintf(stderr," \n");
      fprintf(stderr," |||| |  |  | | |||| |  |  | | |||| \n");
      fprintf(stderr," |||| |  |  | | |||| |  |  | | |||| \n");
      fprintf(stderr," |||| |  |  | | |||| |  |  | | ||||        H\n");
      fprintf(stderr," |||| |  |  | | |||| |  |  | | |||| \n");
      fprintf(stderr," |||| |  |  | | |||| |  |  | | |||| \n");
      fprintf(stderr," |||| |  |  | | |||| |  |  | | |||| \n");
      fprintf(stderr," \n");
      fprintf(stderr," \n");
      fprintf(stderr,"     ______________________ \n");
      fprintf(stderr,"     ______________________ \n");
      fprintf(stderr,"     ______________________ \n");
      fprintf(stderr,"     ______________________ \n");
      fprintf(stderr," \n");
      fprintf(stderr,"     ______________________ \n");
      fprintf(stderr," \n");
      fprintf(stderr," \n");
      fprintf(stderr,"     ______________________              V\n");
      fprintf(stderr," \n");
      fprintf(stderr," \n");
      fprintf(stderr,"     ______________________ \n");
      fprintf(stderr," \n");
      fprintf(stderr,"     ______________________ \n");
      fprintf(stderr,"     ______________________ \n");
      fprintf(stderr,"     ______________________ \n");
      fprintf(stderr,"     ______________________ \n");
      fprintf(stderr," \n");
      fprintf(stderr,"Enter direction of sine wave,  H  or  V : \t");
    }
 
    scanf("%s", orientation);
 
    if (!strcmp("H",orientation) || !strcmp("h",orientation))
      g_scan_orient = LAND;
    else
      g_scan_orient = PORT;
 
    if (g_bartarget) {
 
       fprintf(stderr, "\n");
       print_upright_info(data_filename);
       fprintf(stderr,"Enter target orientation (0,9,-9,1): \n");
       fprintf(stderr,"0 \tupright target\n");
       fprintf(stderr,"9 \trotated 90 deg clockwise\n");
       fprintf(stderr,"-9 \trotated 90 deg counterclockwise\n");
       fprintf(stderr,"1 \trotated 180 deg\n");
 
       scanf("%d", &angle);
       while (angle != 0 && angle != -9 && angle != 9 && angle != 1) {
          fprintf(stderr,"Enter target orientation (0,9,-9,1): \t");
         scanf("%d", &angle);
       }
       angle *= 10;
       if (angle == 10) angle = 180;
       g_upright = angle;
 
    }
  }
 
  fprintf(stderr," \n");
  if (g_printer)
    {
      /* Get the printer ppi */
 
       {
        fprintf(stderr, "printer ppi = \n");
        fprintf(stderr, "    (digital target width in pixels) / (printed target width in inches) \n");
        fprintf(stderr,"Enter printer ppi: ");
        scanf("%u", &g_printer_ppi);
        fprintf(stderr," \n");
      }
    }
 
  if ( g_digital_camera )
  {
    float temp1, temp2;
    /* Get the focal length and object distance if known, else the scale */
    char string[25];
       int num;
 
    fprintf(stderr, "Enter: Lens FL (mm)  Target Distance (inches)\n");
       fprintf(stderr, "       or\n       Scale Factor (TargetWidth/ImageWidth)\n");
       /* Read in any intermediate linefeeds */
       do     fgets(string, 20, stdin);
       while (strncmp(string,"\n",1) == 0);
 
       num = sscanf(string, "%f%f", &temp1, &temp2);
       if (num > 1)
       {
         g_obj_img_scale = (MM_PER_INCH*temp2 - temp1)/(double)temp1;
         g_FL = temp1;
         g_target_distance = temp2;
       }
       else
       {
         g_obj_img_scale = temp1;
         g_FL = 0.0;
       }
       fprintf(stderr, " \n");
  }
  else g_obj_img_scale = 1.0;
 
}
 
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
 
 
void print_development_notice(void)
{
  fprintf(stderr,"\n\n\n\n");
  fprintf(stderr,"********************************************************\n");
  fprintf(stderr,"*                                                      *\n");
  fprintf(stderr,"*                                                      *\n");
  fprintf(stderr,"*                                                      *\n");
  fprintf(stderr,"*                                                      *\n");
  fprintf(stderr,"*                     sineMTF                          *\n");
  fprintf(stderr,"*                                                      *\n");
  fprintf(stderr,"*                     V%s\t                       *\n", VERSION);
  fprintf(stderr,"*                                                      *\n");
  fprintf(stderr,"*         Developed by The MITRE Corporation           *\n");
  fprintf(stderr,"*                                                      *\n");
  fprintf(stderr,"*                                                      *\n");
  fprintf(stderr,"*                                                      *\n");
  fprintf(stderr,"*                                                      *\n");
  fprintf(stderr,"********************************************************\n");
 
  fprintf(stderr,"\n\n\n\n\n\n");
 
 
 
  /* formfeeds don't seem to work on a vt terminal
     fprintf(stderr,"\f");
  */
}
 
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
 
void print_data_rights_notice(void)
{
  fprintf(stderr, "\n\n");
  fprintf(stderr, "*******************************************************************************\n");
  fprintf(stderr, "                               MITRE NOTICE\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "This \"sinemtf\" software was developed by The MITRE Corporation and was produced\n");
  fprintf(stderr, "for the U.S. Government under Contract numbers J-FBI-12-128, J-FBI-86-028, \n");
  fprintf(stderr, "J-FBI-93-039, DAAB07-99-C-C201,DAAB07-00-C-C201, DAAB07-01-C-C201, \n");
  fprintf(stderr, "DAAB07-01-C-N200, DAAB07-02-C-N200, W15P7T-04-C-D001, and W15P7T--5-C-D199 \n");
  fprintf(stderr, "and is subject to the Rights in Data-General clause FAR 52.227-l4, \n");
  fprintf(stderr, "Alt. IV (DEC 2007). \n");
  fprintf(stderr," \n");
  fprintf(stderr, "Redistribution and use in source and binary forms, with or without \n");
  fprintf(stderr, "modification, are permitted provided that the following conditions are met:\n");
  fprintf(stderr, "-Redistribution of source code must retain the above copyright notice, this \n");
  fprintf(stderr, "list of conditions, and the following disclaimer.\n");
  fprintf(stderr, "-Redistribution in binary form must reproduce the above copyright notice, this \n");
  fprintf(stderr, "list of codition and the following disclaimer in the documentation and/or other\n");
  fprintf(stderr, "materials provided with the distribution.\n");
  fprintf(stderr," \n");
  fprintf(stderr, "MITRE IS PROVIDING THIS SOFTWARE \"AS IS\" AND MAKES NO WARRANTY, EXPRESS OR \n");
  fprintf(stderr, "IMPLIED, AS TO THE ACCURACY, CAPABILITY, EFFICIENCY, MERCHANTABILITY, OR\n");
  fprintf(stderr, "FUNCTIONING OF THIS SOFTWARE AND/OR DOCUMENTATION.  IN NO EVENT WILL MITRE BE\n");
  fprintf(stderr, "LIABLE FOR ANY GENERAL, CONSEQUENTIAL, INDIRECT, INCIDENTAL, EXEMPLARY, OR\n");
  fprintf(stderr, "SPECIAL DAMAGES, EVEN IF MITRE HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH\n");
  fprintf(stderr, "DAMAGES, IRRESPECTIVE OF THE CAUSE OF SUCH DAMAGES.\n");
  fprintf(stderr," \n");
  fprintf(stderr,"*******************************************************************************\n");
  fprintf(stderr, "                               PROGRAM DEVELOPMENT\n");
  fprintf(stderr,"Software Development: Steve Bennett, Dave Braunegg, Wendy Etkind,\n");
  fprintf(stderr,"                      Dave Gonthier, Greg Kearnan, Margaret Lepley, \n");
  fprintf(stderr,"                      Norm Nill, Barry Paine, John White\n");
  fprintf(stderr,"Concept/Design, QA, Project Leader: Norm Nill\n");
  fprintf(stderr," \n");
  fprintf(stderr,"*******************************************************************************\n");
  fprintf(stderr," \n");
  /* mtfexit(0); */
}
 
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
 
void print_help(void)
{
  fprintf(stderr,"*************** HELP ***************\n");
  fprintf(stderr," \n");
  fprintf(stderr,"a \tOutput average MTF (in addition to peak MTF)\n");
  fprintf(stderr,"    Average MTF - average of all sample modulations, per frequency\n");
  fprintf(stderr,"    Peak MTF - selects highest modulation sample, per frequency\n");
  fprintf(stderr,"    A sample modulation is computed from an adjacent peak/valley pair in\n");
  fprintf(stderr,"    target reflectance space.\n");
  fprintf(stderr," \n");
  fprintf(stderr,"b \tAuto refine input corner coordinates\n");
  fprintf(stderr,"    Assumes user input corner coords are apprx., refines them to exact corner \n");
  fprintf(stderr,"    coords.\n");
  fprintf(stderr," \n");
  fprintf(stderr,"c \tSelect I/O curve fit type for sine tgt (instead of autofit)\n");
  fprintf(stderr,"    Options for curve fit (target reflectance vs. image gray):\n");
  fprintf(stderr,"    Piecewise linear fit - different straight line between every 2 adjacent \n");
  fprintf(stderr,"                           data points;\n");
  fprintf(stderr,"                           (usually a last resort) \n");
  fprintf(stderr,"    Linear regression - single straight line over all data points.\n");
  fprintf(stderr,"    Smoothed spline - local cubic splines through smoothed data points.\n");
  fprintf(stderr,"    Unsmoothed spline - local cubic splines through all data points; \n");
  fprintf(stderr,"                        if poor fit, auto reverts to S-curve fit.\n");
  fprintf(stderr,"    S-curve - nonlinear regression fit to a logistics S shaped curve.\n");
  fprintf(stderr,"    If program selects curve (no option c), then: \n");
  fprintf(stderr,"       If Not a printer, best curve is fit, in chi-squared sense, from above \n");
  fprintf(stderr,"       options (but excluding piecewise linear fit) \n");
  fprintf(stderr,"       If printer, S curve is fit. \n");
  fprintf(stderr," \n");
#ifdef USE_TIFF
  fprintf(stderr,"d \tCreate diagnostic image (imagename_box.tif)\n");
#else
  fprintf(stderr,"d \tCreate diagnostic image (imagename_box.pgm)\n");
#endif
  fprintf(stderr,"    Creates same size copy of original image in raw format, overlayed with \n");
  fprintf(stderr,"    measurement locations; verifies dimensional data in InputDataFile and \n");
  fprintf(stderr,"    auto corner coords, if applied.\n");
  fprintf(stderr," \n");
  fprintf(stderr,"e \tVerbose output\n");
  fprintf(stderr,"    Prints statistics (sine & density patches), sample modulations, & tabulates\n");
  fprintf(stderr,"    I/O curve.\n");
  fprintf(stderr," \n");
  fprintf(stderr,"f \tReverse image polarity\n");
  fprintf(stderr,"    Reverse image polarity before processing; occasionally needed so program \n");
  fprintf(stderr,"    reads pure white as gray=255, not gray=0 (original image file not changed).\n");
  fprintf(stderr," \n");
  fprintf(stderr,"g \tDigital camera assessment\n");
  fprintf(stderr,"    User inputs focal length & distance-to-target, or, independently computed \n");
  fprintf(stderr,"    scale (scale = target width / width at sensor image plane).\n");
  fprintf(stderr,"    Use digital image size as output by camera (not enlarged or reduced).\n");
  fprintf(stderr,"    A color image must be converted to grayscale before input to program.\n");
  fprintf(stderr," \n");
  fprintf(stderr,"j \tCompare output to FBI's PIV spec\n");
  fprintf(stderr,"    option affects printout content, not actual computations.\n");
  fprintf(stderr," \n");
  fprintf(stderr,"k \tCompare output to FBI's App. F spec\n");
  fprintf(stderr,"    option affects printout content, not actual computations.\n");
  fprintf(stderr," \n");
  fprintf(stderr,"p \tPrinter MTF assessment\n");
  fprintf(stderr,"    Assumes that: \n");
  fprintf(stderr,"     supplied digital sine target (A6/7sine.tif) is printed to hardcopy,\n");
  fprintf(stderr,"     hardcopy print is  scanned via scanner to produce digital image, \n");
  fprintf(stderr,"     digital image is input to sinemtf program, \n");
  fprintf(stderr,"     MTF output of sinemtf is divided by scanner MTF, resulting in printer MTF.\n");
  fprintf(stderr," \n");
  fprintf(stderr,"Bar Target Processing\n");
  fprintf(stderr," Input data file:\n");
  fprintf(stderr,"   'b' replaces 's' in 2nd line of each frequency data set\n");
  fprintf(stderr,"   no density patches in data file \n");
  fprintf(stderr,"   If system I/O is nonlinear, add nonlinearity table.\n");
  fprintf(stderr,"   After last bar pattern info, put {input_val output_grey} pairs\n");
  fprintf(stderr,"   per line. Linear interpolation will be used between. E.g.\n");
  fprintf(stderr,"\t 0.1 \t 5.0\n\t 0.12 \t 7.\n\t 0.16 \t 9.\n\t . \t .\n\t . \t .\n\t . \t .\n\t 1.0 \t 255.\n");
  fprintf(stderr," Locating corner points of target image:\n");
  fprintf(stderr,"   fit an imaginary rectangle to the displayed image, such that each of the \n");
  fprintf(stderr,"   rectangle's 4 sides coincides with the outer edge of the outermost bar \n");
  fprintf(stderr,"   along each side, then UL,UR,LL corners of rectangle are the coordinates\n");
  fprintf(stderr," Normalization: \n");
  fprintf(stderr,"   divide program output modulations by W to get system CTF, \n");
  fprintf(stderr,"   W = modulation of very low freq. bar (< 2.5%% of Nyquist freq.) \n");
  fprintf(stderr," \n");
  fprintf(stderr," \n");
  fprintf(stderr,"More detail in: \n");
  fprintf(stderr,"  MTF_Guide.pdf\n");
  fprintf(stderr,"  MTF_TechRpt.pdf\n");
  fprintf(stderr,"  BarTgt_Guide.pdf\n");
  fprintf(stderr," \n");
}
 
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
 
void print_header(char *image_filename,char *data_filename,
                  unsigned char original_scan_orient)
{
  time_t tm_s,*tm_ptr;
  char tod[82];
  double fabs(double);
 
  tm_ptr = &tm_s;
 
  if (time(tm_ptr) == -1)
    strcpy(tod,"DATE = ????");
  else
    strcpy(tod,ctime(tm_ptr));
 
  MTFPRINT3("%.24s\tv%s\n",tod,VERSION)
  MTFPRINT2("Target Serial Number is:\t%s\n",g_serial_number)
  MTFPRINT2("Image input file name = \t%s\n",image_filename)
  MTFPRINT2("Data input file name = \t\t%s\n",data_filename)
    if (g_bartarget)
  MTFPRINT("\nInput coordinates of target image:")
    else
  MTFPRINT("\nInput coordinates of inner corners of target image:")
    if(g_refine_corners && !g_bartarget)
  MTFPRINT("(auto computed)\n")
    else
  MTFPRINT("\n")
  MTFPRINT3("Upper Left Col =\t\t%g\tUpper Left Row =\t%g",
            g_input_test_pattern_xul,g_input_test_pattern_yul)
  if (g_UL_auto_crnr_err)
    MTFPRINT(" [*]")
  MTFPRINT("\n")
  MTFPRINT3("Upper Right Col =\t\t%g\tUpper Right Row =\t%g",
            g_input_test_pattern_xur,g_input_test_pattern_yur)
  if (g_UR_auto_crnr_err)
    MTFPRINT(" [*]")
  MTFPRINT("\n")
  MTFPRINT3("Lower Left Col =\t\t%g\tLower Left Row =\t%g",
            g_input_test_pattern_xll,g_input_test_pattern_yll)
  if (g_LL_auto_crnr_err)
    MTFPRINT(" [*]")
  MTFPRINT("\n")
  if (g_LL_auto_crnr_err || g_UR_auto_crnr_err || g_UL_auto_crnr_err)
    MTFPRINT(" [*] auto corner refinement failed for this corner, orig data used instead\n")
 
/*
  MTFPRINT("\nCalculated coordinates of outer corners of target image:\n")
  MTFPRINT3("Upper Left Col  =\t\t%.1f\tUpper Left Row  =\t%.1f\n",
            g_orig_test_pattern_xul,g_orig_test_pattern_yul)
  MTFPRINT3("Upper Right Col =\t\t%.1f\tUpper Right Row =\t%.1f\n",
            g_orig_test_pattern_xur,g_orig_test_pattern_yur)
  MTFPRINT3("Lower Left Col  =\t\t%.1f\tLower Left Row  =\t%.1f\n",
            g_orig_test_pattern_xll,g_orig_test_pattern_yll)
  MTFPRINT3("Lower Right Col =\t\t%.1f\tLower Right Row =\t%.1f\n",
            g_orig_test_pattern_xlr,g_orig_test_pattern_ylr)
*/
 
  if (g_scan_image_file_id == 0)
    MTFPRINT3("\nInput image columns = %d, rows = %d\n",
              g_scan_cols, g_scan_rows)
  else
    MTFPRINT4("\nHdr, Cols, Rows of input image = %d, %d, %d\n",
              g_scan_header, g_scan_cols, g_scan_rows)
 
  if (g_bartarget)
    MTFPRINT3("MTF direction = \t%s \tTarget Orientation \t%d deg\n",
            (original_scan_orient == LAND) ? "Horizontal" : "Vertical",
           g_upright)
  else
    MTFPRINT2("MTF direction = \t%s\n",
             (original_scan_orient == LAND) ? "Horizontal" : "Vertical")
 
  if (g_printer)
    MTFPRINT2("Printer resolution (pixels per inch):\t%d (user input)\n", g_printer_ppi);
 
  if (g_digital_camera && g_FL==0.0)
    MTFPRINT2("Target-to-image scale: %.1f ()\n",g_obj_img_scale)
  if (g_digital_camera && g_FL!=0.0)
    MTFPRINT4("Target-to-image scale: %.1f (FL=%.2f mm   Distance=%.1f inches)\n",
             g_obj_img_scale, g_FL, g_target_distance)
  if (original_scan_orient == LAND){
    if (!g_digital_camera)
       MTFPRINT3("Scanner resolution (pixels per inch):\tHORZ = %.3f\tVERT = %.3f\n",
              g_ppi_horz, g_ppi_vert)
    MTFPRINT4("Target Skew Angle:  HORZ = %.3f  VERT = %.3f   |avg.skew| = %.3f degrees\n",
              g_skew_horz*(double)(180.0/PI),
              g_skew_vert*(double)(180.0/PI),
              (fabs(g_skew_avg))*(double)(180.0/PI))
  }
  else{
       if (!g_digital_camera)
       MTFPRINT3("Scanner resolution (pixels per inch):\tHORZ = %.3f\tVERT = %.3f\n",
              g_ppi_vert, g_ppi_horz)
    MTFPRINT4("Target Skew Angle:  HORZ = %.3f  VERT = %.3f   |avg.skew| = %.3f degrees\n",
              g_skew_vert*(double)(180.0/PI),
              g_skew_horz*(double)(180.0/PI),
              (fabs(g_skew_avg))*(double)(180.0/PI))
  }
 
  if(g_bartarget && g_nonlinear) {
    MTFPRINT("NonLinear  I/O curve applied\n")
  }
  MTFPRINT("\n")
 
}
 
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
/* Prints problems depending upon type.  Problems marked as IQS spec related
   are not printed if user has requested no IQS report.  Aliasing problems
   however are always reported.  v5.5 */
 
void print_problems(void)
{
  short i;
 
  if (g_problem_count == 0) {
    if(g_compare)
      MTFPRINT2("\n%s spec check applied.",g_spec)
    else
      MTFPRINT("\nNo spec check applied.")
   MTFPRINT("\nNo problems encountered.\n\n")
  }
  else if ((!g_compare) && g_problem_count-g_IQS_problem_count == 0)
  {
       /* Do nothing when only IQS problems present, but IQS report is
       supressed. */
  }
  else {
    i=0;
 
    /* note, "\f" causes a new page */
    MTFPRINT("\n\t\t\tPROBLEM REPORT\n\n")
 
    if(g_compare)
      MTFPRINT2("%s spec check applied.\n",g_spec)
    else
      MTFPRINT("No spec check applied.\n")
 
    /* Add in footnote about scanner spec, if needed. */
    REPORT_SPEC_FOOTNOTE
 
    while (g_problem_count--) {
      /* Print when IQS comparison is requested *or* when non-IQS problem */
      if ((g_compare) || g_problem_type[i] != IQS)
          MTFPRINT2("%s",g_problems[i])
      i++;
    }
 
    if (too_many_problems != 0)
      {
        MTFPRINT("\nMore problems were found, but number exceeds program limits.\n")
        too_many_problems = 0;
      }
  }
}
 
 
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
 
void put_problem(char *s, int problem_type)
{
 
  if (g_problem_count < MAX_PROBLEMS) {
    g_problem_type[g_problem_count] = (short)problem_type;
    strcpy(&g_problems[g_problem_count++][0],s);
  }
  else
    {
      fprintf(stderr, "Number of problems exceeds count of %d\n", MAX_PROBLEMS);
      too_many_problems = 1;
    }
  if (problem_type == IQS)
    g_IQS_problem_count++;
 
}
 
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
 
 
int read_pgm_header(char *image_filename)
{
  FILE *fp;
  int val, rows, cols, numread;
  char header[256];
 
  fp = fopen(image_filename, "rb");
  if (fp == NULL) return (-2);
 
  do {
    if (fgets(header,255,fp) == NULL)
      return (-1);
  } while (header[0] == '#');
 
  if ( strncmp(header, "P5", 2) != 0 ) return (-1);
 
  do {
    if (fgets(header,255,fp) == NULL)
      return (-1);
  } while (header[0] == '#');
 
  numread = sscanf(header,"%d %d", &cols, &rows);
  g_scan_cols = cols;
  if (numread != 2) {
    /*Only one number on this line, so read the next line to get height. */
    fgets(header,255,fp);
    sscanf(header,"%d", &rows);
  }
  g_scan_rows = rows;
 
  fgets(header,255,fp);
  sscanf(header,"%d", &val);
  if(val <= 255)
    g_bps = 8;
 
  g_scan_header = ftell(fp);
 
  fclose(fp);
  g_pgm = 1;
  g_photometric = 1;
 
  return(0);
}
 
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
 
 
static void write_debug_image(short rotation)
{
  int i, j;
  unsigned int ui, uj, count;
  unsigned char *buf;
  unsigned int width, height;
#ifdef USE_TIFF
  unsigned short * red_table;
  unsigned short * green_table;
  unsigned short * blue_table;
#else
  unsigned char * bmpimg; 
  int imgsize, padsize, sizeall;
  unsigned char bmpheader[14]={0};
  unsigned char bmpinfo[40]={0};
  unsigned char bmppad[3]={0};
#endif

 
  if ((buf = (unsigned char *)calloc((size_t)(g_scan_cols), sizeof(char))) == NULL)
    {
      fprintf(stderr, "\nCannot allocate memory to save diagnostic output image: %s\n\n", g_debug_file_name);
    }
  else
    {
      count = 0;
 
#ifdef USE_TIFF
 
  if((red_table=(unsigned short*)malloc(NUM_COLORS*sizeof(unsigned short)))==NULL || (blue_table=(unsigned short*)malloc(NUM_COLORS*sizeof(unsigned short)))==NULL || (green_table=(unsigned short*)malloc(NUM_COLORS*sizeof(unsigned short)))==NULL ){
	fprintf(stderr,"\nCannot allocate memory for color LUTs\n\n");
  }
      fprintf(stderr, "\nWriting out diagnostic image (TIF): %s      ", g_debug_file_name);

#else
	if((bmpimg=(unsigned char *)malloc(3*g_total_image_width*g_total_image_height*sizeof(unsigned char)))==NULL){
	fprintf(stderr,"\nCannot allocate memory for rgb pixel arrays");
}

      	fprintf(stderr, "\nWriting out diagnostic image (BMP): %s      ", g_debug_file_name);
#endif

      fflush(stderr);
 
      if(rotation==RIGHT || rotation==LEFT)
      {
        width = g_total_image_width;
        height = g_total_image_height;
      }
      else
      {
      width = g_total_image_height;
      height = g_total_image_width;
      }
 
	write_debug_bounds(rotation);

#ifdef USE_TIFF
      TIFFSetField(dtif, TIFFTAG_IMAGEWIDTH, width);
      TIFFSetField(dtif, TIFFTAG_IMAGELENGTH, height);
      TIFFSetField(dtif, TIFFTAG_SAMPLESPERPIXEL, 1);
      TIFFSetField(dtif, TIFFTAG_BITSPERSAMPLE, g_bps);
      TIFFSetField(dtif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
      TIFFSetField(dtif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
      TIFFSetField(dtif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_PALETTE);
      TIFFSetField(dtif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
      i = 8192 / TIFFScanlineSize(dtif);
      TIFFSetField(dtif, TIFFTAG_ROWSPERSTRIP, i==0?1:i);
 
    for(i=0;i<NUM_COLORS;i++){
            if(i==COLORB){
                red_table[i]=((unsigned short)COLORB_RED)<<g_bps;
                green_table[i]=((unsigned short)COLORB_GREEN)<<g_bps;
                blue_table[i]=((unsigned short)COLORB_BLUE)<<g_bps;
            }
            else if(i==COLORA){
                red_table[i]=((unsigned short)COLORA_RED)<<g_bps;
                green_table[i]=((unsigned short)COLORA_GREEN)<<g_bps;
                blue_table[i]=((unsigned short)COLORA_BLUE)<<g_bps;
            }
	    else if(i==COLORC){
		red_table[i]=((unsigned short)COLORC_RED)<<g_bps;
		green_table[i]=((unsigned short)COLORC_GREEN)<<g_bps;
		blue_table[i]=((unsigned short)COLORC_BLUE)<<g_bps;
	    }
            else{
                red_table[i]=green_table[i]=blue_table[i]=(unsigned short)i<<g_bps;
            }
    }
    TIFFSetField(dtif,TIFFTAG_COLORMAP,red_table,green_table,blue_table);
    free(red_table);
    free(green_table);
    free(blue_table);
#else
      /*fprintf(g_ofp,"P5\n%d %d\n255\n",width,height);*/

	/*adapted from code found on stackoverflow*/
	/*byte order is little endian*/
	
	for(i=0;i<g_total_image_width;i++){
		for(j=0;j<g_total_image_height;j++){
			if(g_image_array_modified[i+(j*g_total_image_width)]==COLORA){
				bmpimg[(i+j*g_total_image_width)*3+2]=(unsigned char)COLORA_RED;
				bmpimg[(i+j*g_total_image_width)*3+1]=(unsigned char)COLORA_GREEN;
				bmpimg[(i+j*g_total_image_width)*3]=(unsigned char)COLORA_BLUE; 
			} 
			else if(g_image_array_modified[i+(j*g_total_image_width)]==COLORB){
				bmpimg[(i+j*g_total_image_width)*3+2]=(unsigned char)COLORB_RED;
				bmpimg[(i+j*g_total_image_width)*3+1]=(unsigned char)COLORB_GREEN;
				bmpimg[(i+j*g_total_image_width)*3]=(unsigned char)COLORB_BLUE;
			}
			else if(g_image_array_modified[i+(j*g_total_image_width)]==COLORC){
				bmpimg[(i+j*g_total_image_width)*3+2]=(unsigned char)COLORC_RED;
				bmpimg[(i+j*g_total_image_width)*3+1]=(unsigned char)COLORC_GREEN;
				bmpimg[(i+j*g_total_image_width)*3]=(unsigned char)COLORC_BLUE;
			}
			else{
				bmpimg[(i+j*g_total_image_width)*3+2]=bmpimg[(i+j*g_total_image_width)*3+1]=bmpimg[(i+j*g_total_image_width)*3]=g_image_array_modified[i+(j*g_total_image_width)];
			}
		}
	}

	bmpheader[0]='B';
	bmpheader[1]='M';
	bmpheader[10]=(unsigned char)54;

	bmpinfo[0]=(unsigned char)40;
	bmpinfo[12]=(unsigned char)1;
	bmpinfo[14]=(unsigned char)24;
	/*change bmp tags based on orientation (since width, height, and padsize change accordingly)*/
	if(rotation==RIGHT || rotation==LEFT){
		padsize=(4-(g_total_image_width)%4)%4;
		imgsize=g_total_image_width*g_total_image_height+padsize*3+g_total_image_height*padsize;
		sizeall=imgsize+sizeof(bmpheader)+sizeof(bmpinfo);

		bmpheader[2]=(unsigned char)(sizeall);
		bmpheader[3]=(unsigned char)(sizeall>>8);
		bmpheader[4]=(unsigned char)(sizeall>>16);
		bmpheader[5]=(unsigned char)(sizeall>>24);
	
		bmpinfo[4]=(unsigned char)(g_total_image_width);
		bmpinfo[5]=(unsigned char)(g_total_image_width>>8);
		bmpinfo[6]=(unsigned char)(g_total_image_width>>16);
		bmpinfo[7]=(unsigned char)(g_total_image_width>>24);
		bmpinfo[8]=(unsigned char)(g_total_image_height);
		bmpinfo[9]=(unsigned char)(g_total_image_height>>8);
		bmpinfo[10]=(unsigned char)(g_total_image_height>>16);
		bmpinfo[11]=(unsigned char)(g_total_image_height>>24);
		bmpinfo[24]=(unsigned char)(imgsize);
		bmpinfo[25]=(unsigned char)(imgsize>>8);
		bmpinfo[26]=(unsigned char)(imgsize>>16);
		bmpinfo[27]=(unsigned char)(imgsize>>24);
	}
	else if(rotation==TOP || rotation==BOTTOM){
		padsize=(4-(g_total_image_height)%4)%4;
		imgsize=g_total_image_width*g_total_image_height+padsize*3+g_total_image_width*padsize;
		sizeall=imgsize+sizeof(bmpheader)+sizeof(bmpinfo);

		bmpheader[2]=(unsigned char)(sizeall);
		bmpheader[3]=(unsigned char)(sizeall>>8);
		bmpheader[4]=(unsigned char)(sizeall>>16);
		bmpheader[5]=(unsigned char)(sizeall>>24);
	
		bmpinfo[4]=(unsigned char)(g_total_image_height);
		bmpinfo[5]=(unsigned char)(g_total_image_height>>8);
		bmpinfo[6]=(unsigned char)(g_total_image_height>>16);
		bmpinfo[7]=(unsigned char)(g_total_image_height>>24);
		bmpinfo[8]=(unsigned char)(g_total_image_width);
		bmpinfo[9]=(unsigned char)(g_total_image_width>>8);
		bmpinfo[10]=(unsigned char)(g_total_image_width>>16);
		bmpinfo[11]=(unsigned char)(g_total_image_width>>24);
		bmpinfo[24]=(unsigned char)(imgsize);
		bmpinfo[25]=(unsigned char)(imgsize>>8);
		bmpinfo[26]=(unsigned char)(imgsize>>16);
		bmpinfo[27]=(unsigned char)(imgsize>>24);
	}
	fwrite(bmpheader,1,14,g_ofp);
	fwrite(bmpinfo,1,40,g_ofp);
#endif
 
      switch(rotation)
        /* Rotate the modified image back to its input orientation */
        {
        case RIGHT:
          /* No rotation needed, white patch alreay at top right */
#ifdef USE_TIFF
         for(i=0; i<height; i++) {  /* Write each row separately */
           if( TIFFWriteScanline( dtif, g_image_array_modified+count, i, 0L ) < 0 )
              fprintf(stderr, "\nBad data write on imageline %d\n\n",i);
            count += width;
         }
#else
	/*fwrite(g_image_array_modified, 1, g_total_image_height*g_scan_cols, g_ofp);*/
         for(i=0;i<g_total_image_height;i++){
		fwrite(bmpimg+(g_total_image_width*(g_total_image_height-i-1)*3),3,g_total_image_width,g_ofp);
		fwrite(bmppad,1,padsize,g_ofp);
	}
#endif
          break;
        case LEFT:
          /* Rotate to place white patch at bottom left */
#ifdef USE_TIFF
          j = 0;
         for (i = g_total_image_width*g_total_image_height-1;  i >= 0;  i--) {
          buf[count] = g_image_array_modified[i];
           count++;
           if (count == g_scan_cols) {
             if( TIFFWriteScanline( dtif, buf, j, 0L ) < 0 )
                fprintf(stderr, "\nBad data write on imageline %d\n\n",j);
             count = 0; j++;
           }
         }
#else
             /*fwrite(buf, 1, g_scan_cols, g_ofp);*/
	for(i=0;i<g_total_image_height;i++){
			for(j=g_total_image_width-1;j>=0;j--){
				fwrite(bmpimg+((j+g_total_image_width*i)*3),3,1,g_ofp);
			}
			fwrite(bmppad,1,padsize,g_ofp);
	}	
#endif

         break;
        case TOP:
          /* Rotate to place white patch at top left */
#ifdef USE_TIFF
         for (j = g_total_image_width - 1;  j >= 0;  j--)
           for (ui = 0;  ui < g_total_image_height;  ui++)
           {
             buf[count] = g_image_array_modified[ui*g_total_image_width + j];
             count++;
             if (count == g_scan_cols) {

              if( TIFFWriteScanline( dtif, buf, height-1-j, 0L ) < 0 )
                  fprintf(stderr, "\nBad data write on imageline %d\n\n",height-1-j);
              count = 0;
             }
           }
#else
             /*  fwrite(buf, 1, g_scan_cols, g_ofp);*/
	for(i=0;i<g_total_image_width;i++){
		for(j=0;j<g_total_image_height;j++){
			fwrite(bmpimg+((i+g_total_image_width*j))*3,3,1,g_ofp);
		}
		fwrite(bmppad,1,padsize,g_ofp);
	}
#endif
         break;
        case BOTTOM:
          /* Rotate to place white patch at bottom right */
#ifdef USE_TIFF
         for (uj = 0;  uj < g_total_image_width;  uj++)
           for (i = g_total_image_height - 1;  i >= 0;  i--)
           {
             buf[count] = g_image_array_modified[i*g_total_image_width + uj];
             count++;
             if (count == g_scan_cols) {

              if( TIFFWriteScanline( dtif, buf, uj, 0L ) < 0 )
                  fprintf(stderr, "\nBad data write on imageline %d\n\n",uj);
              count = 0;
             }
           }
#else
            /*  fwrite(buf, 1, g_scan_cols, g_ofp);*/
	for(i=g_total_image_width-1;i>=0;i--){
		for(j=g_total_image_height-1;j>=0;j--){
			fwrite(bmpimg+((i+g_total_image_width*j)*3),3,1,g_ofp);
		}
		fwrite(bmppad,1,padsize,g_ofp);
	}
#endif
         break;
      }
 
      /* Cleanup for LEFT, TOP, and BOTTOM, but this shouldn't be needed */
#ifdef USE_TIFF
      TIFFClose(dtif);
#else
      if (count != 0)
        fwrite(buf, 1, count, g_ofp);
 
      if (fclose(g_ofp) != 0)
        fprintf(stderr, "mtf.c: Could not close g_ofp: %s\n\n", strerror(errno));
#endif
 
      fprintf(stderr, "Done\n");
    }
 
    free(buf);
 
}
 
/*****************************************************************/
 
void mtfexit(int val)
{
#ifdef MSDOS
  while ( getchar() != (int)'\n');
  fprintf(stderr, "\nPress any key to exit ...");
  fflush(stderr);
  while (getchar() == EOF);
  fprintf(stderr, "\n");
#endif
 
  exit(val);
}
 
/******************************************************************************/
/*                                                                  */
/******************************************************************************/
 
/* Read in the image and store it in a big array.  We also flip the image
   into a standard, horizontal format with the white patch at the top
   and right. */
 
void read_in_image(TIFF *tif, short rotation, char *image_name)
{
  short input_row;
  short array_column_start=0;
  short array_row_index;
  short array_column_index;
  short array_column_delta=1;
  short array_row_delta=1;
  short array_row_start=0;
 
  unsigned char *buf;
  short i;
  int j,k;
  float tmp;
  short hmo,wmo;
  float tmpx,tmpy;
 
  if (g_debug) {
    /* given an input image file name of the form "base.ext", here we create
       an output file with the form "base8_box", where "base8" was taken from
       the first 8 characters of the input file name - the first 8 characters
       of "base".  If "base" has less than 8 characters, then all characters
       before the "." are used. 9/2/94 - gtk */
 
    i = strlen(image_name);
    j=i;
    while ((image_name[i] != '/')&&(i>=0))
      {
       if (image_name[i] == '.') j=i;
       i--;
      }
    i++;
    j=j-i;
    if (j > 8) j=8;        /* no more than 8 characters */
    for (k=0;k<j;k++)
      g_debug_file_name[k]=image_name[i+k];
    g_debug_file_name[j]='\0';
#ifdef USE_TIFF
    strncat(g_debug_file_name,"_box.tif",8);
    fprintf(stderr,"\n\nDiagnostic output image:  %s\n\n",g_debug_file_name);
    if((dtif = TIFFOpen(g_debug_file_name, "w")) == NULL) {
       fprintf(stderr, "Could not create %s - try running without 'd' option\n", g_debug_file_name);
       mtfexit(-1);
    }
#else
    strncat(g_debug_file_name,"_box.bmp",8);
    fprintf(stderr,"\n\nDiagnostic output image:  %s\n\n",g_debug_file_name);
    g_ofp = fopen(g_debug_file_name, "w");
    if (g_ofp == NULL) {
      fprintf(stderr,"\n%s%s%s\n","Could not create ",g_debug_file_name,
             " - try running without 'd' option \n");
      mtfexit(-1);
    }
#endif
  }
 
  /* enough to hold a single scan line */
  buf = NULL;
  if (g_scan_image_file_id) {
    buf = (unsigned char*) calloc((unsigned int)g_scan_cols,sizeof(unsigned char));
    g_row_mem = g_scan_cols;
  }
#ifdef USE_TIFF
  else {
    buf = (unsigned char*) calloc((unsigned int)TIFFScanlineSize(tif),sizeof(unsigned char));
    g_row_mem = TIFFScanlineSize(tif);
  }
#endif
  if ( buf == NULL ) {
    fprintf(stderr,
           "\nCan't allocate memory for scanline buffer\n\n");
    mtfexit(-1);
  }
 
  /* enough to hold the whole image */
  g_image_array = (unsigned char*) calloc((int)(g_scan_rows*g_row_mem),sizeof(char));
  if ( g_image_array == NULL ) {
    fprintf(stderr, "\nCan't allocate memory for image array\n\n");
    mtfexit(-1);
  }
 
  if(g_debug) {
    /* enough to hold the whole image */
    g_image_array_modified =
      (unsigned char*) calloc((int)(g_scan_rows*g_scan_cols),sizeof(char));
    if ( g_image_array_modified == NULL ) {
      fprintf(stderr,"\n%s%s\n","Can't allocate memory for",
             " debug image array - try running with debug off");
      mtfexit(-1);
    }
  }
 
  /* As we read-in the test pattern area, we flip/rotate it to a standard
     format; the white patch at the right (top). There are two main cases
     - vertical (TOP or BOTTOM) and horizontal (RIGHT or LEFT), each
     indicating the white patch position in the original scanned image.
     The variables below are set to cause the flip/rotation to occur. These
     variables are used to index into the storage array as we process each
     scan line of the input image. */
  switch(rotation) {
  case TOP:
    /* white patch on top - make top rows into right cols */
    /* In the TOP and BOTTOM cases the height of
       the vertical input image is the width
       of our "normalized", horizontal form */
    /* Start on right hand side of storage array */
    /* g_input_height_mo etc. assigned in the initialize routine */
    array_column_start = g_input_height_mo;
    /* Work from right to left */
    array_column_delta = -1;
    /* We fill in the rightmost column first and work left. */
    /* Each row of input becomes a column in the storage array. */
    /* This is the distance from [0][0] to [1][0]  i.e. from */
    /* any row,col to the row,col directly below it. */
    /* Remember the height of the input is the width of the */
    /* storage array. */
    array_row_delta = 1;
    /* start at top (upper right) */
    array_row_start=0;
    break;
  case BOTTOM:
    /* white patch on bottom - make top rows left cols */
    /* The first pixel of the input image becomes the bottom
       pixel of the first column in the output array. */
    /* For each input row we fill in a storage column from
       bottom to top; working our way from left to right. */
    array_column_start = 0; /* left to right */
    array_column_delta = 1; /* ditto */
    /* bottom to top */
    array_row_delta = -1;
    /* bottom most row */
    array_row_start=g_input_width_mo;
    break;
  case RIGHT:
    /* normal case - white patch at top - read in straight*/
    array_column_start = 0;
    array_column_delta = 1;
    array_row_delta=1;
    array_row_start=0;
    break;
  case LEFT:
    /* white patch on bottom- flip top and btm; rt and lft*/
    /* right hand side of storage array */
    array_column_start = g_input_width_mo;
    /* right to left */
    array_column_delta = -1;
    /* bottom to top */
    array_row_delta=-1;
    /* bottommost row */
    array_row_start = g_input_height_mo;
    break;
  default:
    fprintf(stderr,"\nUnknown rotation!!\n\n");
    mtfexit(-1);
  }
 
  switch(rotation) {
    unsigned short *sbuf;
    unsigned short *inbuf;
 
  case TOP:
  case BOTTOM:
    array_column_index = array_column_start;
    switch(g_bytes_per_pixel) {
    case 1:
      for ( input_row = g_input_row_start;
           input_row <= g_input_last_row; input_row++) {
       read_scan_line(tif,buf,input_row);
       array_row_index=array_row_start;
       for (i=g_input_column_start;
            i<=g_input_last_column;i++) {
         g_image_array[array_column_index +
                     array_row_index *
                     g_input_height] = buf[i];
         array_row_index+=array_row_delta;
       }
       array_column_index+= array_column_delta;
      }
      break;
    case 2:
      sbuf = (unsigned short *) g_image_array;
      inbuf = (unsigned short *) buf;
      for ( input_row = g_input_row_start;
           input_row <= g_input_last_row; input_row++) {
       read_scan_line(tif,buf,input_row);
       array_row_index=array_row_start;
       for (i=g_input_column_start;
            i<=g_input_last_column;i++) {
         sbuf[array_column_index +
                     array_row_index *
                     g_input_height] = inbuf[i];
         array_row_index+=array_row_delta;
       }
       array_column_index+= array_column_delta;
      }
      break;
    default:
      fprintf(stderr,"\nCannot handle more than 2 bytes per pixel!!\n\n");
      mtfexit(-1);
    }
    break;
  case RIGHT:
  case LEFT:
    array_row_index = array_row_start;
    switch(g_bytes_per_pixel) {
    case 1:
      for ( input_row = g_input_row_start;
           input_row <= g_input_last_row; input_row++) {
       read_scan_line(tif,buf,input_row);
       array_column_index=array_column_start;
       for (i=g_input_column_start;
            i<=g_input_last_column;i++) {
         g_image_array[array_row_index *
                     g_input_width +
                    array_column_index] =
           buf[i];
         array_column_index+=array_column_delta;
       }
       array_row_index+=array_row_delta;
      }
      break;
    case 2:
      sbuf = (unsigned short *) g_image_array;
      inbuf = (unsigned short *) buf;
      for ( input_row = g_input_row_start;
           input_row <= g_input_last_row; input_row++) {
       read_scan_line(tif,buf,input_row);
       array_column_index=array_column_start;
       for (i=g_input_column_start;
            i<=g_input_last_column;i++) {
         sbuf[array_row_index *
                     g_input_width +
                     array_column_index] =
           inbuf[i];
         array_column_index+=array_column_delta;
       }
       array_row_index+=array_row_delta;
      }
      break;
    default:
      fprintf(stderr,"\nCannot handle more than 2 bytes per pixel!!\n\n");
      mtfexit(-1);
    }
    break;
  default:
    fprintf(stderr,"\nUnknown rotation!!\n\n");
    mtfexit(-1);
  }
 
  /* switch the corners of the test pattern to match the now rotated image */
  switch(rotation) {
  case TOP:
    hmo = g_scan_rows-1;
    tmpx = g_test_pattern_xul;
    tmpy = g_test_pattern_yul;
    g_test_pattern_xul = hmo - g_test_pattern_yll;
    g_test_pattern_yul = g_test_pattern_xll;
    g_test_pattern_xll = hmo - g_test_pattern_ylr;
    g_test_pattern_yll = g_test_pattern_xlr;
    g_test_pattern_xlr = hmo - g_test_pattern_yur;
    g_test_pattern_ylr = g_test_pattern_xur;
    g_test_pattern_xur = hmo - tmpy;
    g_test_pattern_yur = tmpx;
    /* we turned vertical format into horizontal - switch w + h */
    tmp=g_test_pattern_image_width;
    g_test_pattern_image_width=g_test_pattern_image_height;
    g_test_pattern_image_height=tmp;
    g_total_image_width = g_scan_rows; /* w&h are reversed*/
    g_total_image_height = g_scan_cols;
    break;
  case BOTTOM:
    hmo = g_scan_rows-1;
    wmo = g_scan_cols-1;
    tmpx = g_test_pattern_xul;
    tmpy = g_test_pattern_yul;
    g_test_pattern_xul = g_test_pattern_yur;
    g_test_pattern_yul = wmo - g_test_pattern_xur;
    g_test_pattern_xur = g_test_pattern_ylr;
    g_test_pattern_yur = wmo - g_test_pattern_xlr;
    g_test_pattern_xlr = g_test_pattern_yll;
    g_test_pattern_ylr = wmo - g_test_pattern_xll;
    g_test_pattern_xll = tmpy;
    g_test_pattern_yll = wmo - tmpx;
    /* we turned vertical format into horizontal - switch w + h */
    tmp=g_test_pattern_image_width;
    g_test_pattern_image_width=g_test_pattern_image_height;
    g_test_pattern_image_height=tmp;
    g_total_image_width = g_scan_rows; /*w&h are reversed*/
    g_total_image_height = g_scan_cols;
    break;
  case LEFT:
    hmo = g_scan_rows-1;
    wmo = g_scan_cols-1;
    tmpx = g_test_pattern_xul;
    tmpy = g_test_pattern_yul;
    g_test_pattern_xul = wmo - g_test_pattern_xlr;
    g_test_pattern_yul = hmo - g_test_pattern_ylr;
    g_test_pattern_xlr = wmo - tmpx;
    g_test_pattern_ylr = hmo - tmpy;
    tmpx = g_test_pattern_xur;
    tmpy = g_test_pattern_yur;
    g_test_pattern_xur = wmo - g_test_pattern_xll;
    g_test_pattern_yur = hmo - g_test_pattern_yll;
    g_test_pattern_xll = wmo - tmpx;
    g_test_pattern_yll = hmo - tmpy;
    g_total_image_width = g_scan_cols;
    g_total_image_height = g_scan_rows;
    break;
  case RIGHT:
    g_total_image_width = g_scan_cols;
    g_total_image_height = g_scan_rows;
    break;
  }
 
  /* copy the normalized image into the debug array */
  if (g_debug) {
    switch (g_bytes_per_pixel) {
    case 1:
      for (j=0;j<(int)(g_scan_rows*g_scan_cols);j++)
	if(g_image_array[j]==COLORA){
		g_image_array_modified[j]=COLORA-NUM_HIGH;
	}
	else if(g_image_array[j]==COLORB){
		g_image_array_modified[j]=COLORB+NUM_LOW;
	}
	else if(g_image_array[j]==COLORC){
		g_image_array_modified[j]=COLORC-NUM_HIGH;
	}
	else{
       		g_image_array_modified[j] = g_image_array[j];
	}
      break;
    case 2:
      {
       unsigned short *g_ptr;
       g_ptr = (unsigned short *)g_image_array;
       for (j=0;j<(int)(g_scan_rows*g_scan_cols);j++){  
	g_image_array_modified[j] = (g_ptr[j]>>(g_bps-8));
		if(g_image_array_modified[j]==COLORA){
			g_image_array_modified[j]=COLORA-NUM_HIGH;
		}
		else if(g_image_array_modified[j]==COLORB){
			g_image_array_modified[j]=COLORB+NUM_LOW;
		}
		else if(g_image_array_modified[j]==COLORC){
			g_image_array_modified[j]=COLORC-NUM_HIGH;
		}
	}
      }
      break;
    }
  }
 
  free(buf);
 
}
 
 
/******************************************************************************/
/*                                                                  */
/******************************************************************************/
 
/* Read a scan line into a buffer */
 
void read_scan_line(TIFF *tif,unsigned char *buf,short line_num)
{
  int j;
 
 
  if (g_scan_image_file_id) {
    if ( raw_readscanline(buf,line_num) < 0 ) {
      fprintf(stderr,
             "\nBad data read on imageline %d\n\n",line_num);
      mtfexit(-1);
    }
  }
#ifdef USE_TIFF
  else {
    if ( TIFFReadScanline( tif, buf, line_num, 0L ) < 0 ) {
      fprintf(stderr,
             "\nBad data read on imageline %d\n\n",line_num);
      mtfexit(-1);
    }
  }
#endif
 
  /* IF the tif tag: PHOTOMETRICINTERPRETATION was read as 'white is zero';
     Then the polarity needs to be reversed.  Also, the user may have
     requested a forced polarity reversal.  This can be useful when reading
     RAW images on little-endian or big-endian systems, or when a TIFF
     image has an incorrectly marked photometric interpretation.
     Two polarity reversals (tiff photometric=0 + user request) results in
     no reversal required.
  */
 
  if ( (g_reversepolarity || (g_photometric==0)) &&
      !(g_reversepolarity && (g_photometric==0)) ) {
      unsigned short *sbuf;
 
      switch (g_bytes_per_pixel) {
      case 1:
       for (j=0;j<(int)g_scan_cols;j++)
         buf[j] = g_max_pixel_value-buf[j];
       break;
      case 2:
       sbuf = (unsigned short *)buf;
       for (j=0;j<(int)g_scan_cols;j++)
         sbuf[j] = g_max_pixel_value-sbuf[j];
       break;
      default:
       fprintf(stderr,"\nCannot handle more than 2 bytes per pixel!!\n\n");
       mtfexit(-1);
      }
  }
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
static void write_debug_bounds(short rotation){ 
	double x[4],y[4];
	int c[4],r[4];
	double row,col;
	double sinsa,cossa;	
	int i,j;
	int corner1_x, corner1_y, corner2_x, corner2_y, corner3_x, corner3_y;
	short wmo, hmo;

	sinsa=sin(g_skew_avg);
	cossa=cos(g_skew_avg);
	for(i=0;i<g_number_of_patches;i++){
		x[0]=x[2]=g_test_pattern[i].xul_in_mm;
		x[1]=x[3]=g_test_pattern[i].xul_in_mm+g_test_pattern[i].width_in_mm;

		y[0]=y[1]=g_test_pattern[i].yul_in_mm;
		y[2]=y[3]=g_test_pattern[i].yul_in_mm+g_test_pattern[i].height_in_mm;

		for(j=0;j<4;j++){
			row=(y[j]/MM_PER_INCH)*g_ppi_vert;
			col=(x[j]/MM_PER_INCH)*g_ppi_horz;
			r[j]=NINT(col*sinsa+row*cossa+g_test_pattern_yul);
			c[j]=NINT(col*cossa-row*sinsa+g_test_pattern_xul);
			g_image_array_modified[r[j]*g_total_image_width+c[j]]=COLORA; /*draw points to make sure calculations are done correctly*/
		}
					draw_bound(c[0],r[0],c[1],r[1],COLORA,8);
					draw_bound(c[1],r[1],c[3],r[3],COLORA,8);
					draw_bound(c[3],r[3],c[2],r[2],COLORA,8);
					draw_bound(c[2],r[2],c[0],r[0],COLORA,8);
	}
 
    /* switch the 3 input corners to match the now rotated image,
       use temporary variables so that initial input data isn't lost. */
    switch(rotation) {
    case TOP:
      hmo = g_scan_rows-1;
      corner1_x = hmo - g_input_test_pattern_yul;
      corner1_y = g_input_test_pattern_xul;
      corner2_x = hmo - g_input_test_pattern_yur;
     corner2_y = g_input_test_pattern_xur;
      corner3_x = hmo - g_input_test_pattern_yll;
      corner3_y = g_input_test_pattern_xll;
      break;
    case BOTTOM:
      hmo = g_scan_rows-1;
      wmo = g_scan_cols-1;
      corner1_x = g_input_test_pattern_yul;
      corner1_y = wmo - g_input_test_pattern_xul;
      corner2_x = g_input_test_pattern_yur;
      corner2_y = wmo - g_input_test_pattern_xur;
      corner3_x = g_input_test_pattern_yll;
      corner3_y = wmo - g_input_test_pattern_xll;
      break;
    case LEFT:
      hmo = g_scan_rows-1;
      wmo = g_scan_cols-1;
      corner1_x = wmo - g_input_test_pattern_xul;
      corner1_y = hmo - g_input_test_pattern_yul;
      corner2_x = wmo - g_input_test_pattern_xur;
      corner2_y = hmo - g_input_test_pattern_yur;
      corner3_x = wmo - g_input_test_pattern_xll;
      corner3_y = hmo - g_input_test_pattern_yll;
      break;
    default:
      corner1_x = g_input_test_pattern_xul;
      corner1_y = g_input_test_pattern_yul;
      corner2_x = g_input_test_pattern_xur;
      corner2_y = g_input_test_pattern_yur;
      corner3_x = g_input_test_pattern_xll;
      corner3_y = g_input_test_pattern_yll;
      break;
    }
    /* Highlight the input corner pixels being used for processing */
    /* This is either the manual input, or the auto_refined value.
       Whichever is used for patch location calculations. */
    g_image_array_modified[corner1_x + g_total_image_width*corner1_y]=COLORB;
    g_image_array_modified[corner2_x + g_total_image_width*corner2_y]=COLORB;
    g_image_array_modified[corner3_x + g_total_image_width*corner3_y]=COLORB;


    /* Highlight the computed external corner pixels */
    g_image_array_modified[NINT(g_test_pattern_xul) + 
				g_total_image_width*NINT(g_test_pattern_yul)]=COLORC;
    g_image_array_modified[NINT(g_test_pattern_xur) + 
				g_total_image_width*NINT(g_test_pattern_yur)]=COLORC;
    g_image_array_modified[NINT(g_test_pattern_xll) + 
				g_total_image_width*NINT(g_test_pattern_yll)]=COLORC;
    g_image_array_modified[NINT(g_test_pattern_xlr) + 
				g_total_image_width*NINT(g_test_pattern_ylr)]=COLORC;
 
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/

void draw_bound(int col1, int row1, int col2, int row2, int color, int skip){
	int from, to, x, y, b;	/*b is initial (lowest) y or x*/
	double slope;

	if(col2-col1==0){	/*just in case*/
		if(row1>=row2){
			from=row2;
			to=row1;
		}
		else{
			from=row1;
			to=row2;
		}
		for(y=from;y<=to;y++){	
			if(skip) if(y%skip==0) continue;
			x=col1;	
			g_image_array_modified[y*g_total_image_width+x]=color;
		}	
	}
	else{
		slope=((double)(row2-row1))/((double)(col2-col1));
		if(fabs(slope)>=1){	/*whether to increment in units of y or x for smooth line*/
			if(row2>=row1){
				from=row1;
				to=row2;
			}
			else{
				from=row2;
				to=row1;
			}
			for(y=from;y<=to;y++){
				if(skip) if(y%skip==0) continue;
				if(col2>=col1) b=col1;
				else b=col2;	
				x=((double)(y-from)/slope)+(double)b;
				g_image_array_modified[y*g_total_image_width+x]=color;
			}
		}
		else{
			if(col2>=col1){
				from=col1;
				to=col2;
			}
			else{
				from=col2;
				to=col1;	
			}
			for(x=from;x<=to;x++){
				if(skip) if(x%skip==0) continue;
				if(row2>=row1) b=row1;
				else b=row2;
				y=(double)(slope*(x-from))+b;
				g_image_array_modified[y*g_total_image_width+x]=color;
			}
		}	
	}
}
