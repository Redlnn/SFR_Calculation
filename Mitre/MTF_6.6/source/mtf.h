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
/* Revisions to mtf.h                                                        */
/*                                                                           */
/* V1.1 General Clean-up  brp 11-MAR-94                                      */
/*                                                                           */
/* V1.1a Corrected comment for FIND_FIRST_PEAK  brp 15-MAR-94                */
/*       VERSION is now a string brp 17-MAR-94                               */
/*                                                                           */
/* V1.2 Refined the computation for the size of the "wiggle window" in       */
/*      FIND_PEAKS (it used to be 10% for low freq and 1.5 pixels for high   */
/*      freq.) Also provide a small error box around the predicted valley to */
/*      allow for small errors in the half-period calculation.  brp 21-MAR-94*/
/*                                                                           */
/* v1.3 Changed M_PI, etc. to PI, etc. and made local defines of them for    */
/*      portable code.  gtk 7-JUN-94                                         */
/*      Corrected "debug mode" bug in FIND_PEAKS - see comment in code       */
/*                                                              brp 7-JUN-94 */
/*                                                                           */
/* v2.0 FIND_FIRST_PEAK and FIND_PEAKS now accommodate images split int      */
/*      quadrants.  FIND_PEAKS now looks "ahead" and "back" in a single pass.*/
/*      FIND_PEAKS now examines peaks/valleys that are partially over the    */
/*      "edge" of the 10% safety border - it looks at them up to the border. */
/*                                                              brp 16-AUG-94*/
/*                                                                           */
/* v2.1 Added subroutine headers for round(), magnitude(), dabs(), sign()    */
/*                                                              gtk 23-AUG-94*/
/*                                                                           */
/* V2.3 Added declaration for get_switches().                                */
/*      Changed CORR_LOWER_LIMIT from 0.93 to 0.97.                          */
/*      Revised problem report printout.                                     */
/*      Added REPORT_SPEC_FOOTNOTE                                           */
/*                                                      DJB 6/7/95           */
/*                                                                           */
/* V3.0 Changed CORR_LOWER_LIMIT from 0.97 to 0.98.                          */
/*      Changed GREY_TO_REFLECTANCE to use Lagrangian interpolation, when    */
/*       appropriate.                                                        */
/*      Added declaration of lagrange_grey_to_reflectance(double).           */
/*      Changed REPORT_ALIASING to report upscaling, decimation, or both,    */
/*       where both is referred to as aliasing.                              */
/*      Indentation changes.                                                 */
/*                                                      DJB 7/27/95          */
/*      Don't report aliased frequencies.                                    */
/*                                                      DJB 12/28/95         */
/*      Added WITHIN_IMAGE_TOLERANCE.                                        */
/*      Fix problem with extra carriage returns in printout.                 */
/*                                                      DJB 1/22/96          */
/*      Delete image_array_tmp arg from do_white_patch().                    */
/*                                                      DJB 2/22/96          */
/*                                                                           */
/* V3.1 Changed EFF_FREQ_CUTOFF from 7.0 to 3.0.                             */
/*                                                      DJB 4/23/96          */
/*                                                                           */
/*      Removed copyright and added data rights notice.                      */
/*                                                      DJB 6/25/96          */
/*      Changed GREY_TO_REFLECTANCE to use local linear interpolation instea */
/*      of Lagrangian interpolation.                                         */
/*      Commented out Lagrangian interpolation stuff.                        */
/*                                                      DJB 8/2/96           */
/*      Added MAX_PATCH_LABEL_CHARS.                                         */
/*                                                      DJB 8/5/96           */
/* V4.0 Added support for 1000 ppi.					     */
/*	Added NYQUIST_DIVISOR for new TEST_FOR_ALIAS in sine.c               */
/*	get_orientation() now sets which patch is the white patch            */
/*							WE 9/99		     */
/*                                                                           */
/* V4.1 Added new parameter: PRINTER_TARGET_BASE, used in scaling            */
/*      calculations when printer option P is selected.                      */
/*      PRINTER_TARGET_BASE = 500 is the tgt ppi assumed in construction of  */
/*      the MITRE A-series digital printer targets.  Introducing this        */
/*      parameter allows a single digital target and data file to be used    */
/*      to obtain MTF of any resolution printer.                             */
/*                                                      nbn, 9/00            */
/*                                                                           */
/* V5.5 Made alias printing a part of the permanent problem report           */
/*      which is not removed when option 'n' is selected.                    */
/*                                                      MAL, 01/04           */
/* V5.5.1 Made alias report more explicit                                    */
/*        Use verbose mode 'e' instead of DEBUG_PEAK and DEBUG_PEAK_DETAILS. */
/*                                                      MAL, 05/04           */
/* V5.6   Added '80 columns' warning to Potential FATAL PROBLEMS             */
/*                                                      MAL, 06/04           */
/* V5.7   Made both \r and \n work within SKIP_BLANK.                        */
/*                                                      MAL, 07/04           */
/* V5.7.5 Added coefficients for MTF spec bounds formulas.                   */
/*                                                      MAL, 09/04           */
/* V5.8   Cleaned up code.  Removed unused sections.                         */
/*                                                      MAL, 09/04           */
/* V6.1  Added alias summary						     */
/*                                                      MAL, 09/05           */
/* V6.2  Added non-linear reflectance data for bartarget mode   	     */
/*                                                      MAL, 08/06           */
/* V6.3  Added PIV spec check and made no spec check the default   	     */
/*                                                      MAL, 10/06           */
/* V6.4  Fixed 6.3 bug in AppF CTF 1000 eqn coefficients                     */
/*                                                      MAL, 11/07           */
/* V6.4.01 Compilation without CodeWarrior, and using external TIFF library  */
/*                                                      MAL, 09/12           */
/* V6.4.02  Enabled reading PGM with H/W on separate lines.		     */
/*                                                      MAL, 11/13           */
/* V6.4.03  Added array bound checking for modulation arrays.		     */
/*                                                      MAL, 04/14           */
/* V6.5    Added colors and additional outlines to debug images.	     */
/*         Debug image outputs in BMP format if TIFF is not available.       */
/*                                                      AWU, 07/14           */
/* V6.6    Allow corner patches to be smaller if marked a 'c' instead of 'd'.*/
/*                                                      MAL, 04/16           */
/*****************************************************************************/
#define VERSION "6.6"
#define COLORA 255
#define COLORB 0
#define COLORC 254
#define NUM_HIGH 2
#define NUM_LOW 1
/*num high and low refer to the number of high colors (in 255 range) and number of low colors (in the 0 range).
 */
#define NUM_COLORS 256

#define COLORA_RED 102
#define COLORA_GREEN 194 
#define COLORA_BLUE 165

#define COLORB_RED 252
#define COLORB_GREEN 141
#define COLORB_BLUE 98

#define COLORC_RED 141
#define COLORC_GREEN 160
#define COLORC_BLUE 203

#define PROG_NAME "MTF" /* for DOS */

extern char *strerror(int);
/*
#ifndef linux
extern char *strncat(char *, const char *, size_t);
#endif
*/

/* These are assigned, not defined, for sake of processing speed */

#define ASSIGN_CORNERS                                                  \
        R1=g_test_pattern[i].r[0];                                      \
        R2=g_test_pattern[i].r[1];                                      \
        R3=g_test_pattern[i].r[2];                                      \
        R4=g_test_pattern[i].r[3];                                      \
                                                                        \
        C1=g_test_pattern[i].c[0];                                      \
        C2=g_test_pattern[i].c[1];                                      \
        C3=g_test_pattern[i].c[2];                                      \
        C4=g_test_pattern[i].c[3];

#define G_ASSIGN_CORNERS                                                \
        g_R1=g_test_pattern[i].r[0];                                    \
        g_R2=g_test_pattern[i].r[1];                                    \
        g_R3=g_test_pattern[i].r[2];                                    \
        g_R4=g_test_pattern[i].r[3];                                    \
                                                                        \
        g_C1=g_test_pattern[i].c[0];                                    \
        g_C2=g_test_pattern[i].c[1];                                    \
        g_C3=g_test_pattern[i].c[2];                                    \
        g_C4=g_test_pattern[i].c[3];

#define C_ABS(X,Y) sqrt(X*X+Y*Y)

#define DFT(in,n,z,out,datatype)                                        \
{                                                                       \
        register double w;                                              \
        register double uw;                                             \
        register double dc;                                             \
        register datatype *in_ptr;                                      \
        register int x;                                                 \
        register int u;                                                 \
        register int s1;                                                \
        register int s2;                                                \
        register int stop;                                              \
        register double real,img;                                       \
        register int k;                                                 \
        register double *out_ptr;                                        \
                                                                        \
        s1 = n;                                                         \
        s2 = n + z;     /* z = # of 0's to tack onto input data */      \
                                                                        \
        w=(double)PI*(double)2.0/(double)s2;                            \
                                                                        \
        in_ptr = (datatype *)in;                                        \
        out_ptr = (double *)out;                                         \
                                                                        \
        real=0.0;                                                       \
        for (x=0;x<s1;x++)                                              \
                real += (double)in_ptr[x];                              \
        dc = real/(double)s1;                                           \
        out_ptr[0]=dc;                                                  \
        out_ptr[1]= 0.0;                                                \
                                                                        \
        k=2;                                                            \
        stop=s2/2;                                                      \
        for (u=s2-1;u>=stop;u--) {                                      \
                real=img=0.0;                                           \
                uw=w*(double)u;                                         \
                for (x=0;x<s1;x++) {                                    \
                  real += ((double)in_ptr[x]-dc)*cos(uw*(double)x);     \
                  img  -= ((double)in_ptr[x]-dc)*sin(uw*(double)x);     \
                }                                                       \
                out_ptr[k++]=real;                               \
                out_ptr[k++]=img;                                \
        }                                                               \
}



void get_grey_reflectance_pairs(int number);
void set_grey_reflectance_pair(int i, double grey, double ref);
double GreyToReflectance(double grey_level);

void FindPeaks(double lc, double rc, double ppp, double pphp, 
                int *i_ptr, unsigned char *d_ptr, double fp,
                unsigned int i , int divisor);

#define GET_COORD_PAIR(a,b,c)                                         \
{                                                                       \
                int i;                                                  \
                scanf("%u",a);                                          \
                while(!isdigit(i=getchar())){}                           \
                ungetc(i,stdin);                                        \
                scanf("%u",b);						\
		while(!isdigit(i=getchar())){}                           \
                ungetc(i,stdin);                                        \
                scanf("%u",c); 						\
}

#define GET_FLOAT_PAIR(a,b)                                             \
{                                                                       \
                int i;                                                  \
                scanf("%f",a);                                          \
                while(!isdigit(i=getchar())){}                           \
                ungetc(i,stdin);                                        \
                scanf("%f",b);                                          \
}

#define SKIP_BLANK                                                      \
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;

#define GET_LINE                                                        \
        if (fgets(line,81,fd) == NULL) {                                \
                fprintf(stderr,"Short target data file!!\n");           \
                exit(-1);                                               \
        }

#define EUCLID_DIS(x1,y1,x2,y2)                                         \
        sqrt((double)((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)))


#define APPF 1
#define PIV 2

#define IQS 1			/* Type label for IQS related problems  */

#define REPORT_ALIASING                                                       \
  if (upscaling_detected || decimation_detected) {                            \
    if (upscaling_detected && decimation_detected)                            \
      sprintf(problem_string,"*Aliasing from both upscaling and decimation detected\n");                         \
    else if (upscaling_detected)                                              \
      sprintf(problem_string,"*Aliasing from upscaling detected\n");          \
    else if (decimation_detected)                                             \
      sprintf(problem_string,"*Aliasing from decimation detected\n");         \
                                                                              \
    put_problem(problem_string, 0);                                           \
    put_problem(" (See verbose output for details)\n\n", 0);                 \
  }

#define REPORT_SPEC_FOOTNOTE                                            \
        if (need_spec_footnote) {                                       \
                char problem_string[80];                                \
                /* need_spec_footnote also carries scanner resolution */\
                sprintf(problem_string,                                 \
                        "**Pertains to %s %d ppi scanner spec.\n",      \
                        g_spec, need_spec_footnote);                    \
                put_problem(problem_string, IQS);                       \
                put_problem("\n", IQS);                                 \
        }



#define MTFPRINT(a)                                                     \
{                                                                       \
        fprintf(stdout,a);                                              \
        fprintf(g_mtfout,a);                                            \
}

#define MTFPRINT2(a,b)                                                  \
{                                                                       \
        fprintf(stdout,a,b);                                            \
        fprintf(g_mtfout,a,b);                                          \
}

#define MTFPRINT3(a,b,c)                                                \
{                                                                       \
        fprintf(stdout,a,b,c);                                          \
        fprintf(g_mtfout,a,b,c);                                        \
}

#define MTFPRINT4(a,b,c,d)                                              \
{                                                                       \
        fprintf(stdout,a,b,c,d);                                        \
        fprintf(g_mtfout,a,b,c,d);                                      \
}

#define MTFPRINT5(a,b,c,d,e)                                            \
{                                                                       \
        fprintf(stdout,a,b,c,d,e);                                      \
        fprintf(g_mtfout,a,b,c,d,e);                                    \
}

#define MTFPRINT6(a,b,c,d,e,f)                                          \
{                                                                       \
        fprintf(stdout,a,b,c,d,e,f);                                    \
        fprintf(g_mtfout,a,b,c,d,e,f);                                  \
}

#define MTFPRINT7(a,b,c,d,e,f,g)                                        \
{                                                                       \
        fprintf(stdout,a,b,c,d,e,f,g);                                  \
        fprintf(g_mtfout,a,b,c,d,e,f,g);                                \
}

#define MTFPRINT8(a,b,c,d,e,f,g,h)                                      \
{                                                                       \
        fprintf(stdout,a,b,c,d,e,f,g,h);                                \
        fprintf(g_mtfout,a,b,c,d,e,f,g,h);                              \
}

#define MTFPRINT9(a,b,c,d,e,f,g,h,i)                                    \
{                                                                       \
        fprintf(stdout,a,b,c,d,e,f,g,h,i);                              \
        fprintf(g_mtfout,a,b,c,d,e,f,g,h,i);                            \
}

#define MTFPRINT10(a,b,c,d,e,f,g,h,i,j)                                 \
{                                                                       \
        fprintf(stdout,a,b,c,d,e,f,g,h,i,j);                            \
        fprintf(g_mtfout,a,b,c,d,e,f,g,h,i,j);                          \
}

#define MTFPRINT12(a,b,c,d,e,f,g,h,i,j,k,l)                             \
{                                                                       \
        fprintf(stdout,a,b,c,d,e,f,g,h,i,j,k,l);                        \
        fprintf(g_mtfout,a,b,c,d,e,f,g,h,i,j,k,l);                      \
}

#define TgtFileErrMsg 							\
{									\
  fprintf(stdout,"\n%s\n","***  FATAL PROBLEM IN TARGET DATA FILE  ***");  \
  fprintf(stdout,"%s\n","Not being read as pure text file (open file & re-save as text file), or"); \
  fprintf(stdout,"%s\n","Data is missing or out-of-sequence, or");	\
  fprintf(stdout,"%s\n","More than 80 columns in 1 line, or");	\
  fprintf(stdout,"%s\n\n","Comment lines have extra carriage returns or don't begin with # ");	\
} 

#define TgtFileBarMsg 							\
{									\
  fprintf(stdout,"\n%s\n","***  FATAL PROBLEM IN TARGET DATA FILE  ***");  \
  fprintf(stdout,"%s\n\n","Bar patches are combined with density or sine patches"); \
} 

#define TgtFileCnrMsg 							\
{									\
  fprintf(stdout,"\n%s\n","***  FATAL PROBLEM IN TARGET DATA FILE  ***");  \
  fprintf(stdout,"%s\n\n","If marking corner 'c' patches, there must be exactly 4"); \
} 

/* Moved into skew.c upper[]
#define MTFUL_CPMM 1.05 */ /* upper limit: same for all */
                        /* V5.6 individual upper limits for different 
             		   frequencies removed, since no longer used */

#define MTFLL_SMALL -99999.9999
#define MTFUL_BIG 99999.9999

#define MTF_A_1000 1.000121
#define MTF_B_1000 -0.07799797

#define LAND 0
#define PORT 1

#define TOP 0
#define BOTTOM 2
#define LEFT 1
#define RIGHT 3
#define UNKNOWN 255

#define SINE_PATCH 0
#define DENSITY_PATCH 1
#define DUPLICATE_DENSITY 2

#define NO_SINE_SIGMAS 2.0 /* used in one MTF calc.: (mean_mod+2*SD)/targt_mod*/

#define BEYOND_NYQUIST 999 /* flags a sine patch that has a frequency which
                              exceeds the Nyquist limit of the scanner */

#define EFF_FREQ_CUTOFF 3.0 /* above this skew angle IQS checks are not performed */

#define BORDER 0.1 /* % border around area of interest */

#define WITHIN_IMAGE_TOLERANCE 5 /* Outside corners of grey patches, calculated
                                    from user-supplied inner corners, must be
                                    within the image +/-WITHIN_IMAGE_TOLERANCE */

#define MTFOUTNAME "MTFOUT.txt"
#define TGTDATANAME "TGTDATA"

/* these two values are somewhat arbitrary - they define the dimensions of
   an array used to hold problem messages such as "aliasing occurs in
   sine patch X"  - one problem report per line */
#define MAX_PROBLEMS 60
#define MAX_LINE_LENGTH 132

/* If the linear regression error exceeds this value for a density patch,
   we issue a problem report */
#define DENSITY_MAX_LERR 7.65

#define PRINTER_TARGET_BASE 500 /* mm measurements in TgtDataFiles for MITRE's 
                                 A-series digital printer tgts correspond to 500ppi */
#define PPI_BASE 500 /* Current spec is for 500 ppi scanners */
#define PPI_THRESHOLD 750 /* threshold between 2 baselines */
#define PPI_HR_BASE 1000 /* high resolution spec is for 1000 ppi scanners */
#define PPI_ERROR 5 /* scanners outside PPI_BASE + or - PPI_ERROR are flaged */
#define PPI_HR_ERROR 10 /* high resolution scanners outside PPI_HR_BASE + or - P
                         PI_HR_ERROR are flaged */


/* if ratio of tgt reflectances for 2 density patches is < R_RATIO, then
the patches are averaged together in density.c (they have almost same density).
LOG10(R_RATIO) = threshold density difference between patches.
R_RATIO1 - avg patches that are within 0.040 density; for <17 patches, standard photo sine tgt
R_RATIO2 - avg patches that are within 0.008 density, for >=17 patches, this is set for digital sine tgt 
*/

 #define R_RATIO1 1.0964782
 #define R_RATIO2 1.0185914


#define THROW_OUT 0.1600 /* 16 % */

#define MIN_PPHP .950 /* below the Nyquist - but close enough to still test */

#define NYQUIST_DIVISOR 50.8


#define MAX_SECONDARY_TO_PRIMARY_RATIO 0.27 /* for freq below Nyquist, aliasing is present
                                               if ratio of secondary DFT peak magnitude
                                               to primary peak exceeds this limit. */

#define MAX_SECPRIM_ABOVE_NYQ 0.45   /* for freq marginally above Nyquist, aliasing is present
                                        if ratio of secondary DFT peak magnitude
                                        to primary peak exceeds this limit. 
					- nbn Apr'04 = 0.30 
				        - mal May'04 = 0.45 */


#define MAX_PATCH_LABEL_CHARS 11 /* Max number of characters in test_pattern.patch_label
                                    The string will have one more character, for the \0. */

#define ZERO_BOUND 0.001     /* Value below which a comparison is considered
			         equal to zero.  */
			      

/******************************************************************************/
/*                      prototypes                                            */
/******************************************************************************/


/* mtf.c */

void get_switches(int, char **,char*,char*);	
void get_args(char*,char*,TIFF**,unsigned int *,unsigned int *,unsigned int *);
void print_development_notice(void);
void print_data_rights_notice(void);
void print_help(void);
void print_header(char*,char*,unsigned char);
void print_problems(void);
void put_problem(char*, int);
void usage(char*);
void read_in_image(TIFF*, short,char*);
void read_scan_line(TIFF*,unsigned char*,short);
void mtfexit(int);
void draw_bound(int, int, int, int, int, int);

/* skew.c */

int NINT(double);
short initialize(char*,unsigned int,TIFF *);
void compute_patch_boundaries(void);
void compute_skew_angles(short);
void do_white_patch(unsigned char*,TIFF*,double*,unsigned int*,short,short,short,short);
short get_orientation(TIFF*,unsigned int);
void test_for_image_polarity(unsigned int);
void test_for_bar_target_polarity(unsigned int);
int get_target_data(char*,unsigned int*,unsigned int*,unsigned int*);
short raw_readscanline(unsigned char*,short);
void define_corner_patches(void);
void internal_to_external_corners(double, double, short);
void compute_correct_skews(short);
void approx_ppi(short);
void print_upright_info(char*);

/* density.c */

void compute_densities(unsigned int,unsigned int);
void compute_density_stats(int*,register unsigned int,unsigned int,unsigned int);
short next_ep_pair(register int);
void print_density_data(void);
void regress(void);
void print_histo(unsigned int);
void print_density_stats(void);
void sort_grey_reflectance_pairs(void);
void get_grey_reflectance_pairs(int);
void set_grey_reflectance_pair(int, double, double);
int bad_grey_check(double, int);
double local_linear_interpolate_grey_to_reflectance(double);

/* sine.c */

void compute_sines(void);
int alias_test(int*, double, int);
int compute_modulation(unsigned int);
void print_sine_data(void);
void print_sine_stats(void);

/* spline.c */

double compute_spline_grey_to_reflectance(double);
int smooth_spline(int *, double *);
void invert_spline(void);
void set_up_sort_reflectance_pairs(void);

/* scurve.c */

double compute_scurve_grey_to_reflectance(double);
double scurve_fit(void);
void least_squares(double *, double *, int, int, double *, double *);
int lu_decomp(double **, int, int *, double *);
void lu_solve(double **, double *, int, int *);

/* corner.c */
void refine_corner_coords(TIFF *);

#define MM_PER_INCH 25.4 /* 25.400 millimeters equal 1.000 inches */
#define PI 3.14159264348979323846
#define PI_2 1.57079632679489661923



