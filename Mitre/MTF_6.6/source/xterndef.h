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
/* Revisions to xterndef.h                                                    */
/*                                                                            */
/* V1.1 General Clean-up  brp 11-MAR-94                                       */
/*                                                                            */
/* V2.3 Added need_spec_footnote flag.                                        */
/*                                                      DJB 6/8/95            */
/*                                                                            */
/* V3.0 Added g_printer_ppi, g_printer, g_bad_correlation,                    */
/*      g_report_lagrangian_domain_problem.                                   */
/*                                                      DJB 8/16/95           */
/*      Added g_scan_header.                                                  */
/*                                                      DJB 12/28/95          */
/*      Added g_input_test_pattern_<xul, yul, xur, yur, xll, yll>.            */
/*      Delete g_test_pattern_<xul, yul, xur, yur>_old                        */
/*                                                      DJB 1/30/96           */
/*      Delete g_orient_debug.                                                */
/*                                                      DJB 2/22/96           */
/*                                                                            */
/* V3.1 Removed copyright and added data rights notice.                       */
/*                                                      DJB 6/25/96           */
/*      Changed length of test_pattern.patch_label from 12 to                 */
/*      MAX_PATCH_LABEL_CHARS + 1.                                            */
/*                                                      DJB 8/5/96            */
/* V4.0 Added support for 'n' and 's' switches.                               */
/*                                                      WE 9/99               */
/* V5.8   Cleaned up code.  Removed unused variables.                         */
/*                                                      MAL, 09/04            */
/*                                                                            */
/* V6.6     Added capability to pre-set corner density patches (via 'c')     */
/*          when they are significantly inset from minimal rectangle         */
/*           Changed corner parameters to floating point.                    */
/*                                                      MAL, 04/16           */
/******************************************************************************/

/******************************************************************************/
/*              arrays and file pointers                                      */
/******************************************************************************/
FILE *g_mtfout;
char g_debug_file_name[80];
char image_filename[80];
char data_filename[80];
int g_scan_image_file_id;
int g_ofd; 

unsigned char *g_image_array;
unsigned char *g_image_array_modified;

double *g_modulation_array;
unsigned short *g_mod_peak_array;
unsigned short *g_mod_valley_array;
int g_mod_array_cnt;

char g_problems[MAX_PROBLEMS][MAX_LINE_LENGTH];
short g_problem_count;
char g_problem_type[MAX_PROBLEMS];
short g_IQS_problem_count;

int g_test_pattern_corner[4];

/******************************************************************************/
/*      image, test pattern and patch locations and dimensions                */
/******************************************************************************/

short g_input_column_start;
short g_input_last_column;
short g_input_row_start;
short g_input_last_row;
short g_input_width_mo;
short g_input_width;
short g_input_height_mo;
short g_input_height;

unsigned int g_total_image_width;
unsigned int g_total_image_height;
unsigned int g_bytes_per_pixel;
unsigned int g_max_pixel_value;

float g_orig_test_pattern_xul;
float g_orig_test_pattern_yul;
float g_orig_test_pattern_xur;
float g_orig_test_pattern_yur;
float g_orig_test_pattern_xll;
float g_orig_test_pattern_yll;
float g_orig_test_pattern_xlr;
float g_orig_test_pattern_ylr;

float g_test_pattern_xul;
float g_test_pattern_yul;
float g_test_pattern_xur;
float g_test_pattern_yur;
float g_test_pattern_xll;
float g_test_pattern_yll;
float g_test_pattern_xlr;
float g_test_pattern_ylr;

float g_input_test_pattern_xul;
float g_input_test_pattern_yul;
float g_input_test_pattern_xur;
float g_input_test_pattern_yur;
float g_input_test_pattern_xll;
float g_input_test_pattern_yll;

double g_target_width_in_mm;
double g_target_height_in_mm;

double g_target_width;
double g_target_height;

int g_digital_camera;
double g_obj_img_scale;

float g_test_pattern_image_width;
float g_test_pattern_image_height;

int g_R1;
int g_R2;
int g_R3;
int g_R4;

int g_C1;
int g_C2;
int g_C3;
int g_C4;

int g_rc;
int g_lc;
int g_cr;


double g_l;


/******************************************************************************/
/*              skew and ppi                                                  */
/******************************************************************************/

double g_ppi_horz;
double g_ppi_vert;
int    g_target_res;

double g_patch_width;
double g_patch_height;

double g_skew_horz;
double g_skew_vert;
double g_skew_avg;

double g_high_throw_out;
double g_low_throw_out;

/******************************************************************************/
/*              density regression line                                       */
/******************************************************************************/

double g_slope;
double g_intercept;
double minreflect;
double maxreflect;


/******************************************************************************/
/*              things read in                                                */
/******************************************************************************/

unsigned int g_scan_rows;
unsigned int g_scan_cols;
unsigned int g_scan_header;
unsigned short g_bps;
unsigned short g_photometric;
int g_white_patch;
int g_row_mem;

unsigned char g_scan_orient;

unsigned int g_number_of_patches;
unsigned int number_unique_density_patches;
unsigned int number_of_averages;
unsigned int curve_type;
int **avg_index;

unsigned int g_printer_ppi;

unsigned char g_serial_number[82];
char g_spec[11]; 

unsigned char g_debug;
unsigned char g_print_average_mtf;
unsigned char g_extended;
unsigned char g_compare;
unsigned char g_compare_userset;
unsigned char g_printer;
unsigned char g_streamline_args;
unsigned char g_followcurve;
unsigned char g_user_requestcurve;
unsigned char g_reversepolarity;
unsigned char g_refine_corners;

unsigned short g_bartarget;
unsigned short g_nonlinear;
unsigned short g_corners_read;
int g_upright;

struct spline_reflectance_pair
{
  double grey;
  double reflectance;
  double slope;
  double quadratic;
  double cubic;
  int index;
} *sort_reflectance;

/* we allocate an array of these - one for each patch on the target - a lot
   of read-in info and computed info is stored here */

struct test_pattern {
        /* the first 9 entries are read from the target specific data file */
                        /* the first 6 apply to all patches */
        char patch_label[MAX_PATCH_LABEL_CHARS + 1];
        char new_patch_label[MAX_PATCH_LABEL_CHARS + 1];
        unsigned short patch_type; 

        double xul_in_mm;
        double yul_in_mm;
        double width_in_mm;
        double height_in_mm;
                        /* the next three are given by Sine Patterns Inc.*/
        double patch_density;   /* NA for sine pattern */
        double modulation;      /* NA for density patch */
        double cycles_per_mm;   /* NA for density patch */

	int corner_flag; /* whether patch is flagged as a corner (optional) */
        int rows; /* number of rows used in row averaging (sine.c) */

                /* we get values for next two from mtf.h */
        double mtfll;   /* spec lower limit for mtf */
        double mtful;   /* spec upper limit for mtf */

        /* most of the rest are computed from the scanned-in target */
        double grey_level;
        double final_reflectance;
        double final_grey;
        double Sspline;
        double Scurve;
        double computed_modulation_highest;
        double computed_modulation_average;
        double computed_modulation_average_plus;
        double effective_frequency;
        double mtf_highest;
        double mtf_average;
        double mtf_average_plus;
        unsigned char star;             /* '*' here if alias test fails */
        double reflectance;     /* 10**-density */
        double x[4];    /* "interior" corners of patch - 10% safety margin */
        double y[4];    /* these are in mm */
        int c[4];       /* row and col #s of above */
        int r[4];
        int low;
        int high;
        int total;
        int fliers;
        double first_sigma;
        double first_mean;
        double std_deviation;
        int sine_peak;
        int sine_valley;
        int *high_histo;
        int *low_histo;
		int totcases;
		int decim_over;
		double ratio_sum;
		double ratio_over_sum;
} *g_test_pattern;


/******************************************************************************/
/*              Flags                                                         */
/******************************************************************************/

int need_spec_footnote;

int g_UL_auto_crnr_err;
int g_UR_auto_crnr_err;
int g_LL_auto_crnr_err;
