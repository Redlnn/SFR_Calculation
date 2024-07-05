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
/* Revisions to sine.c                                                       */
/*                                                                           */
/* V1.1 - In find_peaks() and find_peaks_ns() sub_pix_range is now calculated*/
/*      to be 1.0 for 1.0 c/mm and 1.5 c/mm cases at 500 ppi.  This gives    */
/*      better tolerance to corner coordinate errors. brp 15-FEB-94          */
/*                                                                           */
/*      find_peaks() and find_peaks_ns() are now a single macro.             */
/*      added '*' to output for patches where alias test failed.             */
/*      reformatted to fit 80 cols; broke compute_sines() into two routines; */
/*      added macros, a few more comments. brp 9-MAR-94                      */
/*                                                                           */
/* V1.1a Added INTERNAL conditional compile for print_sine_data()            */
/*       internal version prints MTF using Avg Mod and Std Dev  brp 17-MAR-94*/
/*                                                                           */
/* V1.3 used NINT and locally defined PI in place of nint and M_PI           */
/*      made custom TEST_FOR_ALIAS based on ROW_INC for compute_sines_ns()   */
/*      made macros for rc calculations in compute_sines_by_row. Also        */
/*      adjusted left end calculation to only include r1 once for negative   */
/*      skew case.  Changed left column (lc) and right column (rc) calcs     */
/*      to use double casts of any integer values and to use floor and ceil  */
/*      functions to round the rc and lc values respectively.  Changed the   */
/*      TEST_FOR_ALIAS macro to use row1 instead of R1 as the first row. row1*/
/*      is assigned before the macro is called. brp 15-JUN-94                */
/*                                                                           */
/* V2.0 Two cases for no skew - 1.5 to 3 deg and <= 1.5 deg.                 */
/*      Also corrected alias test to accept the peak in either of the two    */
/*      highest freq. bins for samples slightly beyond nyquist.              */
/*                                                                           */
/*      Sine routines are now common regardless of skew!! A common endpoint  */
/*      finder, next_ep_pair is in density.c.  The old compute_sines_ns()    */
/*      and compute_sines() are now a single routine. brp 9-AUG-94           */
/*                                                                           */
/*      The global variables g_high_throw_out and g_low_throw_out, set in    */
/*      compute_density_stats (density.c), are used to exclude scan line     */
/*      values that have been edge enhanced or are noise spikes or are part  */
/*      of a "background" formed at the join of individually scanned quad-   */
/*      rants fused into a single image.                                     */
/*                                                      brp 16_AUG-94        */
/*                                                                           */
/* V2.3 Revised output printout.                                             */
/*      Revised problem report printout.                                     */
/*                                                      DJB 6/7/95           */
/*                                                                           */
/* V3.0 Added printout (when appropriate) regarding Printer MTF.             */
/*      Changed compute_sines and alias_test to report upscaling, decimation,*/
/*       or both, where both is referred to as aliasing.  No testing for     */
/*       decimation if frequency is above Nyquist.                           */
/*      Indentation changes.                                                 */
/*                                                      DJB 7/27/95          */
/*      Don't keep track of aliased frequencies.                             */
/*      If printer MTF (-p), don't compare to 500 ppi scanner spec.          */
/*      Conditionalize on USE_TIFF to not depend on tiff code.               */
/*                                                      DJB 12/28/95         */
/*      When frequency of a sine patch is beyond Nyquist frequency (last DFT */
/*       bin), the corresponding DFT location "folds over" back into the DFT */
/*       bins.  Look for the DFT peak at this folded location, +/- 1.        */
/*      Do not look for secondary peaks (due to decimation) in bin 1.  The   */
/*       value of this bin is often high due to bleeding of the DC component.*/
/*      Fix problem with extra blank lines in printout.                      */
/*                                                      DJB 1/24/96          */
/*                                                                           */
/* V3.1 In alias_test, if peak is not in last bin, check if secondary peak is*/
/*       there.                                                              */
/*      Added DEBUG_PEAKS code.                                              */
/*                                                      DJB 4/23/96          */
/*                                                                           */
/* V3.1 Removed copyright and added data rights notice.                      */
/*                                                      DJB 6/25/96          */
/* V4.0 TEST_FOR_ALIAS changed from sine patterns > 5 cy/mm to               */
/*              sine patterns > 30% of Nyquist frequency where               */
/*              Nyquist frequency = ppi / 50.8                               */
/*      Using new Sine pattern "row averaging" algorithm, check for aliasing */
/*              and finding peaks moved to subroutine elimitating the need   */
/*              for many of the macros                                       */
/*                                                      WE 9/99              */
/* V4.2 alias_and_peak: use the abs of diff between peak and valley incase we*/
/*      are off by 1/2 a cycle; move start pos over by 1/2 cycle if we did   */
/*      get off.                                                             */
/*      FindPeaks: NINT calls either floor or ceil, no need to call floor    */
/*      first, in fact, NINT is better than floor in these cases. Reset      */
/*      trail_low(1,2) only when needed, we were ignoring good data...       */
/*                                                      WE 3/01              */
/*                                                                           */
/* V5.0 Fixed bug in "row averaging" algorithm. Comparisons now in degrees.  */
/*                                                     mal 10/03             */
/*									     */
/* V5.5     Fixed a very obscure sleeper bug introduced in v4.52             */
/*          Comparison to determine if peak and valley need switching        */
/*          in alias_and_peak was made more robust.                          */
/*          Made alias problems report even when option 'n' selected.        */
/*                                                      mal, 01/04           */
/*                                                                           */
/* V5.5.1   Output dump of alias computations is now invoked with verbose    */
/*          mode 'e'. DEBUG_PEAK and DEBUG_PEAK_DETAILS no longer work.      */
/*                                                      mal, 05/04           */
/*                                                                           */
/* V5.6     Clarified extended alias output, moved to end of output.         */
/*          Added IQS check for CTF bartarget data                           */
/*          All extended output written to MTFOUT only, not stdout.          */
/*                                                      mal, 06/04           */
/*                                                                           */
/* V5.7     Clarified/shortened extended alias output.                       */
/*          Made alias check occur for all row averages.                     */
/*          Changed logic for setting decimation/upscaling_detected.         */
/*                                                      mal, 07/04           */
/*                                                                           */
/* V5.7.5   Made alias check use effective frequency for lower limit.        */
/*          Alias check is skipped if fewer than 2 periods present.          */
/*          Effective frequency computation moved to skew.c.                 */
/*          Spec check not performed if skew > 3 degrees, and spec           */
/*          warning message prints effective frequency now.                  */
/*          Effective frequency printout, uses actual eff. freq now.         */
/*                                                      mal, 09/04           */
/*                                                                           */
/* V5.8     Cleaned up code.  Sections no longer used were removed.          */
/*          Fixed unlikely but potential error in alias_count                */
/*	        introduced in v5.7. (Only if <99 rows signal upscaling)      */
/*                                                      mal, 09/04           */
/*                                                                           */
/* V6.0     Added 'g' option handling, for digital cameras                   */
/*                                                      mal, 06/05           */
/*                                                                           */
/* V6.1     Added aliasing summary to verbose mode                           */
/*                                                      mal, 09/05           */
/*                                                                           */
/* V6.2     Added non-linear reflectance data for bartarget mode   	     */
/*                                                      MAL, 08/06           */
/* V6.5     Added marking for last line of the sine computation.   	     */
/*                                                      AWU, 07/14           */
/*****************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <float.h>

#define TIFF char
#include "mtf.h"
#include "xterndcl.h"

FILE *mtftmp;
int limit_g_mod_array_cnt;
int totcases;
int decim_over;
double ratio_sum;
double ratio_over_sum;

/******************************************************************************/
/*                                                                            */
/******************************************************************************/

double GreyToReflectance (double grey_level)
{
  double retval;

  if (g_followcurve==9) {  /* S-curve fit */
    retval = compute_scurve_grey_to_reflectance(grey_level); 
  }
  else if (g_followcurve == 6) { /* Linear regression */
    retval = (grey_level/g_slope - g_intercept/g_slope);
    if (retval>1.0) 
      retval = 1.0;
    else if (retval<0.0) 
      retval = 0.0;
  }
  else if (g_followcurve == 5) /* Piecewise linear */
    retval = local_linear_interpolate_grey_to_reflectance(grey_level);
  else /* (g_followcurve == 7 || g_followcurve == 8)  Spline fit */
    retval = compute_spline_grey_to_reflectance(grey_level);

  return retval;
}

/* 
        Function: FindPeaks

        Work our way across a row using a measure that is one period long. We
        start with the first peak on the scan line (found by FIND_FIRST_PEAK),
        and move from peak to peak, computing the modulation of each peak and
        its corresponding valleys, 1/2 period away. These modulations are stored
        in array that is latter used in the patch's MTF computation after all
        rows in the patch have been processed.

        We use a real index since the period is real; rounding to an int only
        when we index into the array.  We allow some wiggle room at each end of
        the measure because ppp is not exact. 
*/
void FindPeaks(double lc, double rc, double ppp, double pphp, 
                                int *i_ptr, unsigned char *d_ptr, double fp,
                           unsigned int i , int divisor)
{                                                                            
  register int *a_ptr;                                                   
  register int trial_low1;                                                    
  register int trial_low2;                                                    
  register int trial_high;                                                    
  register int tl;                                                            
  register int diff;                                                          
  register int max1;                                                          
  register int max2;                                                          
  register int low1=0;                                                          
  register int high1=0;                                                         
  register int low2=0;                                                          
  register int high2=0;                                                         
  double fm,fk,fsum;                                                          
  double d_high,d_low;                                                        
  double max_ref,min_ref;                                                     
  double sub_pix=0.0,sub_pix_range,sub_pix_inc,h_period_del; 
  double dtmp, dtmp0, dtmp1;
  int itmp;
  void *new_array;
                                                                              
  a_ptr = (int *)i_ptr;                                                  
                                                                              
  /* sub_pix_range range is a - to + window centered on a "predicted          
     peak" and through which we determine the "true peak". The true           
     peak will usually be slightly different from the predicted position      
     because the computed ppp are not exact. */                               
                                                                              
  h_period_del =.01*pphp;       /* compensate for small hppi err */           
                                                                              
  /* @495 ppi you get .974409 pphp for 10c/mm; 1.948819 ppp.                  
     We use a base of .75 because it seems to work well and                   
     1.948819/20 + .75 < .974409 keeping us (in theory) within + or -         
     a half period in wiggling the peak in our error window. */               
                                                                              
  dtmp = ppp/20.0;
  sub_pix_range = dtmp+0.75;                                              
  if (sub_pix_range < 1.5) {                                                  
    dtmp = sub_pix_range/2.0;
    sub_pix_inc = 0.9999999999*dtmp;                           
  }                                                                           
  else {                                                                      
    sub_pix_range = (double)NINT(sub_pix_range);                              
    sub_pix_inc = 1.0;                                                        
  }                                                                           
  for (fk=fp;fk<=rc;fk+=ppp) {                                                
    max1=max2=0;                                                              
      /* big enough initial value so that "diff" will be negative for those   
         cases where we never find a trial low. Should yield 32767 on a 16    
         bit compiler, 98303 on a 32 bit (normal case) - both are big, but    
         98303 is large enough to be larger than the sum of many bytes when   
         we sum rows for the non-skewed case.  
	 NEW: just use system maximum value
      */
    trial_low1=trial_low2 = INT_MAX;                   
    for (sub_pix=-sub_pix_range;sub_pix<=sub_pix_range;                       
         sub_pix+=sub_pix_inc) {                                              
      dtmp = fk + sub_pix;
      dtmp0 = NINT(dtmp);
      if (dtmp0 < lc)                                           
        continue;                                                             
      if (dtmp0 > rc)                                           
        break;                                                                
      /* first consider each potential "highest pixel" in the "wiggle window" 
         around our predicted peak at fk */                                   
      itmp = (int)dtmp0;
      trial_high = (int)*(a_ptr + itmp);
      /*trial_high = (int)*(a_ptr + NINT(floor(fk + sub_pix)));*/
      /* we may have a fracture in the center of the image for images that    
         were scanned on small scanners and pasted together.  If we detected  
         any out of range pixels in compute_sines, we've set them to a neg-   
         ative value - we won't process them here. */                         
      if (trial_high < 0) {                                                   
        trial_high = 0;                                                       
        continue;                                                             
      }                                                                       
      if (g_debug) {
        *(d_ptr + itmp) = 0; 
        /**(d_ptr+NINT(floor(fk+sub_pix)))=0;*/
      }
      /* corresponding to each candidate "highest pixel" is a valley pixel    
         which is 1/2 period up from the peak */                              
      fsum = fk+sub_pix+pphp;                                                 
      /* look a little on either side of the predicted valley to compensate   
         for small horizontal ppi errors that would throw the hppi value      
         slightly off.  h_period_del is 1% of the pphp. In some cases this    
         may be enough to carry us into an adjacent pixel to the predicted    
         valley pixel (we still consider the predicted pixel). */             
      dtmp = fsum-h_period_del;
      dtmp1 = fsum+h_period_del;
      for (fm=dtmp;fm<=dtmp1;fm+=h_period_del) {
        /* fk has moved far enough across so that next valley up is past rc */ 
        dtmp0 = NINT(fm);
        if (dtmp0 > rc)                                                   
          break;                                                              
        itmp = (int)dtmp0;
        tl = (int)*(a_ptr + itmp);
        /*tl = (int)*(a_ptr + NINT(floor(fm)));*/
        /* we're at a split in the image - see comments about trial_high above*/
        if (tl < 0) 
          continue;  
        if (g_debug) {
          *(d_ptr + itmp) = 255;
          /**(d_ptr+NINT(floor(fm)))=g_max_pixel_value;*/
        }
        if (tl < trial_low1)                                                  
          trial_low1=tl;                                                      
      }                                                                       
      /* if we never find a trial high (trial high = 0) or never find a trial 
         low (trial low = 98303 - then diff will be negative and not greater  
         than max1 */                                                         
      diff = trial_high-trial_low1;                                           
      if (diff > max1) {                                                      
        max1=diff;                                                            
        high1 = trial_high;                                                   
        low1 = trial_low1;                                                    
      }                                                                       
      /* look back at previous valley from peak at fk */                      
      fsum = fk+sub_pix-pphp;                                                 
      dtmp = fsum-h_period_del;
      dtmp1 = fsum+h_period_del;
      for (fm=dtmp;fm<=dtmp1;fm+=h_period_del) {
        /* fk is not yet far enough to right to have a previous valley */     
        dtmp0 = NINT(fm);
        if (dtmp0 < lc)                                                   
          continue;                                                           
        itmp = (int)dtmp0;
        tl = (int)*(a_ptr + itmp);
        /*tl = (int)*(a_ptr + NINT(floor(fm)));*/
        if (tl < 0) 
          continue;
        if (g_debug) {
          *(d_ptr + itmp) = 255;
          /**(d_ptr+NINT(floor(fm)))=g_max_pixel_value;*/
        }
        if (tl < trial_low2)                                                  
          trial_low2=tl;                                                      
        /* we keep the two trial_lows/maxs separate so that each peak valley  
           pair (back and ahead) will be considered when computing the average 
           mtf */                                                             
      }                                                                       
      diff = trial_high-trial_low2;                                           
      if (diff > max2) {                                                      
        max2=diff;                                                            
        high2 = trial_high;                                                   
        low2 = trial_low2;                                                    
      }                                                                       
    }         
    /* This avoids overrunning the end of these arrays. */
    if (limit_g_mod_array_cnt <= g_mod_array_cnt + 1) { 
      fprintf(stderr,"\nLimit_g_mod_array_cnt = %d is too small.  Increasing by 16.\n\n", limit_g_mod_array_cnt );
      limit_g_mod_array_cnt += 16;
      if ((new_array = realloc (g_modulation_array, limit_g_mod_array_cnt*sizeof(double)))== NULL) {
		fprintf(stderr,"\nCan't re-allocate memory for g_modulation_array!!\n\n");
		mtfexit(-1);
      }
      else g_modulation_array = (double *)new_array;
      if ((new_array = realloc (g_mod_peak_array, limit_g_mod_array_cnt*sizeof(short)))== NULL){
		fprintf(stderr,"\nCan't re-allocate memory for g_mod_peak_array!!\n\n");
		mtfexit(-1);
      }
      else g_mod_peak_array = (unsigned short *)new_array;
      if ((new_array = realloc (g_mod_valley_array, limit_g_mod_array_cnt*sizeof(short)))== NULL){
		fprintf(stderr,"\nCan't re-allocate memory for g_mod_valley_array!!\n\n");
		mtfexit(-1);
      }
      else g_mod_valley_array = (unsigned short *)new_array;
    }
    if (max1 != 0) {                                                          
      d_high = (double)high1/(double)divisor;                                 
      d_low = (double)low1/(double)divisor;
      if (g_bartarget && g_nonlinear) {
	bad_grey_check(d_high, 1);
	bad_grey_check (d_low, 1);
      }
      max_ref = GreyToReflectance(d_high);
      min_ref = GreyToReflectance(d_low);                                   
      g_modulation_array[g_mod_array_cnt] =
        (max_ref-min_ref)/(max_ref+min_ref);
      g_mod_peak_array[g_mod_array_cnt] = NINT(d_high);
      g_mod_valley_array[g_mod_array_cnt++] = NINT(d_low);
      if(g_extended){
	g_test_pattern[i].high_histo[NINT(d_high)]++;                
	g_test_pattern[i].low_histo[NINT(d_low)]++;                 
      }
    }                                                                         
    if (max2 != 0) {                                                          
      d_high = (double)high2/(double)divisor;                                 
      d_low = (double)low2/(double)divisor;                                   
      if (g_bartarget && g_nonlinear) {
	bad_grey_check (d_high, 1);
	bad_grey_check (d_low, 1);
      }
      max_ref = GreyToReflectance(d_high);                                  
      min_ref = GreyToReflectance(d_low);                                   
      g_modulation_array[g_mod_array_cnt] =                                 
        (max_ref-min_ref)/(max_ref+min_ref);                                   
      g_mod_peak_array[g_mod_array_cnt] = NINT(d_high);
      g_mod_valley_array[g_mod_array_cnt++] = NINT(d_low);
      if(g_extended){
	g_test_pattern[i].high_histo[NINT(d_high)]++;                   
	g_test_pattern[i].low_histo[NINT(d_low)]++;                  
      }
    }                                                                         
  }                                                                           
}

/* 
   Perform aliasing test (if frequency is large enough) and then
   find peaks and valleys.   Return value indicates whether upscaling
   or decimation was detected. 1 for upscaling, 100 for decimation.
*/
short alias_and_peak(int i, double lower_alias_limit, int *row_array,
                     short prev_alias_count, 
                     double ppp, double pphp,
                     int start_c, int end_c, short right_col, int s_row,
                     int number_of_rows_averaged)
{
   double fp;                    /* first peak in a row */
   short alias_count = 0;
   unsigned char *image_debug_ptr;

   /* Test all cases above the lower alias limit with large enough pphp. v5.7 
      Make sure to use effective_frequency for check rather than
      target file frequency. v5.7.5 */
   if (g_test_pattern[i].effective_frequency >= lower_alias_limit &&
       (pphp > MIN_PPHP) ) 
     alias_count = alias_test(row_array, ppp, (int)(end_c-start_c));

   /* If there are already 99 instances of upscaling, then it is impossible
      to add any more upscaling to alias_count  v5.8 */
   if (prev_alias_count%100 == 99) 
     alias_count = (alias_count/100) * 100;

   /* If there are already 654 instances of decimation, then it is impossible
      to add any more decimation to alias_count v5.8 */
   if (prev_alias_count/100 == 654) 
     alias_count %= 100;
   

/* Find the 1st peak on a line from a sine patch.  The location of the 1st
   peak is a float, since the pixels per period (ppp) and half period (pphp)
   are floats as well.  We look for the 1st peak in a window that is 1 period
   long. This window is sub-divided into bins with a 1.0 pixel resolution. For
   each bin we sum bin + (bin + period) + (bin + 2 periods) ...  as well as a
   separate sum for (bin + half-period) + (bin + 3/2period) ... The bin
   having the largest average difference between two such sums is the 1st peak.
   Note that row_array can be char or int (skew or no-skew); datatype tells
   which. */
#define datatype int
   {
     double lc = 0.0, rc = right_col;
     double fm,fk;
     struct score {
             double peak;
             double valley;
             double count;
     } score;
     double winner,avg=0.0;
     double pos;
     datatype *row_ptr;

     row_ptr = row_array;
     pos=-1.0;
     winner = 0.0;
/* The increment can be made less than 1.00 in order to do sub pixel
   binning. This works, and indeed, the fp returned is often fractional
   instead of integer.  However, the refinement is absorbed in the
   overall peak finding process.  Even when the HPPI was intentionally
   made incorrect, there seemed to be little difference in the MTFs
   when fractional increments were used. */
     for(fk=lc;fk<lc+ppp;fk+=1.00) {
        score.peak=0.0;
        score.valley=0.0;
        score.count=0.0;
        for (fm=fk;fm<rc-pphp;fm+=ppp) {
           if ((double)row_ptr[NINT(fm)] >= 0 &&
               (double)row_ptr[NINT(fm+pphp)] >= 0) {
              score.peak   += (double)row_ptr[NINT(fm)];
              score.valley += (double)row_ptr[NINT(fm+pphp)];
              score.count++;
           }
        }
        if (score.count > 0.0) {
/* we use an average because some pixels at the right of the window
   may have one fewer peak-valley pairs to sum with than those on the
   left. */
           avg = (score.peak-score.valley)/score.count;
           avg = fabs(avg);
	   /* Only do this comparison when score data available MAL 01/15/04 */
	   if (avg > winner) {
	      winner = avg;
              pos = fk;
	      /* Bug introduced in v4.2.  If row_ptr[pos]==-1 then 
		 this causes the wrong behavior.  Instead check the 
		 score sums to see if peak/valley need reversing. v5.5 */
	      /* peak = (double)row_ptr[NINT(pos)]; */
	      /* valley = (double)row_ptr[NINT(pos+pphp)]; */
	      /* if (valley > peak) */
              if (score.valley > score.peak) 
		 pos += pphp; /* we are off by half a peak */
	   }
        }
      }

      fp = pos;
   }

   image_debug_ptr = NULL;
   if(g_debug)
     image_debug_ptr = 
       &g_image_array_modified[s_row*g_total_image_width+start_c];
   
   if (fp >= 0.0) {
       FindPeaks((double)0.0,(double)right_col,ppp,pphp,row_array,
                 image_debug_ptr,
		 fp,i,number_of_rows_averaged);
   }
   return (alias_count);
}

/* Windows compilers without rint (MinGW is OK) */
#if defined (_WIN32) & !defined(__GNUC__)  
double rint(x)
double x;
{
  int y = (int)x; 
  double half = (double)y + 0.5;
  double ans;

  if (x >= half)
    ans = (double)y + 1.0;
  else
    ans = (double)y;
  return (ans);
}
#else
extern double rint(double);
#endif

int determine_number_of_rows_averaged(double ppi, double nyquist,
                                      double cycles_per_mm, double abs_skew)
{
  int ans;
  double F, R, p;

  /* Convert to degrees v5.0 - mal */
  abs_skew *= 180.0/PI;

  /* F =  normalized frequency = cycles_per_mm / Nyquist
                               = cycles_per_mm * 50.8 / ppi */
  F = cycles_per_mm / nyquist;

  /* R = number of rows (rounded to whole number) */
  if ((0.0 <= abs_skew) && (abs_skew <= 1.0)) {
    p = pow(F,-1.015904);
    R = 5.9970579 * p;
  } else if ((1.0 < abs_skew) && (abs_skew <= 2.0)) {
    p = pow(F,-1.029771);
    R = 2.8639964 * p;
  } else if ((2.0 < abs_skew) && (abs_skew <= 3.0)) {
    p = pow(F,-1.03381);
    R = 1.8645913 * p;
  } else if ((3.0 < abs_skew) && (abs_skew <= 5.0)) {
    p = pow(F,-1.06002);
    R = 1.0397801 * p;
  } else { /* (abs_skew > 5.0) */
    R = 1.0;
  }

  if (R > (0.1 * ppi)) {
    R = 0.1 * ppi;
  }

  ans = (int)rint(R);
  if (ans < 1)
    ans = 1;
  return (ans);
}

void compute_sines(void)
{
  register int start_c=0, end_c=0;
  register int s_row=0;

  char problem_string[82];
  char freq_type[10];
  int modulation_problems;
  unsigned int i;
  short j;
  int fk;
  short right_col=0;
  short alias_count;
  int decimation_detected = 0; 
  int upscaling_detected = 0;
  int decimation_detected_this_pass, upscaling_detected_this_pass;

  static int row_array_data[4096];
  int *row_array;

  double ppp;                   /* pixels_per_period */
  double pphp;                  /* pixels_per_half_period */
  double ppi;
  double cycles_per_mm;
  double min_pixels_per_row;
  double nyquist;
  double abs_skew_avg = fabs(g_skew_avg);
  double cos_skew_avg = cos(g_skew_avg);
  double inside_area = 1.0 - 2.0*(double)BORDER;
  int number_of_rows_averaged;
  short save_row_count, row_count;
  int tmp,ii;
  int bogus_cnt=0;
  int num_non_singleton_averages=0;

  /* lower_alias_limit = 30% of nyquist frequency in cy/mm 
    nyquist_frequency in cy/mm = g_ppi_horz / NYQUIST_DIVISOR; */
  double lower_alias_limit = (g_ppi_horz / NYQUIST_DIVISOR) * 0.30;

  /* This revised row averaging performs the maximum amount of noise smoothing
     that is consistent with a maximum  1/2 percent reduction in MTF modulation
     values. */
  int testing_for_alias;

  /* Set frequency type string for correct printout labelling and adjust 
  	 lower_alias_limit if needed. */
  if ( g_digital_camera ) 
  {
  	lower_alias_limit *= g_obj_img_scale;
  	sprintf(freq_type,"    image");
  }
  else sprintf(freq_type,"effective");

  /* Create a temporary file to hold extended alias debug info temporarily. */
  if (g_extended) 
    mtftmp = fopen("_MTF_TMP", "w");

  if (g_scan_orient == LAND) 
    ppi = g_ppi_horz;
  else 
    ppi = g_ppi_vert;

  /* Nyquist in cy/mm = ppi / 50.8 */
  nyquist =  ppi / NYQUIST_DIVISOR;

  for (i=0;i<4096;i++) row_array_data[i]=-99;
  row_array=&row_array_data[64];

  modulation_problems = 0;
  decimation_detected = upscaling_detected = 0;

  for (i=0;i<g_number_of_patches;i++) {

    if (g_test_pattern[i].patch_type != SINE_PATCH) 
        continue;


    if(g_extended){
      g_test_pattern[i].high_histo = (int *)calloc(g_max_pixel_value+1,sizeof(int));
      g_test_pattern[i].low_histo = (int *)calloc(g_max_pixel_value+1,sizeof(int));
	  totcases = 0;
	  decim_over = 0;
	  ratio_sum = 0.0;
	  ratio_over_sum = 0.0;
    }
    
    cycles_per_mm = g_test_pattern[i].cycles_per_mm;

    number_of_rows_averaged = determine_number_of_rows_averaged(ppi, nyquist,
                                                                cycles_per_mm,
                                                                abs_skew_avg);
    g_test_pattern[i].rows = number_of_rows_averaged;

    if (number_of_rows_averaged > 1) num_non_singleton_averages++;
    
    if (g_extended)
      fprintf(mtftmp, "\nTgt cy/mm: %f     %s cy/mm = %f\n",
	      cycles_per_mm, freq_type, g_test_pattern[i].effective_frequency);

    decimation_detected_this_pass = 0;
    upscaling_detected_this_pass = 0;
    testing_for_alias = 0;
    alias_count = 0;
    g_test_pattern[i].star=' ';

    /* we store the corner coords into easy to use variables with short names.
       The variables are global to allow access by the enpoint finder
       (next_ep_pair) without passing a lot of parameters. */
    G_ASSIGN_CORNERS          /*g_R1 = g_test_pattern[i].r[0],etc. */

    /* g_l is the unskewed 80% patch width in pixels. */
    /* inside_area = 80% = 0.8 */
    g_l = (g_test_pattern[i].width_in_mm/MM_PER_INCH)*g_ppi_horz * inside_area;
      
    /* pixels per period */
    ppp = g_ppi_horz / (cos_skew_avg * cycles_per_mm * MM_PER_INCH);
          
    /* we got to have this many pixels in a row before we can look at it */
    min_pixels_per_row = ppp;
      
    /* pixels/half-period */
    pphp = ppp / 2.0;
    if (pphp < MIN_PPHP) {
      g_test_pattern[i].computed_modulation_highest = BEYOND_NYQUIST;
      g_test_pattern[i].mtf_highest = 0.0;
      continue;
    }
      
      
    /* safe (overkill) count of how much memory to allocate for modulation
       storage */
    if (g_skew_avg < 0)
      limit_g_mod_array_cnt = (g_R3-g_R2)*4*NINT((double)(g_C4-g_C1)/(ppp*cos_skew_avg));
    else
      limit_g_mod_array_cnt = (g_R4-g_R1)*4*NINT((double)(g_C2-g_C3)/(ppp*cos_skew_avg));

    /* modulation for entire patch stored here - each row can have many
       modulations, one for each peak/valley in row */

    g_modulation_array = (double*)calloc((int)(limit_g_mod_array_cnt),sizeof(double));

    if ( g_modulation_array == NULL ) {
      fprintf(stderr,"\nCan't allocate memory for g_modulation_array!!\n\n");
      mtfexit(-1);
    }
      
    g_mod_peak_array = (unsigned short *)
      calloc((int)(limit_g_mod_array_cnt),sizeof(unsigned short));

    if ( g_mod_peak_array == NULL ) {
      fprintf(stderr,"\nCan't allocate memory for g_mod_peak_array!!\n\n");
      mtfexit(-1);
    }
      
    g_mod_valley_array = (unsigned short *)
      calloc((int)(limit_g_mod_array_cnt),sizeof(unsigned short));

    if ( g_mod_valley_array == NULL ) {
      fprintf(stderr,"\nCan't allocate memory for g_mod_valley_array!!\n\n");
      mtfexit(-1);
    }
      

    g_mod_array_cnt=0;
      
    g_cr=-1;
    save_row_count = row_count = 0;
    /* work thru successive scan lines of a patch, getting new left and right
       endpoints for the scan line. next_ep_pair returns 0 when it runs out
       of scan lines  - see density.c to look at next_ep_pair code. */
    while(next_ep_pair(i)) {
	/*uncomment below to see a filled block of where rows are not averaged*/
	/*if(g_cr+g_test_pattern[i].rows>g_test_pattern[i].r[3] && g_debug) draw_bound(g_lc,g_cr,g_rc,g_cr,COLORC,0);*/

      if ((double)(g_rc - g_lc) < min_pixels_per_row) {
        save_row_count = row_count;
        row_count=0;
        continue;
      }
      j=0;
      /* row_count is reset to 0 once we have processed a the number of rows
         averaged (of scan lines). */
      if (row_count == 0) {

	/*rows not used for sine calculations*/
	if(g_cr+g_test_pattern[i].rows>g_test_pattern[i].r[3] && g_debug) draw_bound(g_lc,g_cr,g_rc,g_cr,COLORC,2);

        /* always use the same start and end points for successive rows (for
           the low angle cases where we sum rows.)  The 10% border gives us
           enough leeway to do this even though some successive endpoint may
           be a unit more or less than the 1st endpoint. */
        start_c = g_lc;
        end_c   = g_rc;
        s_row   = g_cr;
        right_col=end_c-start_c;
        /* we'll get a bogus MTF if we ever overrun our borders in FIND_PEAKS */
        /* row_array starts at row_array_data[64] */
        for (ii = 0; ii < 65; ii++) {
          row_array_data[ii]=-99;
          row_array_data[ii+right_col+63]=-99;
        }
	if (right_col+64+64>4096) {
	  MTFPRINT2("ERROR: Static row_array_data is too small. Need %d\n",
		    right_col+128);
	}
	
	switch(g_bytes_per_pixel) {
	  unsigned short *sbuf;
	  
	case 1:
	  for (fk = start_c; fk <= end_c; fk++) {
	    tmp = (int)g_image_array[g_cr*g_total_image_width+fk];
	    if (tmp > (int)g_high_throw_out || tmp < (int)g_low_throw_out) {
	      /* Logic is being changed so that: 
		 -99 = Out of bounds
		 -1 = Bogus value encountered during sum.
		 Nothing will ever be summed into any value set to -1. */
	      tmp = -1;
	      bogus_cnt++;
	    }
	    row_array[j++] = tmp;
	  }
	  break;
	  
	case 2:
	  sbuf = (unsigned short *)g_image_array;
	  for (fk = start_c; fk <= end_c; fk++) {
	    tmp = (int)sbuf[g_cr*g_total_image_width+fk];
	    if (tmp > (int)g_high_throw_out || tmp < (int)g_low_throw_out) {
	      /* Logic is being changed so that: 
		 -99 = Out of bounds
		 -1 = Bogus value encountered during sum.
		 Nothing will ever be summed into any value set to -1. */
	      tmp = -1;
	      bogus_cnt++;
	    }
	    row_array[j++] = tmp;
	  }
	  break;
	}
      }
      /* for the multiple lines case we want to sum in successive lines */
      else { /* row_count != 0 */
	switch(g_bytes_per_pixel){
	  unsigned short *sbuf;
	  
	case 1:
	  for (fk = start_c; fk <= end_c; fk++) {
	    tmp = (int)g_image_array[g_cr*g_total_image_width+fk];
	    if (row_array[j] < 0) /* Never touch negative values */
	      j++;
	    else if (tmp > (int)g_high_throw_out || tmp < (int)g_low_throw_out)
	      {
		bogus_cnt++;
		row_array[j++] = -1;
	      }
	    else 
	      row_array[j++] += tmp;
	  }
	  break;
	case 2:
	  sbuf = (unsigned short *)g_image_array;
	  for (fk = start_c; fk <= end_c; fk++) {
	    tmp = (int)sbuf[g_cr*g_total_image_width+fk];
	    if (row_array[j] < 0) /* Never touch negative values */
	      j++;
	    else if (tmp > (int)g_high_throw_out || tmp < (int)g_low_throw_out)
	      {
		bogus_cnt++;
		row_array[j++] = -1;
	      }
	    else 
	      row_array[j++] += tmp;
	  }
	  break;
	}
      }
      row_count++;
        
      /* if we have a complete set of summed rows then go get the modulations.
         The macros below are defined above.  They in turn, use macros defined
         in mtf.h. */
      if ((number_of_rows_averaged != 0) &&
          !(row_count % number_of_rows_averaged)) {


	/*the below will display where the calculations are done on the averaged rows*/
	/*if(g_debug) draw_bound(g_lc,g_cr,g_rc,g_cr,COLORC,4);*/

        testing_for_alias = 1;
        row_count = 0;
        alias_count += alias_and_peak(i,lower_alias_limit,row_array,
                                      alias_count, 
                                      ppp,pphp,start_c,end_c,right_col,s_row,
                                      number_of_rows_averaged);
      }
    }

    /* force test: number_of_rows_averaged must be larger than number of rows in
       measurement box height */
    if (testing_for_alias == 0) {
      alias_count += alias_and_peak(i,lower_alias_limit,row_array,
                                    alias_count,
                                    ppp,pphp,start_c,end_c,right_col,s_row,
                                    number_of_rows_averaged);
      g_test_pattern[i].rows = save_row_count;
    }
  
    /* this actually does the mtf computations using the array of modulations
       computed from the scan lines within a patch */
    modulation_problems += compute_modulation(i);


    free(g_modulation_array);


    /* Upscaling is reported in the 1's
       Decimation is reported in the 100's */
    /* Check if at least 1 of rows checked for aliasing have upsampling v5.7
       Previous versions checked for 1/3 of rows having aliasing.*/
    j = alias_count % 100;
    if ((j !=0 )) {
      upscaling_detected = 1;
      upscaling_detected_this_pass = 1;
    }

    /* v5.7 Report if at least one decimation case is found.
       Previous versions checked for 1/3 of rows having aliasing.*/
    if (alias_count >= 100) {
      decimation_detected = 1;
      decimation_detected_this_pass = 1;
    }
    if (decimation_detected_this_pass || upscaling_detected_this_pass)
      g_test_pattern[i].star = '*';

	if(g_extended) {
		g_test_pattern[i].totcases = totcases;
		g_test_pattern[i].decim_over = decim_over;
		g_test_pattern[i].ratio_sum = ratio_sum;
		g_test_pattern[i].ratio_over_sum = ratio_over_sum;
	}
  } /* for (i=0;i<g_number_of_patches;i++) */

  if (modulation_problems != 0)
    put_problem("\n", IQS);
  
  /* Print warning when no row averaging is performed. v5.7.5 */
  if (num_non_singleton_averages == 0){
    
    put_problem("caution: \"peak MTF\" may be noisy, \"avg MTF\" may be more accurate\n",0);  
    put_problem("(row averaging noise reduction not performed in peak MTF computation\n",0);  
    put_problem("because skew angle > 5 degrees)\n\n",0);  
  }
  
 /* 
 	 following segment for testing only:
 	 prints g_high_throw_out as determined by various
 	 definitions, and prints # of bogus sine peak or valley values 
  {
    int tmp, n;
    tmp = g_test_pattern[g_white_patch].high;
    n = 0;
    while (tmp > 0) {
      tmp /= 2;
      n++;
    }
    tmp = pow(2.0, (float)n) - 1;

    MTFPRINT5("\ng_high_throw_out = %d  %g  %g  (%d bpp) \n",
	      g_test_pattern[g_white_patch].high,
	      (g_test_pattern[g_white_patch].grey_level + tmp)/2.0,
	      (g_test_pattern[g_white_patch].grey_level + g_max_pixel_value)/2.0,
	      n);
    MTFPRINT2("Bogus_cnt in compute_sines = %d \n", bogus_cnt);
  }
*/ 

    
  REPORT_ALIASING
    
  print_sine_data();

  if (g_extended) {
  
    fprintf(mtftmp,"\n\nAliasing Summary\n");
    fprintf(mtftmp,"\nCy/mm      Cases     Avg      fraction of cases              Avg of cases\n");
    fprintf(mtftmp,  "                     Ratio      > TestValue                   > TestValue\n");
	
    for (i=0;i<g_number_of_patches;i++) {
		if (g_test_pattern[i].patch_type != SINE_PATCH) 
			continue;
			if (g_test_pattern[i].totcases > 0) 
			{
				fprintf(mtftmp,   "%6.3f   %7d    %6.3f             %6.3f", 
					g_test_pattern[i].cycles_per_mm,
					g_test_pattern[i].totcases,
					g_test_pattern[i].ratio_sum/(double)g_test_pattern[i].totcases,
					(double)g_test_pattern[i].decim_over/(double)g_test_pattern[i].totcases
				);
				if (g_test_pattern[i].decim_over > 0)
					fprintf(mtftmp, "                   %6.3f", 
					g_test_pattern[i].ratio_over_sum/(double)g_test_pattern[i].decim_over);
				fprintf(mtftmp, "\n");
			}
    }
	
	fprintf(mtftmp, "\nAll data for secondary ratios (not tertiary)\n");

    fclose(mtftmp);
  }

}


/******************************************************************************/
/*                                                                            */
/******************************************************************************/

#define BIN_TO_FREQUENCY(_bin)  \
 ((float)(_bin) * g_ppi_horz /  \
                  ((float)(fcount/2) * cos(g_skew_avg) * 2.0 * MM_PER_INCH))

/* Check a row of a sine patch for aliasing.  We look for aliasing by taking
   a sample (from the line) that contains a (very close to) whole number of
   periods of the sine wave represented by the patch.  We run a DFT on the
   sample.  If the sample is unaliased the strongest peak should occur in the bin
   that represents the number of periods in the sample (this is latter referred 
   to as the "fundamental" frequency of the sample.) */
/* Return a 1 if upscaling (formerly "aliasing") is detected.
   Return a 100 if decimation is detected. */

int alias_test(int *row_array, double ppp, int fcount)
{
  double out_array[1024];
  int fundamental, folded_fundamental;
  int max_index = 0;
  double max, secondary, a, b, c;
  short j;
  double exact_fund;
  int number_zeros;
  int fundamental_beyond_nyquist;
  int secondary_index = 0;
  /* Only used for debugging */
  int tertiary_index = 0;
  double tertiary = 0.0;

  if (ppp == 0.0)
    return (0);

  fcount = (int)((double)fcount/ppp);

  /* V5.7.5 changed threshold to 2, to prevent problems which
     occur when fundamental = 1 */
  if (fcount < 2) return(0);
  fcount = NINT((double)fcount*ppp);
 
  number_zeros = 0;
  
  DFT(row_array,fcount,number_zeros,out_array,int)
    
  /* our biggest bin should be at "fundamental" */
  exact_fund = (double)fcount/ppp;
  fundamental= NINT(exact_fund);

  /* Find primary and secondary peaks */
  a = b = c = FLT_MAX;        /* This initialization ignores "peaks" at bin 1 */
  max = secondary = 0.0;
  for (j = 2;  j <= fcount;  j += 2)
    {
      a = b;
      b = c;
      c = C_ABS(out_array[j],out_array[j+1]);
      if ((a < b) && (b > c)) {
        /* b is a peak (we assume no plateaus) */
        if (b > max)
          /* b is new max */
          {
	    if (g_extended) {
	      tertiary_index = secondary_index;
	      tertiary = secondary;
	      secondary_index = max_index;
	    }
            secondary = max;
            max_index = (j-2)/2;
            max = b;
          }
        else if (b > secondary)
          /* b is not a new max, but is a new secondary */
          {
	    if (g_extended) {
	      tertiary_index = secondary_index;
	      tertiary = secondary;
	      secondary_index = (j-2)/2;
	    }
            secondary = b;
          }
	else if (g_extended) {
	  if (b > tertiary) {
	    /* b is not a new max nor a new secondary, but is a new tertiary */
            tertiary_index = (j-2)/2;
            tertiary = b;
          }
	}
      }
    }
    
  /* Check for peak at end */
  if (c > b) {

    /* b is a peak (we assume no plateaus) */
    if (c > max)
      /* c is new max */
      {
	if (g_extended) {
	  tertiary_index = secondary_index;
	  tertiary = secondary;
	  secondary_index = max_index;
	}
        secondary = max;
        max_index = (j-2)/2;
        max = c;
      }
    else if (c > secondary)
      /* c is not a new max, but is a new secondary */
      {
	if (g_extended) {
	  tertiary_index = secondary_index;
	  tertiary = secondary;
	  secondary_index = (j-2)/2;
	}
        secondary = c;
      }
    else if (g_extended) {
      if (c > tertiary) {
      /* b is not a new max nor a new secondary, but is a new tertiary */
        tertiary_index = (j-2)/2;
        tertiary = c;
      }
    }
  }


  if (g_extended) {
    totcases++;
	ratio_sum += secondary/max;
    fprintf(mtftmp,"max @ %5.2f   sec @ %5.2f   ratio %.3f",
            BIN_TO_FREQUENCY(max_index),
            BIN_TO_FREQUENCY(secondary_index),
	    (secondary / max));
  }

  /* Check for aliasing.  Return a 1 if aliasing found, 0 if not. */
  /* Upscaling is detected if the DFT peak is not at the correct frequency.
     Decimation is detected if there is a strong secondary DFT peak. */

  /* First check the location of the primary peak */
  fundamental_beyond_nyquist = 0;
  if (ppp >= 2.0)
    /* Fundamental within Nyquist -- strongest peak must be within +/-1 of 
       the fundamental */
    {
      if ((max_index < (fundamental - 1)) || ((fundamental + 1) < max_index))
        {
	  if (g_extended) {
	    fprintf(mtftmp,"    UPSCALING\n              ter @ %5.2f   ratio %.3f      fund. @ %5.2f\n",
		    BIN_TO_FREQUENCY(tertiary_index), (tertiary / max),
                    BIN_TO_FREQUENCY(fundamental));
	  }
          return(1);
        }
    }
  else
    /* Fundamental beyond Nyquist -- it is past the end of the DFT bins and 
       "folds over" back into the bins.  Strongest peak must be within +/-1 
       of this folded frequency. */
    {
      fundamental_beyond_nyquist = 1;
      /* Divide fcount to force truncation--fcount may be odd and we need to 
	 match the results of the "for" loop, above */
      folded_fundamental = fcount/2 - fundamental + fcount/2;
      if ((max_index < (folded_fundamental - 1)) || 
	  ((folded_fundamental + 1) < max_index))
        {
	  if (g_extended) {
	    fprintf(mtftmp,"    UPSCALING\n              ter @ %5.2f   ratio %.3f      folded_fund. @ %5.2f\n",
		    BIN_TO_FREQUENCY(tertiary_index), (tertiary / max),
                    BIN_TO_FREQUENCY(folded_fundamental));
	  }
          return(1);
        }
    }

  /* Primary peak is at the correct frequency.  However, aliasing also occurs
     if the secondary peak is too strong relative to the primary peak. */

  if ((!fundamental_beyond_nyquist && 
       ((secondary / max) > MAX_SECONDARY_TO_PRIMARY_RATIO))
      || ((secondary / max) > MAX_SECPRIM_ABOVE_NYQ))   /* added Apr'04, nbn*/
   {
     if (g_extended) {
	   decim_over++;
	   ratio_over_sum += secondary/max;
       fprintf(mtftmp,"    DECIMATION\n              ter @ %5.2f   ratio %.3f\n",
	       BIN_TO_FREQUENCY(tertiary_index), (tertiary / max));
		   }
     return(100);
   }
  else {
    if (g_extended) 
      fprintf(mtftmp,"\n              ter @ %5.2f   ratio %.3f\n",
	      BIN_TO_FREQUENCY(tertiary_index), (tertiary / max));
    return(0);
  }
}


/******************************************************************************/
/*                                                                            */
/******************************************************************************/

/* Compute modulations and mtf's for a sine patch.  This is done using an array
   of individual modulations for all the peak/valley pairs that were sampled
   for the patch. */
/* Return non-zero if any problems were reported. */
int compute_modulation(unsigned int i)
{
  double sigma,sigmas,mean,xsmm;
  double xsqrd;
  double mac,macmo;
  double max_mod;
  char tf_str[4];

  int j, max_mod_pos=-1;
  char problem_string[82];
  int return_value;

  if(g_bartarget) sprintf(tf_str, "CTF");
  else sprintf(tf_str, "MTF");

  return_value = 0;
  if (g_mod_array_cnt) {
    mean = 0.0;
    xsqrd =  0.0;
    max_mod = -99999.99;
    for (j=0;j<g_mod_array_cnt;j++) {
      mean += g_modulation_array[j];
      if (g_modulation_array[j] > max_mod) {
        max_mod = g_modulation_array[j];
	max_mod_pos = j;
      }
      xsqrd += g_modulation_array[j] * g_modulation_array[j];
    }
    mac = (double)g_mod_array_cnt;
    macmo = mac-1.0;
    mean /= mac;
    xsmm = xsqrd/macmo - mac*mean*mean/macmo;
    if (xsmm > 0.0) sigma = sqrt(xsmm);
    else sigma = 0.0;
    sigmas = sigma*(double)NO_SINE_SIGMAS;
    
    g_test_pattern[i].first_mean = mean;
    g_test_pattern[i].std_deviation = sigma;
    g_test_pattern[i].computed_modulation_average_plus =mean+sigmas;
    g_test_pattern[i].computed_modulation_average = mean;
    g_test_pattern[i].computed_modulation_highest = max_mod;
    g_test_pattern[i].fliers=0;
    for (j=0;j<g_mod_array_cnt;j++)
      if ((g_modulation_array[j] - mean) > sigmas)
        g_test_pattern[i].fliers++;
  }
  else {
    g_test_pattern[i].computed_modulation_highest = 0;
  }
  g_test_pattern[i].mtf_highest =
    g_test_pattern[i].computed_modulation_highest /
      g_test_pattern[i].modulation;
  g_test_pattern[i].mtf_average =
    g_test_pattern[i].computed_modulation_average /
      g_test_pattern[i].modulation;

  g_test_pattern[i].sine_peak = g_mod_peak_array[max_mod_pos];
  g_test_pattern[i].sine_valley = g_mod_valley_array[max_mod_pos];
  
  /*
  MTFPRINT10("SinePatch %d:  first_mean=%g  std=%g   mod_avg_plus=%g\n mod_avg=%g mod_high=%g fliers=%d\n mtf_high=%g mtf_avg=%g\n\n", 
	   i,
	   g_test_pattern[i].first_mean,
	   g_test_pattern[i].std_deviation,
	   g_test_pattern[i].computed_modulation_average_plus,
	   g_test_pattern[i].computed_modulation_average,
	   g_test_pattern[i].computed_modulation_highest,
	   g_test_pattern[i].fliers,
	   g_test_pattern[i].mtf_highest,
	   g_test_pattern[i].mtf_average)
  */

  if (!g_printer) {
    double freq_cutoff = (double)EFF_FREQ_CUTOFF*PI/180.0; 
    double abs_skew_avg = fabs(g_skew_avg); 

    if ((abs_skew_avg <= freq_cutoff) && 
	(g_test_pattern[i].mtf_highest > g_test_pattern[i].mtful ||
	 g_test_pattern[i].mtf_highest < g_test_pattern[i].mtfll)) 
      {
	sprintf(problem_string,
		"Computed %s at %.3f cy/mm is outside spec** range of %.3f - %.2f\n",
		tf_str, 
		g_test_pattern[i].effective_frequency,
		g_test_pattern[i].mtfll,
		g_test_pattern[i].mtful);
	put_problem(problem_string, IQS);
	return_value = 1;
	need_spec_footnote = g_target_res;
      }
  }
  
  return(return_value);
}


/******************************************************************************/
/*                                                                            */
/******************************************************************************/

/* output to the report */

void print_sine_data(void)
{
  unsigned int i;

  char tf_str[4];
  char freq_type[10];
  
  if(g_digital_camera) sprintf(freq_type,"Image    ");
  else                 sprintf(freq_type,"Effective");

  if(g_bartarget) sprintf(tf_str, "CTF");
  else sprintf(tf_str, "MTF");
  
  MTFPRINT("\n")
  if (g_print_average_mtf) {
    MTFPRINT4("Tgt Freq.   %s     Target        peak %s       avg %s       Std Dev\n", freq_type, tf_str, tf_str)
    MTFPRINT2("cy/mm       Freq.         Mod                                        avg %s\n", tf_str)

    for (i=0;i<g_number_of_patches;i++) {
      if (g_test_pattern[i].patch_type == SINE_PATCH &&
          g_test_pattern[i].computed_modulation_highest !=
          BEYOND_NYQUIST) {
        MTFPRINT9("%6.3f%12.3f%13.3f%15.3f%1c%14.3f%1c%14.3f\n",
                  g_test_pattern[i].cycles_per_mm,
                  g_test_pattern[i].effective_frequency,
                  g_test_pattern[i].modulation,
                  g_test_pattern[i].mtf_highest,
                  g_test_pattern[i].star,
                  g_test_pattern[i].mtf_average,
                  g_test_pattern[i].star,
                  g_test_pattern[i].std_deviation)
      }
      else if (g_test_pattern[i].computed_modulation_highest ==
               BEYOND_NYQUIST)
        MTFPRINT2("%6.3f\t\tSampled Beyond Nyquist Frequency\n", g_test_pattern[i].cycles_per_mm)
    }
  }
  else {
    MTFPRINT3("Tgt Freq.   %s     Target       peak %s\n", freq_type, tf_str)
    MTFPRINT2("%s\n","cy/mm       Freq.         Mod")

    for (i=0;i<g_number_of_patches;i++) {
      if (g_test_pattern[i].patch_type == SINE_PATCH &&
          g_test_pattern[i].computed_modulation_highest !=
          BEYOND_NYQUIST) {
        MTFPRINT6("%6.3f%12.3f%13.3f%15.3f%1c\n",
                  g_test_pattern[i].cycles_per_mm,
                  g_test_pattern[i].effective_frequency,
                  g_test_pattern[i].modulation,
                  g_test_pattern[i].mtf_highest,
                  g_test_pattern[i].star)
      }
      else if (g_test_pattern[i].computed_modulation_highest ==
               BEYOND_NYQUIST)
        MTFPRINT2("%6.3f\t\tSampled Beyond Nyquist Frequency\n", g_test_pattern[i].cycles_per_mm)
    }
  }

  if (g_printer)
    MTFPRINT4("Printer %s = (%s output of this run)/(%s of scanner used to scan print)\n", tf_str, tf_str, tf_str);

}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/

/* if the -e option is selected */

void print_sine_stats(void)
{
  unsigned int i;
  char buffer[101];
  char freq_type[10];
  
  if (g_digital_camera) sprintf(freq_type,"    image");
  else					sprintf(freq_type,"effective");

  fprintf(g_mtfout, "\n\n\n\tSINE PATCHES: IMAGE GRAY LEVEL STATISTICS\n\n");
  
  fprintf(g_mtfout, "* Sine peak & adjacent valley selected for Peak MTF\n\n");

  for (i=0;i<g_number_of_patches;i++) {
    if (g_test_pattern[i].patch_type == SINE_PATCH) {
      fprintf(g_mtfout, "Tgt cy/mm: %f     %s cy/mm = %f\n",
	      g_test_pattern[i].cycles_per_mm, freq_type, 
	      g_test_pattern[i].effective_frequency);
      
      if (g_test_pattern[i].computed_modulation_highest == BEYOND_NYQUIST) 
	continue;
      fprintf(g_mtfout, "rows averaged for peak MTF: %d             sinepatch: %s\n",
                g_test_pattern[i].rows, g_test_pattern[i].patch_label);
      fprintf(g_mtfout, "Avg Image Mod = %f     Std Dev = %f",
	      g_test_pattern[i].first_mean,
	      g_test_pattern[i].std_deviation);
      
      print_histo(i);
    }
  }

  mtftmp = fopen("_MTF_TMP", "r");
  fprintf(g_mtfout, "\nALIASING REPORT\n");

  fprintf(g_mtfout, "Aliasing from Decimation if: (sec/max) > TestValue\n");
  fprintf(g_mtfout," TestValue below Nyquist = %f\n", MAX_SECONDARY_TO_PRIMARY_RATIO);
  fprintf(g_mtfout," TestValue above Nyquist = %f\n", MAX_SECPRIM_ABOVE_NYQ);
  fprintf(g_mtfout,"Aliasing from Upscaling if max frequency is not expected frequency\n\n");

  while (fgets(buffer, 100, mtftmp) != NULL) {
    fprintf(g_mtfout,buffer);
  }
  fclose(mtftmp);
  fprintf(g_mtfout,"\n");

  /* Now remove the temporary file */
  remove("_MTF_TMP");

  
}

