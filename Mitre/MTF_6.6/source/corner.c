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

/******************************************************************************/
/*									      */
/* Revisions to corner.c						      */
/*									      */
/* V5.0 Original generated to refine corner point positions                   */
/*                                                     mal, 10/03             */
/*                                                                            */
/* V5.5     Cleaned up compiler warnings.                                     */
/*                                                     mal, 01/04             */
/*                                                                            */
/******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
	
typedef void TIFF;
#include "mtf.h"
#include "xterndcl.h"


/******************************************************************************/
/*									      */
/******************************************************************************/
/* Determine col(row) that contains the most edge pixels */
int hist_edge(int *edges, int *hist, int width, int length)
{
  int i, mindex, max;
  int number_of_points;

  number_of_points=0;
  mindex=-1;

  /* Zero out the histogram space */
  for (i=0;i<length; i++) hist[i] = 0;

  /* Histogram edge position data */
  for (i=0;i<width; i++) {
    if (edges[i] > 0) {
      number_of_points++;
      hist[edges[i]]++;
    }
  }
  
  /* Find histogram max value and position */
  max = -1;
  for (i=0;i<length; i++) {
    if ( max <= hist[i] ) {
      max = hist[i];
      mindex = i;
    }
  }
 
  return (mindex);
  
}

/******************************************************************************/
/*									      */
/******************************************************************************/
/* See if current max position is significantly larger than 5 previous maxima. 
   If so, it's a significant edge, if not it's noise. */
int check_preceeding_stats(int *buf, int mi, int max) 
{
  int cnt, high, i;
  int *bptr;
  
  if( mi > 5 ) {
    /* Check local stats to see if max2 is indeed edge-like */
    cnt = 0;
    bptr = buf+mi-4;
    high = 0;
    for (i=mi-4; i>=0 && cnt<5 ; i--, bptr--){	
      if( abs(*bptr) > high ) { 
	high = abs(*bptr); 
	cnt++;
      }
    }

    if (high < 3.0*max/4) /* This is a true strong edge maximum, so keep it  */
      return 1;
  }
  return 0;
}

/******************************************************************************/
/*									      */
/******************************************************************************/
/* Locate first significant edge in a line of pixels */
int find_first_edge (int *line, int *buf, int n) 
{
  int diff, i, max, max2, mi, mi2;
  int *ptr1, *ptr, *bptr;
  int samesign, valid;

  /* Compute difference at neighboring pixels.  */
  ptr = line;
  ptr1 = line+1;
  bptr = buf;
  for(i=0; i<n-1; i++) {
    diff = *ptr++ - *ptr1++;
    *bptr++ = diff;
  }
  *bptr=0; /* Set last difference (unknown) to zero */
  
  /* Scan diff array supressing any values that are not 
     local absolute maxima. Record the two highest local maxima. */
  ptr = buf;
  bptr = buf+1;
  for(i=0; i<n-1; i++) {
    samesign = (*ptr * *bptr>=0)?1:0;
    if (samesign && (abs(*ptr) < abs(*bptr))) *ptr=0;
    ptr++;
    bptr++;
  }
  ptr--;
  bptr--;
  i--;
  max = max2 = 0;
  mi = mi2 = -1;
  for(; i>=0; i--) {
    samesign = (*ptr * *bptr>=0)?1:0;
    if (samesign && (abs(*ptr) > abs(*bptr))) *bptr=0;
    else {
      if(abs(*bptr) > max) {	
	max2 = max; mi2=mi;
	max = abs(*bptr); mi=i+1;
      }
      else if(abs(*bptr) > max2) { max2 = abs(*bptr); mi2=i+1; }
    }
    ptr--;
    bptr--;
  }
  if (abs(*bptr)>max) {
	max2 = max; mi2=mi;
	max = abs(*bptr); mi=0;
  }
  else if(abs(*bptr) > max2) { max2 = abs(*bptr); mi2=0; }

  /* Check the two maxima candidates to find best */

  /* If largest is left-most then discard second largest */
  if (mi < mi2) {
    /* Mi2 is not a valid candidate */
    mi2 = mi;
    max2 = max;
  }

  /* Check the stats preceeding two maxima until a significant edge
     is found */
  valid = check_preceeding_stats(buf, mi2, max2);
  if( valid ) {
    mi = mi2;
    max=max2;
  }
  else 
    valid = check_preceeding_stats(buf, mi, max);

  /* If a significant edge was found, then move back towards the left,
     to find area where data becomes 'flat'. */
  if (valid) {
    while(mi>0 && abs(line[mi] - line[mi-1]) >= max/2) mi--;
    return (mi);
  }
  
  else return(-1);
  
}

/******************************************************************************/
/*									      */
/******************************************************************************/
/* Data array contains small portion of the image centered on the input corner
coordinates, oriented so that the panel of interest is to the upper left. 
Search for location of the left/uppermost significant corner */
int find_corner(int **data, int width, int height, int num, int *nx, int *ny)
{
  int i,j;
  int *vedge, *hedge;
  int *rowvavg, *colhavg, *rowptr, *colptr;
  int *buf;
  int maxv, maxh, length, tmp;
  int *cptr, *c2ptr;
    
  vedge = (int *) calloc(height,sizeof(int));
  hedge = (int *) calloc(width,sizeof(int));
  
  rowvavg = (int *) calloc(width,sizeof(int));
  colhavg = (int *) calloc(height,sizeof(int));

  length = (height < width) ? width: height;
  buf = (int *) calloc(length,sizeof(int));

  /* Average parallel to anticipated edge direction */
  /* (vertical or horizontal) */
  /* Then find pixel just prior to first significant edge */

  /* First vertical edge pixel position */
  cptr = data[0];
  for (i=0; i<5; i++){
    rowptr = rowvavg;
    for (j=0; j<width; j++)
      *rowptr++ += *cptr++;
  }
  vedge[2] = find_first_edge (rowvavg, buf, width);

  /* First horizontal edge pixel position */
  for (i=0; i<height; i++) {
    cptr = data[i];
    for (j=0; j<5; j++)
      colhavg[i] += *cptr++;
  }
  hedge[2] = find_first_edge (colhavg, buf, height);

  /* Rest of the vertical edge pixel positions */
  cptr = data[0];
  c2ptr = data[5];
  for (i=3; i<height-2; i++) {
    rowptr = rowvavg;
    for (j=0; j<width; j++)
      *rowptr++ += *c2ptr++ - *cptr++;
    vedge[i] = find_first_edge (rowvavg, buf, width);
  }
  
  /* Rest of the horizontal edge pixel positions */
  for (j=3; j<width-2; j++) {
    colptr = colhavg;
    for (i=0; i<height; i++)
      *colptr++ += data[i][j+2]-data[i][j-3];
    hedge[j] = find_first_edge (colhavg, buf, height);
  }

  /* Determine how far edge info extends in each direction */
  maxv = vedge[2];
  for (i=3; i<height-2; i++) 
      if (vedge[i] > maxv) maxv = vedge[i];
  
  maxh=hedge[2];
  for (i=3; (i<=maxv)&&(i<width-2); i++) 
      if (hedge[i] > maxh) maxh = hedge[i];
  
  /* Find the predominant edge position in each direction */
  *nx = hist_edge(vedge, buf, maxh, length);
  *ny = hist_edge(hedge, buf, *nx, length);
  if (*nx < 5 || *ny < 5) return 0;
  *nx = hist_edge(vedge, buf, *ny, length);
  if (*nx < 5) return 0;
  if (*ny > 12) {
    tmp = hist_edge(vedge+*ny-10, buf, 10, length);
    if (abs(*nx - tmp) < 3) *nx = tmp;
  }
  if (*nx > 12) {
    tmp = hist_edge(hedge+*nx-10, buf, 10, length);
    if (abs(*ny - tmp) < 3) *ny = tmp;
  }

  free(hedge);
  free(vedge);
  free(buf);

  return 1;
  
}

/******************************************************************************/
/*									      */
/******************************************************************************/
/* This is done on the non-normalized, file storage version of the image.
   We search a box of width 30% density patch width, centered at the user 
   corner entry. Find the ideal corner.  */

void refine_corner_coords(TIFF *tif)
{
  int i,j,k, cnt, outi;
  unsigned char *buf;
  int **data;
  int width, height;
  int newx, newy;
  int valid;

  int search_width, search_height;
  float corner_xest[3];
  float corner_yest[3];

  /* Only search area 30% size of actual patch.
     User needs to give approximate value at least this close.
  */
  search_width = NINT(g_ppi_horz*0.3*g_patch_width);
  search_height= NINT(g_ppi_vert*0.3*g_patch_height);

  /* Number to corners to aid later loop */
  corner_xest[0] = NINT(g_input_test_pattern_xul);
  corner_yest[0] = NINT(g_input_test_pattern_yul);
  
  corner_xest[1] = NINT(g_input_test_pattern_xur);
  corner_yest[1] = NINT(g_input_test_pattern_yur);
  
  corner_xest[2] = NINT(g_input_test_pattern_xll);
  corner_yest[2] = NINT(g_input_test_pattern_yll);

  /* Add extra 5 rows/cols of uniform data */
  width = search_width+5;
  height = search_height+5;

  data = (int **) malloc(height*sizeof(int *));
  data[0] = (int *) malloc(width*height*sizeof(int));
  for (i=1; i<height; i++) 
    data[i] = data[i-1] + width;

  buf = (unsigned char *) malloc(g_scan_cols*sizeof(short));

  /* For each corner, read in neighboring image data to a fixed
     orientation.  (Outer corner patch is always to upper left.) */
  for (k=0; k<3; k++) {
    int begx,endx, begy,endy;
    
    begx=corner_xest[k]-search_width/2;
    begy=corner_yest[k]-search_height/2;

    if (k==0 || k==2)  begx -= 5;
    if (k==0 || k==1)  begy -= 5;

    endx = width+begx-1;
    endy = height+begy-1;

    /* Read either 8-bit or 16-bit data into integer data array */
    switch( g_bytes_per_pixel) {
      unsigned short *sbuf;
      
    case 1:
      for (i=begy, cnt=0 ; cnt<height; cnt++, i++) {
	read_scan_line(tif,buf,(short)i);
	if (k==2) outi = endy-i;
	else outi = cnt;
	if (k!=1)
	  for(j=begx; j<=endx; j++)
	    data[outi][j-begx]=buf[j];
	else
	  for(j=begx; j<=endx; j++)
	    data[outi][endx-j]=buf[j];
      }
      break;

    case 2:
      sbuf = (unsigned short *) buf;
      for (i=begy, cnt=0 ; cnt<height; cnt++, i++) {
	read_scan_line(tif,buf, (short)i);
	if (k==2) outi = endy-i;
	else outi = cnt;
	if (k!=1)
	  for(j=begx; j<=endx; j++)
	    data[outi][j-begx]=sbuf[j];
	else
	  for(j=begx; j<=endx; j++)
	    data[outi][endx-j]=sbuf[j];
      }
      break;
    }

    /* Determine relative automatic corner location */
    valid = find_corner(data, width, height, k, &newx, &newy);

    if (valid) {
      /* Get absolute coordinates for new corners */
      corner_xest[k] = begx + newx+0.5;
      corner_yest[k] = begy + newy+0.5;
      if (k==1) corner_xest[k] = begx + (width - newx-0.5);
      if (k==2) corner_yest[k] = begy + (height - newy-0.5);
    }
    else { /* Something went wrong, report the error */
      if(k==0)
	g_UL_auto_crnr_err=1;
      else if(k==1)
	g_UR_auto_crnr_err=1;
      else 
	g_LL_auto_crnr_err=1;
    }
  }

  free(data[0]);
  free(data);

  /* Reset the corner position info */
  g_input_test_pattern_xul = corner_xest[0] ;
  g_input_test_pattern_yul = corner_yest[0] ;
  g_input_test_pattern_xur = corner_xest[1] ;
  g_input_test_pattern_yur = corner_yest[1] ;
  g_input_test_pattern_xll = corner_xest[2] ;
  g_input_test_pattern_yll = corner_yest[2] ; 

}

