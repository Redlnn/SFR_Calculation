function [status, dat, fitme, esf] = sfrmat2(io, a, del);
%  sfrmat2:  Slanted-edge and color mis-registration analysis
%  [status, dat, fitme, esf] = sfrmat2(io, a, del);
%       From a selected edge area of an image, the program computes
%       the ISO slanted edge MTF. Input file can be single or
%       three-record file. Format can be JPEG or TIF. The image is 
%       displayed and a region of interest (ROI) can be chosen, or
%       the entire field will be selected by not moving the mouse
%       when defining an ROI. Either a vertical or horizontal edge
%       feature can be selected. 
%  Input arguments:
%       io (optional) integer [1-4] controls computation, and 
%       reporting of output
%       a (optional) an nxm or nxmx3 array of data     (used with io=4)
%       del (optional) sampling interval in mm or dpi  (used with io=4)
%          If dx < 1 it is assumes to be sampling pitch in mm
%          If io = 4 (see below, no GUI) and del is not specified,
%          it is set
%          equal to 1, so frequency is given in cy/pixel.
%  Returns: 
%       status = 0 if normal execution
%       dat = computed sfr data in 2 (freq, sfr) or 5 columns.
%       fitme = coefficients for the linear equations for the fit to
%               edge locations for each color-record
%       esf   = supersampled edge spread function vector
%  Usage:
%       
%       sfrmat2(0) = ISO Plugin type output (no color registration,
%                     single luminance SFR), with prompts for I/O 
%       sfrmat2(1) = ISO output + edge location(s)
%
%       sfrmat2  
%       sfrmat2(2) = (default) R,G,B,Lum SFRs + edge location(s)
%
%       sfrmat2(3) = only edge location(s)
%       [status, dat, fitme] = sfrmat2(4, array, a, del) 
%                  returns SFRs and edge locations to the variables
%                  [[status, dat, fitme] without plotting or other
%                  output. Input 'array' is an nxm or nxmx3 array of
%                  data. No Graphical User Interface (GUI) is used with
%                  this option.
%  Needs: ahamming, cent, centroid, clipping, deriv1, findedge,
%         getoecf, getroi, inbox1, inbox3, inputdlg, project, results,
%         rotatev, splash, isarray
% 
%  Author: Peter Burns, peter.burns@kodak.com
%          12 August 2003
%  Copyright (c) International Imaging Industry Association
%******************************************************************
status = 0;
defpath = path;            % save original path
home = pwd;                % add current directory to path                   
addpath([home]);
name =    'sfrmat';
version = '2.1';
when =    '12 AUGUST 2003';

%close all;          %close all figure window before starting (disabled)
% Suppresses interpreting of e.g. filenames
 set(0, 'DefaultTextInterpreter', 'none'); 
pflag = 0;
if nargin==0;
 io = 2;   % default option
elseif nargin>3;
 disp('Incorrect number or arguments. There should be 0-3');
 status = 1;
 return;
if nargin<3
 del = 1;
end
elseif io<0;
 disp('Bad input arguments, first argument, io, must be  1, 2, 3, 4');
 status = 1;
 return;
elseif io>4;  %%
 disp('Bad input arguments, first argument, io, must be  1, 2, 3, 4');
 status = 1;
 return;
elseif io==4;
 if nargin<2;
 disp('Bad input arguments, if io=4 must you must supply (second argument) data');
 status = 1;
 return;
 end;
end;

%% 

if io ~= 4;
  swin = splash(name, version, when);
  % Select file name for reading
  % edit the next line to change the default path for input file selection
  def ='*.*';

 [filename, pathname] = uigetfile(def,'Select input data file (JPEG, TIF, BMP, HDF, PCX, XWD)');
 close(swin);
 if (filename==0)
  disp('No file selected. To try again type: > sfrmat2');
  status = 1;
  return;
 end;
 filename=[pathname,filename];
 readstat = 0;
 % Check for common extensions for JPEG, TIFF images. 
 % This may not be needed since IMREAD does this sort of guessing
 temp = filename(length(filename)-2:length(filename));   %last 3 characters
                            
 if           temp=='iff'
          fmt = 'tif';
       elseif temp=='IFF'
          fmt = 'tif';
       elseif temp =='TIF'
          fmt = 'tif';
       elseif temp =='tif'
          fmt = temp;
       elseif temp=='JPG'
          fmt='jpg';
       elseif temp=='PEG'
          fmt='jpg';
       elseif temp=='jpg'
          fmt=temp;
       elseif temp=='peg'
          fmt='jpg';
       elseif temp=='bmp'
          fmt=temp;
       elseif temp=='BMP'
          fmt='bmp';
       elseif temp=='hdf'
          fmt=temp;
       elseif temp=='HDF'
          fmt='hdf';
       elseif temp=='pcx'
          fmt=temp;
       elseif temp=='PCX'
          fmt='pcx';
       elseif temp=='xwd'
          fmt=temp;
       elseif temp=='XWD'
          fmt='xwd';
       else
        disp(' ');
        disp(' ********  WARNING: Unknown file format  **************'); 
        disp(' sfrmat programme looks for extensions: tif, tiff, jpeg');
        disp(' jpg, bmp, hdf, pcx, xwd, TIF, TIFF, JPEG, JPG, BMP, ...');
        disp(' Matlab function IMREAD will infer format and attempt to');
        disp(' read. If this fails, type "help imread" for info. on ');
        disp(' supported formats.');
        readstat = 1;
       return;
 end;

 if readstat ==0;
  atemp = imread(filename, fmt);  
 else;
  atemp = imread(filename);
 end;

[nlin npix ncol] = size(atemp);
 
% input sampling and luminance weights
 if ncol==1;
	 del = inbox1;
 else [del, weight] = inbox3; 
 end;

cname = class(atemp);
if cname(1:5) == 'uint1';   % uint16
 smax = 2^16-1;
elseif cname(1:5) == 'uint8';
 smax = 255;
else
 smax = 1e10;
end

[a, roi] = getroi(atemp);               % extract Region of interest
clear atemp                             % *******************************
[nlow, nhigh, cstatus] = clipping(a, 0, smax, 0.005);
if cstatus ~=1; 
 disp('Fraction low data');
 disp(nlow);
 disp('Fraction high data');
 disp(nhigh);
end;

oename = 0;
[oename, oepath] = uigetfile([pathname,'*.*'],['Select OECF data file, ',num2str(ncol),' record(s)']);

if oename==0;
 disp(' No OECF transformation chosen');
 oename='none';
 else;
 
 [a, oestatus] = getoecf(a, oepath, oename);   % Transforms a using OECF LUT from file chosen
 if oestatus ~=0;
  disp('');
  disp('*  Number of data colors does not match OECF LUT file,');
  disp('*  Care to try again?');
  [oename, oepath] = uigetfile([oepath,'*.*'],'Select other OECF data file or CANCEL');
 
  if oename==0;
   disp(' No OECF transformation chosen');
   oename='none';
   else;
   %oename = [oepath,oename];	
   [a, oestatus] = getOECF(a, oepath, oename);   % Transforms a using OECF LUT from file chosen
   if oestatus ~=0;
    status = 1;
    return;
   end;
  end;
 end;
end;
else;                     % when io = 4
  a= double(a);
%                  % default sampling and color weights
 if del > 1
  del = 25.4/del;  % Assume input was in DPI convert to pitch in mm
 end;
  weight = [0.3
            0.6 
           0.1];
end; 
[nlin npix ncol] = size(a);

% Form luminance record using the weight vector for red, green and blue
if ncol ==3;
 lum = zeros(nlin, npix);
 for i=1:nlin;
 for j=1:npix;
 lum(i,j) = weight(1,1)*a(i,j,1)+weight(2,1)*a(i,j,2)+weight(3,1)*a(i,j,3);
 end;
 end;
  cc = zeros(nlin, npix, 4);
  for i=1:nlin;
  for j=1:npix;
   cc(i,j,:) = [ a(i, j, 1), a(i, j, 2), a(i, j, 3), lum(i,j)];
  end;
  end;
  a = cc;
  clear cc;
  clear lum;
  ncol = 4; 
end;
% rotate horizontal edge so it is vertical
 [a, nlin, npix, rflag] = rotatev(a);

loc = zeros(ncol, nlin);

fil1 = [0.5 -0.5];
fil2 = [0.5 0 -0.5];
% Need 'positive' edge for good centroid calculation
 tleft  = sum(sum(a(:,      1:5,  1),2));
 tright = sum(sum(a(:, npix-5:npix,1),2));
  if tleft>tright;
   fil1 = [-0.5 0.5];
   fil2 = [-0.5 0 0.5];
  end
% Test for low contrast edge;
 test = abs( (tleft-tright)/(tleft+tright) );
 if test < 0.2;
  disp(' ** WARNING: Edge contrast is less that 20%, this can');
  disp('             lead to high error in the SFR measurement.');
 end; 

fitme = zeros(ncol, 2);
slout = zeros(ncol, 1);

% smoothing window for first part of edge location estimation - 
%  to used on each line of ROI
 win1 = ahamming(npix, (npix+1)/2);    % Symmetric window

for color=1:ncol;                      % Loop for each color

 %
 pname = ' ';
  if ncol~=1;
   pname =[' Red '
           'Green'
           'Blue '
           ' Lum '];
  end;

%%
 c = deriv1(a(:,:,color), nlin, npix, fil1);

% compute centroid for derivative array for each line in ROI. NOTE WINDOW array 'win'
 for n=1:nlin
  loc(color, n) = centroid( c(n, 1:npix )'.*win1) - 0.5;   % -0.5 shift for FIR phase
 end;
% clear c

 fitme(color,:) = findedge(loc(color,:), nlin);
 place = zeros(nlin,1);
 for n=1:nlin;
  place(n) = fitme(color,2) + fitme(color,1)*n;
  win2 = ahamming(npix, place(n));
 loc(color, n) = centroid( c(n, 1:npix )'.*win2) -0.5;
 end;

 fitme(color,:) = findedge(loc(color,:), nlin);
 fitme(color,:); % used previously to list fit equations

%%

end;                               % End of loop for each color

disp('  ');

summary{1} = [' ']; % initialize

if io > 0;

 ncolout = ncol;                   % output edge location listing
 if ncol == 4;
  ncolout = ncol - 1;
 end

 midloc = zeros(ncolout,1);
 summary{1} = ['Edge location, slope']; % initialize

 for i=1:ncolout;
  slout(i) = - 1./fitme(i,1);      % slope is as normally defined in image coods.
  if rflag==1,                     % positive flag it ROI was rotated
   slout(i) =  - fitme(i,1);
  end;

% evaluate equation(s) at the middle line as edge location
   midloc(i) = fitme(i,2) + fitme(i,1)*((nlin-1)/2);

  summary{i+1} = [midloc(i), slout(i)];
 end

disp('Edge location(s) and slopes = ' ), disp( [midloc(1:ncolout), slout(1:ncolout)]); 
 if ncol>2;
  summary{1} = ['Edge location, slope, misregistration (second record, G, is reference)'];
  misreg = zeros(ncolout,1);
  for i=1:ncolout;
   misreg(i) = midloc(i) - midloc(2);
   summary{i+1}=[midloc(i), slout(i), misreg(i)];
  end;

 disp('Misregistration, with green as reference (R, G, B, Lum) = ');
 for i = 1:ncolout
  fprintf('%10.4f\n', misreg(i));
 end;
 end

end                            % end of check if io > 0

% Full linear fit is available as variable fitme. Note that the fit is for
% the projection onto the X-axis,
%       x = fitme(color, 1) y + fitme(color, 2)
% so the slope is the inverse of the one that you may expect


% *********************************** Stop before sfr calculations here
if io ==3;
 dat = 0;                      % avoids error message
 return;
end;

%%

nbin = 4;

nn =   floor(npix *nbin);
mtf =  zeros(nn, ncol);
nn2 =  nn/2 + 1;

%%

freq = zeros(nn, 1);
for n=1:nn;
  freq(n) = nbin*(n-1)/(del*nn);
end;

freqlim = 1;
 % limits plotted sfr to 0- 1 cy/pxel freqlim = 2 for all data

nn2out = round(nn2*freqlim/2);

nfreq = n/(2*del*nn);    % half-sampling frequency

win = ahamming(nbin*npix,(nbin*npix+1)/2);      % centered Hamming window

%%
% **************                      Large SFR loop for each color record

esf = zeros(nn,ncol);

for color=1:ncol
% project and bin data in 4x sampled array
 [point, pstatusp] = project(a(:,:,color), loc(color, 1), fitme(color,1), nbin);

 esf(:,color) = point;  

% compute first derivative via FIR (1x3) filter fil
  c = deriv1(point', 1, nn, fil2);
  c = c';
 mid = centroid(c);
  temp = cent(c, round(mid));              % shift array so it is centered
  c = temp;
  clear temp;

% apply window (symmetric Hamming)
 c = win.*c;

% Transform, scale %%

 temp = abs(fft(c, nn));
 mtf(1:nn2, color) = temp(1:nn2)/temp(1);
%% 
end;

dat = zeros(nn2out, ncol+1);
for i=1:nn2;
 dat(i,:) = [freq(i), mtf(i,:)];
end;

if io ==4;         
 return
end

% Plot SFRs on same axes
if ncol >1;
  sym{1} = []; 
  sym{1} = '--r';
  sym{2} = '-g';
  sym{3} = '-.b';
  sym{4} = '*k';
  ttext = [filename,' r = - -  g = -  b = -.-  lum = *'];
 else;
  ttext = filename;
  sym{1} = 'k';
end

screen = get(0, 'ScreenSize');

 defpos = get(0, 'DefaultFigurePosition');
 set(0, 'DefaultFigurePosition', [15 25 0.6*screen(3) 0.4*screen(4)]);
   
 
if del==1;
   funit =  'cy/pixel';
else funit = 'cy/mm';
end;

 plot( freq( 1:nn2out), mtf(1:nn2out, 1), sym{1});
 hold on;
  title(ttext);
  xlabel(['     Frequency, ', funit]);
  ylabel('SFR');
  if ncol>1;
   for n = 2:ncol;
   plot( freq( 1:nn2out), mtf(1:nn2out, n), sym{n});
  end;
end;

line([nfreq ,nfreq],[.05,0]),
text(.95*nfreq,+.08,'half-sampling'),

 hold off;
zoom on;


defname = [pathname,'*.*'];
   [outfile,outpath]=uiputfile(defname,'File name to save results');
   foutfile=[outpath,outfile];
   if size(foutfile)==[1,2],
      if foutfile==[0,0],
         disp('Saving results: Cancelled')
         return;
      end;
  
    else;
    results(dat, filename, roi, oename, summary, foutfile);
   end;

% Clean up
% Reset figure position and size
  set(0, 'DefaultFigurePosition', defpos);
% Reset text interpretation
  set(0, 'DefaultTextInterpreter', 'tex')
  path(defpath);             % Restore path to previous list
  cd([home]);                % Return to working directory
 
disp(' * sfrmat2 finished  *');

return;
