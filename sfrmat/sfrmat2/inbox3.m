function [del, weights] = inbox3;
% [del, weights] = inbox3 dialogue box for input of data sampling, weights
% for red, green and blue signals for luminanace calculation for SFR 
% calculation
%  Usage: [del, weights] = inbox3
%  del     = sampling interval in mm
%  weights = 3x1 array of luminance weights  for R, G, B
%             defaults: 0.3  0.6  0.1
% Calls inputdlg function, supplied with the toolbox
% matlab/uitools. If you have problems, check which version of
% inputdlg.m (or corresponding inputdlg.p) is being called. You need
% version 1.48 or later.
% Peter Burns 19 April 2000
% Copyright (c) International Imaging Industry Association 1999 - 2003

title='  Data sampling & weights ';
prompt={'Data sampling in dpi',' or mm', 'Luminance weights  for R, G, B' };
   def={'-', '1', ['0.3'
                     '0.6'
                     '0.1']};
lineNo=[1, 1, 3]';
AddOpts.Resize='on';
AddOpts.WindowStyle='normal';
AddOpts.Interpreter='tex';
answer=inputdlg(prompt, title, lineNo, def, AddOpts);

sflag = 0;
if length(char(answer(1)))~=1;
 sflag = 1;
 elseif char(answer(1))~='-';
 sflag = 1;
end;
if sflag==0
  del =  str2num(char(answer(2)));  
 else;
  del =  str2num(char(answer(1)));
  del = 25.4/del;
end; 

weights = str2num(char(answer(3)));
if sum(weights, 1) > 1.0;
 disp(' ***  WARNING: Sum of Luminance weights is greater than 1  ***');
end;
