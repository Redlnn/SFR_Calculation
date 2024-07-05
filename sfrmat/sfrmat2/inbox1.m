function [del] = inbox1;
% [del] = inbox1 dialogue box for data sampling for SFR calculation
%  Usage: [del] = inbox1
%  del = sampling interval in mm
% Calls inputdlg function, supplied with the toolbox
% matlab/uitools. If you have problems, check which version of
% inputdlg.m (or corresponding inputdlg.p) is being called. You need
% version 1.48 or later.
% Peter Burns 5 Aug. 2002
% Copyright (c) International Imaging Industry Association

title='  Data sampling ';
prompt={'Data sampling in dpi',' or mm' };
   def={'-', '1'};
lineNo=[1];
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
