function [array, status] = getoecf(array, oepath,oename);
% [array, status] = getoecf(array, oepath,oename)  Read and apply oecf
% Reads look-up table and applies it to a data array
%   array = data array (nlin, pnix, ncolor)
%   oepath = table pathname, e.g. /home/sfr/dat
%   oename = tab-delimited text file for table (256x1, 256,3)
%   array = returns transformed array
%   status = 0  OK, 
%          = 1 bad table file
% Peter Burns 5 Aug. 2002
% Copyright (c) International Imaging Industry Association

status = 0;
stuff = size(array);
nlin = stuff(1);
npix = stuff(2);
if size(stuff)==[1 2];
   ncol = 1;
else;
   ncol = stuff(3);
end;

temp = [oepath,oename];
oedat =load(temp);
%oedat = oename;
dimo = size(oedat);
if dimo(2) ~=ncol;
   status = 1;
   return;
end;
if ncol==1;
   for i=1: nlin;
      for j = 1: npix;
        array(i,j) = oedat( array(i,j)+1, ncol);
      end;
   end;
else;
   for i=1: nlin;
      for j = 1: npix;
         for k=1:ncol;
            array(i,j,k) = oedat( array(i,j,k)+1, k);
         end;
      end;
   end;
end;
