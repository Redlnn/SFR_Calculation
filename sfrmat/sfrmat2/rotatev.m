function [a, nlin, npix, rflag] = rotatev(a)
% [a, nlin, npix, rflag] = rotatev(a)    Rotate array
%
% Rotate array so that long dimensions is vertical (line) drection
% a = input array(npix, nlin, ncol)
% nlin, npix are after rotation if any
% flag = 0 no roation, = 1 rotation was performed
%
% Peter Burns 5 Aug. 2002
% Copyright (c) International Imaging Industry Association

dim = size(a);
nlin = dim(1);
npix = dim(2);
if size(dim)==[1 2];
  ncol =1;
 else;
  ncol = dim(3);
 end;

 rflag =0;
if npix>nlin;
 rflag =1;
 b    = zeros(npix, nlin, ncol);
 temp = zeros(npix, nlin);
                    
 for i=1:ncol;
  temp = a(:, :, i)';
  b(:,:,i) = temp;
 end;
 a = b;
  temp=nlin;
  nlin = npix;
  npix = temp;
end;
