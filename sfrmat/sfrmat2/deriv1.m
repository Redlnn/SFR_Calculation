function  [b] = deriv1(a, nlin, npix, fil)
% [b] = deriv1(a, nlin, npix, fil)   First derivative of array
%  Computes first derivative via FIR (1xn) filter
%  Edge effects are suppressed and vector size is preserved
%  Filter is applied in the npix direction only
%   a   = (nlin, npix) data array
%   fil = array of filter coefficients, eg [[-0.5 0.5]
%   b   = output (nlin, npix) data array
% Peter Burns 5 Aug. 2002
% Copyright (c) International Imaging Industry Association

 b = zeros(nlin, npix);
 nn = length(fil);
 for i=1:nlin;
  temp = conv(fil, a(i,:));
  b(i, nn:npix) = temp(nn:npix);    %ignore edge effects, preserve size
  b(i, nn-1) = b(i, nn);
 end


