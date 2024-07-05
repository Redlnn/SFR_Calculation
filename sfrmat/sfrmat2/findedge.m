function  [slope, int] = findedge(cent, nlin)
% [slope, int] = findedge(cent, nlin)  Fits linear equation to data
% Fit linear equation to data, written to process edge location array
%   cent = array of (centroid) values
%   nlin = length of cent
%   slope and int are from the least-square fit
%    x = int + slope*cent(x)
%  Note that this is the inverse of the usual cent(x) = int + slope*x
%  form
% Peter Burns 5 Aug. 2002
% Copyright (c) International Imaging Industry Association

 index=[0:nlin-1];
 [slope int] = polyfit(index, cent, 1);            % x = f(y)
return
