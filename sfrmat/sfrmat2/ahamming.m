function [data] = ahamming(n, mid)
% [data] = ahamming(n, mid)  Generates asymmetrical Hamming window
%  array. If mid = (n+1)/2 then the usual symmetrical Hamming array
%  is returned
%   n = length of array
%   mid = midpoint (maximum) of window function
%   data = window array (nx1)
% Peter Burns 5 Aug. 2002
% Copyright (c) International Imaging Industry Association

data = zeros(n,1);

wid1 = mid-1;
wid2 = n-mid;
wid = max(wid1, wid2);
for i = 1:n;
	arg = i-mid;
	data(i) = 0.54 + 0.46*cos( pi*arg/wid );
end;
