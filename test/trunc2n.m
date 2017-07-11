%
% NAME
%   trunc2n -- truncate to n significant digits
%
% SYNOPSIS
%   y = trunc2n(x, n)
%
% INPUTS
%   x  -  value to truncate
%   n  -  number of significant digits
%
% OUTPUT
%   y  -  x truncated to n significant digits
%

function y = trunc2n(x, n)

k = 10.^floor(log10(abs(x)) - n + 1);

y = floor(x ./ k) .* k;

