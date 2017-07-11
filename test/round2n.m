%
% NAME
%   round2n -- round to n significant digits
%
% SYNOPSIS
%   y = round2n(x, n)
%
% INPUTS
%   x  -  value to round
%   n  -  number of significant digits
%
% OUTPUT
%   y  -  x rounded to n significant digits
%

function y = round2n(x, n)

k = 10.^floor(log10(abs(x)) - n + 1);

y = round(x ./ k) .* k;

