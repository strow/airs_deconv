%
% sup_gauss -- higher-order Gaussian function
%
% SYNOPSIS
%   y = sup_gauss(x, b, w, p)
%
% INPUTS
%   x = input
%   b = function center
%   w = full-width half max (FWHM)
%   p = higher-order exponent (default 1)
%   
% OUTPUTS
%   y = weight, max 1
%
% DISCUSSION
%
%   basic gaussian (from Wikipedia):
%     y = exp(-(x - b)^2 / (2*c^2));
%     b = center
%     c = std
%
%   Inflection points are at x = b - c and x = b + c,
%   FWHM = 2*sqrt(2*ln(2))*c for a standard Gaussian.
%
%   higher-order gaussian (from Wikipedia):
%     y = exp(-((x - b)^2 / (2*c^2))^p);
%     b = center
%     c = nominal std (for p = 1)
%     p = higher-order exponent
%
%   Increasing p flattens the top and increases the side slope.  
%   We want to specify the function in terms of center, FWHM, and
%   exponent p.  This is done with a scaling factor q, derived by
%   solving exp(-((x - b)^2 / (2*c^2))^p) = 1/2 for x.
%     
% AUTHOR
%  H. Motteler, 16 May 2018
%

function y = sup_gauss(x, b, w, p)

% default to standard Gaussian
if nargin == 3
  p = 1;
end

% scaling factor for standard Gaussian
% q = 2 * sqrt(2*log(2))

% scaling factor for higher-order Gaussian
q = 2 * sqrt(2) * log(2)^(1/(2*p));

% express c (nominal std) in terms of FWHM
c = w / q;

% compute the higher order Gaussian
y = exp(-((x - b).^2 / (2*c^2)).^p);

