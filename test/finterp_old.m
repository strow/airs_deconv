%
% NAME
%   finterp - fourier interpolation 
%
% SYNOPSIS
%   function [rad2, frq2] = finterp(rad1, frq1, dv2, opt);
%
% INPUTS
%   rad1 - radiance array, in column order
%   frq1 - rad1 frequencies, a monotonic vector
%   dv2  - rad2 channel spacing
%   opt  - optional parameters 
%
% OUTPUTS
%   rad2 - output radiance vector
%   frq2 - rad2 frequencies
%
% OPTS FIELDS
%   dv1  - specifiy dv1, independent of frq1
%   roff - rolloff at band edges, in frq1 units    
%   tol  - tolerance for the rational approximation
%          of dv1/dv2; the default is 1e-6
%
% NOTES
%   rad1 and rad2 may be arrays, but frq1 and frq2 are vectors;
%   each column of rad1 has the frequency scale frq1, and each
%   column of rad2 has the frequency scale frq2
%
%   The input frq1 must be monotonic, and it can can contain gaps.
%   The output frq2 will approximately span min(frq) to max(frq1),
%   and will not contain any gaps.
%
%   In many cases frq1 is a regular grid, and dv1 is frq1(2)-frq1(1),
%   but that is not required.  The rad1 data is interpolated to a dv1
%   spacing before further processing.  This is useful for cases where
%     1. frq1 been subject to roundoff or other error and dv1 is known,
%     2. frq1 is sparse or otherwise not regular,
%     3. we want to interpolate, as in taking kcarta 0.0025 1/cm data 
%        to a vlaser spaced grid.
%   The output grid frq2 will approximately span min(frq1) to max(frq1)
%   and will not contain any gaps.
%
%   tol should not be set to less than the uncertainty of dv1/dv2
%
% AUTHOR
%   H. Motteler, 15 Sept 99

function [rad2, frq2] = finterp(rad1, frq1, dv2, opts);

% defaults
dv1  = frq1(2) - frq1(1);  % get dv1 from frq1 spacing
roff = 0;		   % band edge rolloff in frq1 units
tol  = 1e-6;		   % uncertainty of dv1/dv2
info = 0;		   % flag to print FFT sizes

% option to override defaults with gopts fields
if nargin == 4
  optvar = fieldnames(opts);
  for i = 1 : length(optvar)
    vname = optvar{i};
    if exist(vname, 'var')
%     fprintf(1, 'finterp: setting %s\n', vname)
      eval(sprintf('%s = opts.%s;', vname, vname));
    else
      fprintf(1, 'finterp: unknown option variable %s\n', vname);
    end
  end
end

% guarantee frq1 is a column vector
frq1 = frq1(:);

% check that input args conform
[nrow,ncol] = size(rad1);
if nrow ~= length(frq1)
  error('rad1 rows must match frq1 length')
end

% frequency domain filter at band edges
if roff ~= 0
  rpts = round(roff / dv1);  % rolloff points
  f0 = (1+cos(pi+(0:rpts-1)*pi/(rpts)))/2;
  filt = [f0, ones(1,nrow-2*rpts), fliplr(f0)]';
  for i = 1 : ncol
    rad1(:,i) = rad1(:,i) .* filt;
  end
end

% shift v1 and v2 to exact multiples of dv
v1 = min(frq1);
v1ind = ceil(v1/dv1) + 1;
v1 = (v1ind - 1) * dv1;
v2 = max(frq1);
v2ind = floor(v2/dv1) + 1;
v2 = (v2ind - 1) * dv1;

dv1grid = (v1ind - 1 : v2ind - 1) * dv1;

% try to find vmax >= v2 whose indices n1 and n1 have a 
% plausible factorization, so the fft is reasonably fast

% express dv1/dv2 as a ratio of integers
[a1, a2] = rat(dv1/dv2, tol); 

% find n1 and n2 as multiples of a2 and a1
j = max(nextpow2(v2/(a2*dv1)), 0);
n1 = a2 * 2^j + 1;
n2 = a1 * 2^j + 1;
vmax = (n1 - 1) * dv1;

if info
  fprintf(1, 'old finterp: n1 = %7d, n2 = %5d, dx = %6.3e, dv = %6.3e\n', ...
              n1-1, n2-1, 1/(2*vmax), dv1);
end

% output frequency scale
v1ind2 = round(v1/dv2) + 1;
v2ind2 = round(v2/dv2) + 1;
frq2 = (v1ind2 - 1 : v2ind2 - 1) * dv2;

% initialize output array
rad2 = zeros(length(frq2), ncol);

% loop on radiance columns
for i = 1 : ncol

  % embed input data in a [0-vmax] interval
  warning off
  rad1a = interp1(frq1, rad1(:,i), dv1grid, 'linear');
  warning on
  rad1b = zeros(n1, 1);
  rad1b(v1ind:v2ind) = rad1a;

  % transform spectral data to an interferogram
  intf1 = ifft([rad1b; flipud(rad1b(2:n1-1,1))]);

  % extend (with zero fill) if n1 < n2, or truncate if n1 > n2
  intf2 = zeros(n2, 1);
  intf2(1:min(n1,n2)) = intf1(1:min(n1,n2));

  % transform the new interferogram back to radiances
  rad2a = fft([intf2; flipud(intf2(2:n2-1,1))]);

  % select radiances from the input band
  rad2(:,i) = rad2a(v1ind2:v2ind2);

end

