%
% NAME
%   mksconv1 -- build a sparse AIRS convolution matrix
%
% SYNOPSIS
%   [sconv, sfreq, tfreq] = mksconv1(sfile, cfreq, dvs)
%
% INPUTS
%   sfile  - AIRS HDF SRF tabulation file
%   cfreq  - desired AIRS channel frequencies
%   dvs    - sfreq step size, default is 0.0025
% 
% OUTPUTS
%   sconv  - m x n sparse convolution matrix
%   sfreq  - n-vector, sconv column frequency grid
%   tfreq  - m-vector, sconv row frequency grid
% 
% DESCRIPTION
%   mksconv1 takes a vector of AIRS channel frequencies and returns
%   an m x n sparse convolution matrix.  sfreq is chosen internally
%   to span the SRFs of the requested channels.  If r is radiances
%   at the sfreq grid, then c = sconv * r is AIRS channel radiances.
%
%   cfreq is desired channel frequencies and tfreq is nominal center
%   frequencies of the tabulated SRFs.  If these differ by more than
%   0.02 1/cm, a warning is given.  mksconv1 does not sort tfreq by
%   frequency.
%
% AUTHOR
%   H. Motteler, 20 June 2013
%

function [sconv, sfreq, tfreq] = mksconv1(sfile, cfreq, dvs)

% default to a 0.0025 1/cm input grid
if nargin < 3
  dvs = 0.0025;
end

% read the srf data
chanid = double(hdfread(sfile, 'chanid'));
tfreq  = double(hdfread(sfile, 'freq'));
fwgrid = double(hdfread(sfile, 'fwgrid'));
srfval = double(hdfread(sfile, 'srfval'));
width  = double(hdfread(sfile, 'width'));

% match channel sets
tfreq = tfreq(:);
cfreq = cfreq(:);
ix = interp1(tfreq, (1:length(tfreq))', cfreq, 'nearest');
vdif = max(abs(tfreq(ix) - cfreq));
if vdif > 0.02
  fprintf(1, 'mksconv1: warning, max(abs(tfreq - cfreq)) = %g\n', vdif)
end

% trim tabulated SRF data to common channels
tfreq = tfreq(ix);
srfval = srfval(ix, :);
width = width(ix, :);
[nchan, nspts] = size(srfval);

% tgrid is an nchan x nspts array of frequency points for srfval
tgrid = (width * fwgrid) + tfreq * ones(1, nspts);

% get spanning frequencies for this SRF set
v1 = min(tgrid(:));  
v2 = max(tgrid(:));  
v1 = ceil(v1 / dvs) * dvs;
v2 = floor(v2 / dvs) * dvs;

% get the sconv column grid
n = round((v2 - v1) / dvs);
sfreq = v1 + (0 : n)' * dvs;
nfreq = n + 1;
% isequal(v2, sfreq(end))

% initialize sparse matrix lists
si = []; sj = []; sd = [];

% loop on channels, interpolate SRFs to sfreq
for i = 1 : nchan

  t1 = tgrid(i,1);
  t2 = tgrid(i,nspts);
  jx = find(t1 <= sfreq & sfreq <= t2);
  k = length(jx);

  si = [si; ones(k, 1) * i];
  sj = [sj; jx];

% stmp = interp1(tgrid(i,:), srfval(i,:), sfreq(jx), 'spline');
  stmp = interp1(tgrid(i,:), srfval(i,:), sfreq(jx), 'linear');
  stmp = stmp ./ sum(stmp);

  sd = [sd; stmp];
end

% build a sparse matrix from the lists
sconv = sparse(si, sj, sd, nchan, nfreq, length(sd));

