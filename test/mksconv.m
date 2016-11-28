%
% NAME
%   mksconv -- build a sparse AIRS convolution matrix
%
% SYNOPSIS
%   [sconv, sfreq, tfreq] = mksconv(sfile, cind, dvs)
%
% INPUTS
%   sfile  - AIRS HDF SRF tabulation file
%   cind   - channel indices in SRF tabulation
%   dvs    - sfreq step size, default is 0.0025
% 
% OUTPUTS
%   sconv  - m x n sparse convolution matrix
%   sfreq  - n-vector, sconv column frequency grid
%   tfreq  - m-vector, sconv row frequency grid
% 
% DESCRIPTION
%   mksconv takes an m-vector of AIRS channel indices and returns 
%   an m x n sparse convolution matrix.  cind is channel indices and
%   tfreq corresponding frequencies of the tabulated SRFs.  sfreq is
%   chosen internally to span the SRFs of the requested channels.
%
%   if r is radiances at the sfreq grid, then c = sconv * r is
%   channel radiances at the tfreq grid
%
%   this version of mksconv works from channel indices.  mksconv1
%   works from a channel frequency list.  mksconv2 also works from 
%   a frequency list, and returns sconv and tfreq in sorted order
%
% AUTHOR
%   H. Motteler, 22 Nov 2016
%

function [sconv, sfreq, tfreq] = mksconv(sfile, cind, dvs)

% default to a 0.0025 1/cm input grid
if nargin < 3
  dvs = 0.0025;
end

% read the SRF tabulation
[tfreq, tgrid, srfval] = srf_read(sfile);
tfreq = tfreq(cind);
tgrid = tgrid(cind, :);
srfval = srfval(cind, :);
[nchan, nspts] = size(srfval);

% get spanning frequencies for this SRF set
v1 = min(tgrid(:));  
v2 = max(tgrid(:));  
v1 = ceil(v1 / dvs) * dvs;
v2 = floor(v2 / dvs) * dvs;

% get the sconv column grid
n = round((v2 - v1) / dvs);
sfreq = v1 + (0 : n)' * dvs;
nfreq = n + 1;
% isclose(v2, sfreq(end))

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

