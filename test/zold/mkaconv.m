%
% NAME
%   mkaconv -- build a sparse AIRS convolution matrix
%
% SYNOPSIS
%   [aconv, tfreq] = mkaconv(sfile, cind, sfreq)
%
% INPUTS
%   sfile  - AIRS HDF SRF tabulation file
%   cind   - channel indices in SRF tabulation
%   sfreq  - n-vector, sconv col frequency grid 
% 
% OUTPUTS
%   aconv  - m x n sparse convolution matrix
%   tfreq  - m-vector, sconv row frequency grid
% 
% DESCRIPTION
%   Similar to mksconv but takes the column frequency grid sfreq
%   rather than a size as an input parameter.  If r is radiance at
%   the sfreq grid, then c = aconv * r is channel radiances at the
%   tfreq grid.
%
% AUTHOR
%   H. Motteler, 21 Sep 2016
%

function [aconv, tfreq] = mkaconv(sfile, cind, sfreq)

% sconv columns
nfreq = length(sfreq);

% read the srf data
[tfreq, tgrid, srfval] = srf_read(sfile);
tfreq = tfreq(cind);
tgrid = tgrid(cind, :);
srfval = srfval(cind, :);
[nchan, nspts] = size(srfval);

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
aconv = sparse(si, sj, sd, nchan, nfreq, length(sd));

