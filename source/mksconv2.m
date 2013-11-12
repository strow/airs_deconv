%
% NAME
%   mksconv2 -- build a sparse AIRS convolution matrix
%
% SYNOPSIS
%   [sconv, sfreq, tfreq] = mksconv2(sfile, cfreq, dvs)
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
%   mksconv2 takes a vector of AIRS channel frequencies and returns
%   an m x n sparse convolution matrix.  sfreq is chosen internally
%   to span the SRFs of the requested channels.  if r is radiances
%   at the sfreq grid, then c = sconv * r is channel radiances.
%
%   cfreq is desired channels and tfreq corresponding frequencies of
%   the tabulated SRFs.  These are sorted and the closest pairs are
%   matched, if the difference is less than 0.04 1/cm.  tfreq may be
%   shorter than cfreq if there is no match for some channels.
%
%   mksconv2 calls seq_match and the 2008 version of h4sdread.
%
% AUTHOR
%   H. Motteler, 20 June 2013
%

function [sconv, sfreq, tfreq] = mksconv2(sfile, cfreq, dvs)

% default to a 0.0025 1/cm input grid
if nargin < 3
  dvs = 0.0025;
end

% read the srf data
[alist, fattr] = h4sdread(sfile);

for i = 1 : length(alist)
  switch alist{i}{1}
    case 'chanid', chanid = double(alist{i}{2})';
    case 'freq',   tfreq   = double(alist{i}{2})';
    case 'fwgrid', fwgrid = double(alist{i}{2})';
    case 'srfval', srfval = double(alist{i}{2})';
    case 'width',  width  = double(alist{i}{2})';
  end
end
clear alist fattr

% sort and match channel sets
cfreq = sort(cfreq(:));
[fsrt, isrt] = sort(tfreq);
[ci, fi] = seq_match(cfreq, fsrt, 0.04);
ix = isrt(fi);
% isequal(fsrt(fi), tfreq(ix))

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

