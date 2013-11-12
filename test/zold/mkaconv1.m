%
% NAME
%   mkaconv -- build an nchan x nchan AIRS convolution matrix
%
% SYNOPSIS
%   function mkaconv(sfile, afile);
%
% INPUTS
%   sfile  - HDF-format SRF filename
%   afile  - matlab output filename (optional)
% 
% OUTPUT is a Matlab data file with the variables
%   Amat   - AIRS convolution matrix
%   afrq   - sorted & trimmed AIRS channel set
%   aind   - index of afrq in regular AIRS channel set
% 
% NOTES
%   derived from mksconv
%
%   mkaconv trims too-close neighbors from the AIRS L1b channel and
%   tabulates AIRS SRFs in an nchan x nchan matrix Amat where rows
%   are SRFs and columns are AIRS channel frequencies.  The trimming
%   is needed for Amat to have a reasonable condition number.
%
%   AIRS channels have a naturally variable spacing.  The channel
%   trimming is done by taking a linear function that gives a lower
%   bound on acceptable channel spacing, as a function of frequency,
%   and dropping channels that violate the bound
%
%   the main application of the matrix is to deconvolve AIRS data.
%   If c is AIRs data at the afrq grid, then d = inv(Amat)*c is the
%   deconvolved channel data, still at the AIRS grid
%
%   uses h4sdread.m from /asl/matlab/h4tools to read the SRF file
%
% AUTHOR
%  H. Motteler, 28 Aug 2012
%

% function mkaconv(sfile, afile)

% just hard code some params, for now
nargin = 2;
sfile = '/asl/matlab/srftest/srftables_m140f_withfake_mar08.hdf';
afile = 'airs/SRFtest';

addpath /asl/matlab/h4tools

% default SRF source file
if nargin < 1
  sfile = '/asl/data/airs/srf/srftablesV10.hdf';
end

% default matlab output file
if nargin < 2
  [path,name,ext] = fileparts(sfile);
  afile = [name,'.mat'];
end

% read the srf data
[alist, fattr] = h4sdread(sfile);

for i = 1 : length(alist)
  switch alist{i}{1}
    case 'chanid', chanid = double(alist{i}{2})';
    case 'freq',   freq   = double(alist{i}{2})';
    case 'fwgrid', fwgrid = double(alist{i}{2})';
    case 'srfval', srfval = double(alist{i}{2})';
    case 'width',  width  = double(alist{i}{2})';
  end
end
clear alist fattr

[nchan, nspts] = size(srfval);

%-----------------------------------------------
% drop close neighbors from the L1b channel set
%-----------------------------------------------

% start with the L1B channel set
nchan = 2378;

% sort by channel frequency, fsrt = freq(isrt)
[fsrt, isrt] = sort(freq(1:nchan));

% set a reasonable lower bound on channel step size as a linear
% function a*x + b of frequency
a = 4e-4;
b = -0.04;

% drop channels that are too close to neighbors
ftmp = zeros(nchan,1);
itmp = zeros(nchan,1);

ftmp(1) = fsrt(1);
itmp(1) = isrt(1);
i = 1;
j = 2;

while 1

  while j <= nchan && fsrt(j) <= ftmp(i) + a*ftmp(i) + b
    j = j + 1;
  end

  % after the inner loop either j > nchan and we're done, 
  % or fsrt(j) > ftmp(i) + a*ftmp(i) + b, and we have a good 
  % next channel freq
  
  if j > nchan
    break
  end

  ftmp(i+1) = fsrt(j);
  itmp(i+1) = isrt(j);

  i = i + 1;
  j = j + 1;
end

% switch to the trimmed grid
nchan = i;
afrq = ftmp(1:nchan);
aind = itmp(1:nchan);

chanid = chanid(aind);
srfval = srfval(aind, :);
width = width(aind);

%-----------------------------------
% build an nchan x nchan SRF matrix
%-----------------------------------

% srfreq is the nchan x nspts frequency grids for srfval
srfreq = (width * fwgrid) + afrq * ones(1, nspts);

Amat = zeros(nchan, nchan);

for ci = 1 : nchan

   v1 = srfreq(ci, 1);
   v2 = srfreq(ci, nspts);

   i1 = min(find(v1 <= afrq));
   i2 = max(find(afrq <= v2));

%  stmp = interp1(srfreq(ci,:), srfval(ci,:), afrq(i1:i2), 'linear');
   stmp = interp1(srfreq(ci,:), srfval(ci,:), afrq(i1:i2), 'spline');

   Amat(ci, i1:i2) = stmp ./ sum(stmp);

%  plot(afrq(i1:i2), stmp)
%  pause

   % sanity check for valid array indices
   if i1 < i2
     continue
   end

   keyboard

end

% see how we did
% cond(Amat)

% save the results
eval(sprintf('save %s Amat afrq aind', afile));

