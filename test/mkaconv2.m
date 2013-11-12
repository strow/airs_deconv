%
% NAME
%   mkaconv2 -- build an nchan x nchan AIRS convolution matrix
%
% SYNOPSIS
%   function mkaconv2(cfreq, sfile, afile);
%
% INPUTS
%   cfreq  - AIRS channel frequencies
%   sfile  - AIRS HDF-format SRF filename
%   afile  - matlab output filename (optional)
% 
% OUTPUT is a Matlab data file with the variables
%   Amat   - AIRS convolution matrix
%   afrq   - sorted trimmed AIRS channel set
%   aind   - index of afrq in cfreq
%   sind   - index of afrq in SRF tabulation
% 
% NOTES
%   the channel frequencies cfreq are sorted and matched with
%   tabulated SRF center frequencies,
%
%   channel neighbors that are too close are dropped and SRFs are
%   tabulated in an nchan x nchan matrix Amat, where rows are SRFs
%   and columns AIRS channel frequencies.  The trimming improves 
%   the condition number, especially for the 1b channel set.
%
%   AIRS channels have a variable spacing.  The channel trimming 
%   is done by taking a linear function that gives a lower bound 
%   on acceptable channel spacing, as a function of frequency, and
%   dropping channels that violate the bound
%
%   the main application of Amat is the deconvolution of AIRS data.
%   If c is AIRs data at the afrq grid, then d = inv(Amat)*c is the
%   deconvolved channel data, still at the AIRS grid
%
%   derived from mksconv, which reads an HDF SRF file and builds a
%   sparse matrix for applying AIRS convolutions to high res data
%
% AUTHOR
%  H. Motteler, 28 Aug 2012
%

function mkaconv2(cfreq, sfile, afile)

% default SRF source file
if nargin < 2
  sfile = '/asl/data/airs/srf/srftablesV10.hdf';
end

% default matlab output file
if nargin < 3
  [path,name,ext] = fileparts(sfile);
  afile = [name,'.mat'];
end

addpath /home/motteler/mot2008/hdf/h4tools

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

[nsrf, nspts] = size(srfval);

%-------------------------------------------------
% drop close neighbors from the input channel set
%-------------------------------------------------

%  cfreq  - input AIRS channel set
%  fsrt   - cfreq in guaranteed sorted order
%  isrt   - indices for fsrt in cfreq
%  iSRF   - indices for fsrt in tabulated SRFs
%  freq   - tabulated SRF channel center frequencies

[fsrt, isrt] = sort(cfreq);
nchan = length(fsrt);
iSRF = interp1(freq, 1:length(freq), fsrt, 'nearest');

% check that channel sets are a reasonable match
df = max(abs(fsrt - freq(iSRF)));
if df > 0.001
  fprintf(1, 'mkaconv: max cfreq vs SRF chan diff %g\n', df);
end

% set a lower bound (from plot_chans1.m) on channel step size 
% as a linear function a*x + b of frequency x
a = 4e-4;
b = -0.04;

% set up the channel trimming loop
ftmp = zeros(nchan,1);
jsrt = zeros(nchan,1);
jSRF = zeros(nchan,1);

ftmp(1) = fsrt(1);
jsrt(1) = isrt(1);
jSRF(1) = iSRF(1);
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
  jsrt(i+1) = isrt(j);
  jSRF(i+1) = iSRF(j);

  i = i + 1;
  j = j + 1;
end

% switch to the trimmed grid
nchan = i;
afrq = ftmp(1:nchan);   % final frequency grid
sind = jSRF(1:nchan);   % index of afrq in SRF tabulation
aind = jsrt(1:nchan);   % index of afrq in source 1c grid

chanid = chanid(sind);    % AIRS channel ID's
srfval = srfval(sind, :); % SRF values
width = width(sind);      % SRF widths

% index sanity check
df2 = max(abs(afrq - freq(sind)));
if df2 > df
  fprintf(1, 'mkaconv: frequency grid error, df2=%g\n', df2);
  keyboard
end

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

   % i1 should be less than i2, but check...
   if i1 >= i2
     fprintf(1, 'mkaconv: SRF index error\n');
     keyboard
   end

%  stmp = interp1(srfreq(ci,:), srfval(ci,:), afrq(i1:i2), 'linear');
   stmp = interp1(srfreq(ci,:), srfval(ci,:), afrq(i1:i2), 'spline');

   % normalize rows 
   Amat(ci, i1:i2) = stmp ./ sum(stmp);

%  plot(afrq(i1:i2), stmp)
%  pause

end

% see how we did
fprintf(1, 'mkaconv: cond(Amat) = %g\n', cond(Amat));

% save the results
eval(sprintf('save %s Amat afrq aind sind', afile));

