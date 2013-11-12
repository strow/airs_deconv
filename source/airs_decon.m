%
% NAME
%   airs_decon -- deconvolve AIRS radiances
%
% SYNOPSIS
%   [brad, bfrq] = airs_decon(arad, afrq, sfile, bfile, dvb)
%
% INPUTS
%   arad   - AIRS channel radiances, m x k array
%   afrq   - AIRS channel frequencies, m-vector
%   sfile  - AIRS HDF4 SRF tabulation file
%   bfile  - cache for deconvolution matrix
%   dvb    - deconvolution frequency step
%
% OUTPUTS
%   brad   - deconvolved radiances
%   bfrq   - deconvolution frequency grid
%
% DISCUSSION
%   airs_decon deconvolves AIRS radiances to an intermediate grid.  
%   This is done by taking the pseudo-inverse of the SRF convolution
%   matrix and applying this to the AIRS channel radiances.
%
%   Calculating the pseudo-inverse matrix for the deconvolution is
%   relatively slow so this data is saved in a file.  The matrix is
%   a function of the SRF tabulation file, channel frequencies, and
%   deconvolution grid step, and is updated if any of these change.
%
% AUTHOR
%  H. Motteler, 30 Oct 2013
%

function [brad, bfrq] = airs_decon(arad, afrq, sfile, bfile, dvb)

% default deconvolution matrix cache
if nargin < 4
  bfile = 'bconv.mat';
end

% default deconvolution frequency step
if nargin < 5
  dvb = 0.1;
end

% check if cache exists and is up-to-date
bok = 0;
if exist(bfile) == 2
  d = load(bfile, 'sfile', 'afrq', 'dvb');
  if strcmp(d.sfile, sfile) && isequal(d.afrq, afrq) && isequal(d.dvb, dvb)
    load(bfile)
    bok = 1;
  end
end

% if not, recalculate and save the pseudo-inverse
if ~bok
  fprintf(1, 'recalculating the pseudo-inverse...\n')
  [bconv, bfrq, fx] = mksconv2(sfile, afrq, dvb);
  binv = pinv(full(bconv));
  save(bfile, 'binv', 'bconv', 'bfrq', 'sfile', 'afrq', 'dvb')
end

% apply the deconvolution
brad = binv * arad;

