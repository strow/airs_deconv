%
% NAME
%   airs_decon -- deconvolve AIRS radiances
%
% SYNOPSIS
%   [brad, bfrq] = airs_decon(arad, afrq, sfile, dvb)
%
% INPUTS
%   arad   - AIRS channel radiances, m x k array
%   afrq   - AIRS channel frequencies, m-vector
%   sfile  - AIRS HDF4 SRF tabulation file
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
%   Calculating the pseudo-inverse is slow, so this data is saved 
%   in a persistent struct.  The pseudo-inverse is a function of the
%   SRF tabulation file, channel frequencies, and deconvolution grid
%   step, and is updated if any of these change.
%
% AUTHOR
%  H. Motteler, 30 Oct 2013
%

function [brad, bfrq] = airs_decon(arad, afrq, sfile, dvb)

% default deconvolution frequency step
if nargin < 4
  dvb = 0.1;
end

% persistent struct to cache the pseudo-inverse matrix
persistent bc

% check that cache exists and is valid for current params
bcOK = ~isempty(bc) && ...
         isfield(bc, 'sfile') && strcmp(bc.sfile, sfile) && ...
         isfield(bc, 'afrq') && isequal(bc.afrq, afrq) && ...
         isfield(bc, 'dvb') && isequal(bc.dvb, dvb);

% if cache is empty or invalid, recreate the pseudo-inverse
if ~bcOK
  bc = struct;
  fprintf(1, 'calculating the pseudo-inverse...\n')
  [bconv, bfrq, fx] = mksconv2(sfile, afrq, dvb);
  bc.binv = pinv(full(bconv));
  bc.sfile = sfile;
  bc.afrq = afrq;
  bc.bfrq = bfrq;
  bc.dvb = dvb;
end

% apply the cached deconvolution
brad = bc.binv * arad;

% return cached frequency grid
bfrq = bc.bfrq;

