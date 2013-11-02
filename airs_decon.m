%
% NAME
%   airs_decon -- deconvolve AIRS radiances
%
% SYNOPSIS
%   [rad3, bfreq] = airs_decon(sfile, cfreq, rad2, bfile, dvb)
%
% INPUTS
%   sfile  - SRF tabulation file
%   rad2   - channel radiances
%   cfreq  - channel frequencies
%   bfile  - cache for deconvolution matrix
%   dvb    - deconvolution frequency step
%
% OUTPUTS
%   rad3   - deconvolved radiances
%   bfreq  - deconvolution frequency grid
%
% DISCUSSION
%    derived from cris_test1, calls mksconv2
%
% AUTHOR
%  H. Motteler, 30 Oct 2013
%

function [rad3, bfreq] = airs_decon(sfile, cfreq, rad2, bfile, dvb)

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
  d = load(bfile, 'sfile', 'cfreq', 'dvb');
  if strcmp(d.sfile, sfile) && isequal(d.cfreq, cfreq) && isequal(d.dvb, dvb)
    load(bfile)
    bok = 1;
  end
end

% if not, recalculate and save the pseudo-inverse
if ~bok
  fprintf(1, 'recalculating the pseudo-inverse...\n')
  [bconv, bfreq, fx] = mksconv2(sfile, cfreq, dvb);
  binv = pinv(full(bconv));
  save(bfile, 'binv', 'bconv', 'bfreq', 'sfile', 'cfreq', 'dvb')
end

% apply the deconvolution
rad3 = binv * rad2;

