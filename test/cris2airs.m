%
% NAME
%   cris2airs -- translate CrIS high res to AIRS radiances
%
% SYNOPSIS
%   [arad, afrq] = ...
%       cris2airs(rLW, rMW, rSW, vLW, vMW, vSW, sfile, cfrq, opt1)
%
% INPUTS
%   rLW    - CrIS LW radiance, nLW x k array
%   rMW    - CrIS MW radiance, nMW x k array
%   rSW    - CrIS SW radiance, nSW x k array
%   rLW    - CrIS LW frequency, nLW vector
%   rMW    - CrIS MW frequency, nMW vector
%   rSW    - CrIS SW frequency, nSW vector
%   sfile  - AIRS HDF4 SRF tabulation file
%   cfrq   - desired AIRS channel frequencies
%   opt1   - optional input parameters
%
% opt1 fields
%   dvb     - deconvolution grid step (default 0.1 cm-1)
%   resmode - ccast resolution mode (default hires2)
%
% OUTPUTS
%   arad   - simulated AIRS radiances, n x k array
%   afrq   - AIRS channel frequencies, n-vector
%
% DISCUSSION
%   cris2aris translates Cris to AIRS radiances by deconvolving
%   CrIS to an intermediate grid and then reconvolving with the 
%   AIRS channel response functions
%
%   cfrq is desired channel frequencies and afrq is nominal center
%   frequencies of the tabulated SRFs.  If these differ by more than
%   0.02 1/cm, a warning is given.  iasi2airs does not sort afrq by
%   frequency.
%
% COPYRIGHT
%   Copyright 2013-2015, Atmospheric Spectroscopy Laboratory.  
%   This code is distributed under the terms of the GNU GPL v3.
%
% AUTHOR
%   H. Motteler, 12 Sep 2015
%

% test env
addpath /home/motteler/cris/iasi_decon
addpath /home/motteler/mot2008/hdf/h4tools
addpath /home/motteler/cris/ccast/source

f1 = '/asl/data/cris/ccast/sdr60_hr/2015/231/SDR_d20150819_t0821345.mat';
d1 = load(f1);
rLW = squeeze(d1.rLW(:, 5, 15, 21:30));
rMW = squeeze(d1.rMW(:, 5, 15, 21:30));
rSW = squeeze(d1.rSW(:, 5, 15, 21:30));
vLW = d1.vLW; vMW = d1.vMW; vSW = d1.vSW;
clear d1
sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
cfrq = load('freq2645.txt');
opt1 = struct;
opt1.resmode = 'hires2'; 
nargin = 8;

% [arad, afrq] = cris2airs(rLW, rMW, rSW, vLW, vMW, vSW, sfile, cfrq, opt1);

% defaults
dvb = 0.1;            % deconv grid step size
wlaser = 773.1301;    % nominal value is OK here
resmode = 'hires2';   % ccast resolution mode

% process input options
if nargin == 9
  if isfield(opt1, 'dvb'), dvb = opt1.dvb; end
  if isfield(opt1, 'wlaser'), wlaser = opt1.wlaser; end
  if isfield(opt1, 'resmode'), resmode = opt1.resmode; end
end

% check that array sizes match
vLW = vLW(:); vMW = vMW(:); vSW = vSW(:); 
[m, nobs] = size(rLW);
if m ~= length(vLW)
  error('rLW and vLW do not conform')
end
[m, nobs] = size(rMW);
if m ~= length(vMW)
  error('rMW and vMW do not conform')
end
[m, nobs] = size(rSW);
if m ~= length(vSW)
  error('rSW and vSW do not conform')
end

% get the AIRS convolution matrix
[sconv, sfrq, tfrq] = mksconv1(sfile, cfrq, dvb);

% get cris LW user grid
[instLW, userLW] = inst_params('LW', wlaser, opt1);

% interpolate radiance to intermediate grid
[rLWb, vLWb] = finterp(rLW, vLW, dvb);

% match intermediate grids
[ix, jx] = seq_match(vLWb, sfrq);

% get AIRS LW channel subset
kx = find(userLW.v1 <= tfrq & tfrq <= userLW.v2);

% convolve to AIRS radiances
rLWa = sconv(kx, jx) * rLWb(ix, :);
vLWa = tfrq(kx);

% get cris MW user grid
[instMW, userMW] = inst_params('MW', wlaser, opt1);

% interpolate radiance to intermediate grid
[rMWb, vMWb] = finterp(rMW, vMW, dvb);

% match intermediate grids
[ix, jx] = seq_match(vMWb, sfrq);

% get AIRS MW channel subset
kx = find(userMW.v1 <= tfrq & tfrq <= userMW.v2);

% convolve to AIRS radiances
rMWa = sconv(kx, jx) * rMWb(ix, :);
vMWa = tfrq(kx);

% get cris SW user grid
[instSW, userSW] = inst_params('SW', wlaser, opt1);

% interpolate radiance to intermediate grid
[rSWb, vSWb] = finterp(rSW, vSW, dvb);

% match intermediate grids
[ix, jx] = seq_match(vSWb, sfrq);

% get AIRS SW channel subset
kx = find(userSW.v1 <= tfrq & tfrq <= userSW.v2);

% convolve to AIRS radiances
rSWa = sconv(kx, jx) * rSWb(ix, :);
vSWa = tfrq(kx);

% concatenate bands
arad = [rLWa; rMWa; rSWa];
afrq = [vLWa; vMWa; vSWa];

plot(afrq, real(rad2bt(afrq, arad)));

